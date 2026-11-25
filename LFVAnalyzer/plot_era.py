#!/usr/bin/env python3
"""
Automate per-era plotIt execution by dynamically selecting histogram include files
based on workflow stage (process, TauFF calculation/application, DNN, yields) without
manually editing config_<era>.yml.
"""
# Usage examples:
#   python plot_era.py -I process_1111_v1/muon/v15_2024
#   python plot_era.py -f -I process_1111_v1 -C muon -Y v15_2024
#   python plot_era.py -f -I process_1111_v1 -C muon
#   python plot_era.py -f -I process_1111_v1_fake_lss/muon/v15_2024 -C muon
#   python plot_era.py -I process_1111_v1_FF/muon/v15_2024 -F
#   python plot_era.py -I dnn_out/muon/v15_2024_postprocess -D

import argparse
import os
import re
import shlex
import shutil
import subprocess
import sys
import tempfile

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(SCRIPT_DIR)
PLOTIT_ROOT = os.path.join(REPO_ROOT, 'plotIt')
PLOTIT_BIN = os.path.join(PLOTIT_ROOT, 'plotIt')

ERAS = {
    '2022':     {'dir_candidates': ['v12_2022', 'v2022', '2022']},
    '2022EE':   {'dir_candidates': ['v12_2022EE', 'v2022EE', '2022EE']},
    '2023':     {'dir_candidates': ['v12_2023', 'v2023', '2023']},
    '2023BPix': {'dir_candidates': ['v12_2023BPix', 'v12_2023_BPix', 'v2023_BPix', 'v2023BPix', '2023BPix', '2023_BPix']},
    '2024':     {'dir_candidates': ['v15_2024', 'v2024', '2024']},
}

ERA_ALIASES = {
    '2022': '2022',
    'v12_2022': '2022',
    'v2022': '2022',
    '2022ee': '2022EE',
    '2022_ee': '2022EE',
    'v12_2022ee': '2022EE',
    'v2022ee': '2022EE',
    '2023': '2023',
    'v12_2023': '2023',
    'v2023': '2023',
    '2023bpix': '2023BPix',
    '2023_bpix': '2023BPix',
    'v12_2023bpix': '2023BPix',
    'v12_2023_bpix': '2023BPix',
    'v2023bpix': '2023BPix',
    'v2023_bpix': '2023BPix',
    '2024': '2024',
    'v15_2024': '2024',
    'v2024': '2024',
}

FAKE_REGIONS = ['lss', 'los', 'tss', 'tos']
FAKE_REGION_PATH_RE = re.compile(r'(?:^|/)[^/]+_fake_(lss|los|tss|tos)(?:/|$)')
ERA_TOKEN_RE = re.compile(r'(2023_?BPix|2022_?EE|2024|2023|2022)', re.IGNORECASE)
SIGNAL_TYPE_RE = re.compile(r'^\s+type\s*:\s*[\'"]?signal[\'"]?(?:\s*(?:#.*)?)?$')


def fail(msg, code=1):
    print(f"[ERROR] {msg}")
    sys.exit(code)


def resolve_channel(explicit_channel, abs_input):
    if explicit_channel:
        return explicit_channel
    lower = abs_input.lower()
    has_mu = 'muon' in lower
    has_el = 'electron' in lower
    if has_mu and has_el:
        fail(f"Ambiguous channel in input path '{abs_input}' (contains both 'muon' and 'electron'); specify -C/--channel.")
    if has_mu:
        return 'muon'
    if has_el:
        return 'electron'
    fail(f"Could not auto-detect channel ('muon' or 'electron') from '{abs_input}'; specify -C/--channel.")


def resolve_era(explicit_era, abs_input):
    if explicit_era:
        key = explicit_era.strip().lower()
        if key not in ERA_ALIASES:
            fail(f"Unsupported era '{explicit_era}'. Supported: 2022, 2022EE, 2023, 2023BPix, 2024 (with optional v12_/v15_ prefix).")
        return ERA_ALIASES[key]

    # Prefer basename matches first so parent directory names cannot confuse detection
    for candidate_str in (os.path.basename(abs_input), abs_input):
        matches = ERA_TOKEN_RE.findall(candidate_str)
        canonical = {ERA_ALIASES[m.lower()] for m in matches if m.lower() in ERA_ALIASES}
        if len(canonical) == 1:
            return next(iter(canonical))
        if len(canonical) > 1:
            fail(f"Ambiguous era tokens {sorted(canonical)} in '{candidate_str}'; specify -Y/--year explicitly.")

    fail(f"Could not auto-detect era from '{abs_input}'; specify -Y/--year explicitly.")


def resolve_era_subdir(channel_dir, era, postfix):
    """Find <cand>_postprocess{postfix} first, then fall back to raw <cand>{postfix} / <cand>."""
    candidates = ERAS[era]['dir_candidates']
    for cand in candidates:
        post_dir = os.path.join(channel_dir, f'{cand}_postprocess{postfix}')
        if os.path.isdir(post_dir):
            return post_dir
    for cand in candidates:
        if postfix:
            raw_post_dir = os.path.join(channel_dir, f'{cand}{postfix}')
            if os.path.isdir(raw_post_dir):
                return raw_post_dir
        raw_dir = os.path.join(channel_dir, cand)
        if os.path.isdir(raw_dir):
            return raw_dir
    return None


def resolve_fake_targets(raw_input, channel, explicit_year, regions, postfix):
    """Expand -I <version> in -f mode into concrete (region, era, input_dir, output_dir) targets."""
    if not raw_input.strip() or raw_input.strip().rstrip('/') in ('', '.'):
        fail("In -f (--calcff) version mode, provide a version prefix with -I (e.g. -I process_1111_v1) or a full fake region directory path.")

    version_prefix = os.path.abspath(raw_input.rstrip('/'))
    region_chan_dirs = {}
    for region in regions:
        chan_dir = os.path.join(f'{version_prefix}_fake_{region}', channel)
        if not os.path.isdir(chan_dir):
            fail(f"Expected fake channel directory does not exist: '{chan_dir}'")
        region_chan_dirs[region] = chan_dir

    targets = []
    if explicit_year:
        era = resolve_era(explicit_year, '')
        for region in regions:
            chan_dir = region_chan_dirs[region]
            era_dir = resolve_era_subdir(chan_dir, era, postfix)
            if not era_dir:
                fail(f"Missing era directory for era='{era}' under '{chan_dir}' (searched *_postprocess{postfix} and raw era dirs).")
            out_dir = os.path.join(chan_dir, f'figure_{era}{postfix}')
            targets.append((region, era, era_dir, out_dir))
        return targets

    for era in ERAS:
        per_region = {}
        missing_regions = []
        for region in regions:
            chan_dir = region_chan_dirs[region]
            era_dir = resolve_era_subdir(chan_dir, era, postfix)
            if era_dir:
                per_region[region] = (era_dir, os.path.join(chan_dir, f'figure_{era}{postfix}'))
            else:
                missing_regions.append(region)

        if not missing_regions:
            for region in regions:
                era_dir, out_dir = per_region[region]
                targets.append((region, era, era_dir, out_dir))
        elif per_region:
            found_regions = [r for r in regions if r in per_region]
            fail(f"Incomplete fake region directories for era '{era}' ({channel}): found in {found_regions}, missing in {missing_regions}.")

    if not targets:
        fail(f"No Run 3 era directories found for channel '{channel}' under '{version_prefix}_fake_*'.")
    return targets


def resolve_plots_include(args):
    if args.plots:
        plots = [p.strip() for p in args.plots.split(',') if p.strip()]
        if not plots:
            fail("--plots was specified but contained no valid histogram YAML filenames.")
        return plots

    if args.tauFF:
        return ['histos_yield_S5.yml'] if args.yield_only else ['histos_FFapply.yml', 'histos_yield_S5.yml']
    if args.calcTauFF:
        if args.yield_only:
            fail("-y (--yield-only) cannot be combined with -f (--calcff); histos_FFcalc.yml has no separate yield-only YAML.")
        return ['histos_FFcalc.yml']
    if args.DNN:
        if args.yield_only:
            fail("-y (--yield-only) cannot be combined with -D (--dnn); histos_dnn.yml has no separate yield-only YAML.")
        return ['histos_dnn.yml']
    return ['histos_yield.yml'] if args.yield_only else ['histos_yield.yml', 'histos.yml', 'histos_reco.yml']


def extract_files_include(base_cfg_path):
    """Extract active YAML filenames from the top-level 'files: include: [...]' block."""
    with open(base_cfg_path, 'r') as f:
        lines = f.readlines()

    in_files = False
    found_files_section = False
    included_files = []

    for line in lines:
        stripped = line.strip()
        if in_files:
            if not stripped or stripped.startswith('#'):
                continue
            if line[0].isspace():
                m = re.match(r'^include\s*:\s*\[(.*)\]\s*(?:#.*)?$', stripped)
                if m:
                    for item in m.group(1).split(','):
                        cleaned = item.strip().strip('\'"')
                        if cleaned:
                            included_files.append(cleaned)
                continue
            in_files = False

        if not line[:1].isspace() and re.match(r'^files\s*:', line):
            if found_files_section:
                fail(f"Multiple top-level 'files:' sections found in '{base_cfg_path}'.")
            in_files = True
            found_files_section = True

    if not found_files_section:
        fail(f"No top-level 'files:' section found in '{base_cfg_path}'.")
    if not included_files:
        fail(f"No active 'include: [...]' found under 'files:' in '{base_cfg_path}'.")
    return included_files


def strip_signal_blocks(files_yml_path):
    """Remove all active top-level entries with 'type: signal' from a files*.yml file."""
    with open(files_yml_path, 'r') as f:
        lines = f.readlines()

    blocks = []
    current_block = None

    for lineno, raw_line in enumerate(lines, 1):
        stripped = raw_line.strip()
        if not stripped or stripped.startswith('#'):
            continue
        if not raw_line[0].isspace():
            if current_block is not None:
                blocks.append(current_block)
            current_block = [raw_line]
        else:
            if current_block is None:
                fail(f"Unexpected indented YAML line before any top-level entry in '{files_yml_path}' line {lineno}: {raw_line.rstrip()}")
            current_block.append(raw_line)

    if current_block is not None:
        blocks.append(current_block)

    kept_blocks = []
    removed_count = 0
    for block in blocks:
        if any(SIGNAL_TYPE_RE.match(line) for line in block[1:]):
            removed_count += 1
        else:
            block_text = ''.join(block)
            kept_blocks.append(block_text if block_text.endswith('\n') else block_text + '\n')

    if not kept_blocks:
        fail(f"No non-signal file entries remain after stripping signals in '{files_yml_path}'.")

    return '\n'.join(kept_blocks), removed_count


def _format_include_block(section_name, includes):
    include_literal = ', '.join(f"'{p}'" for p in includes)
    return f"{section_name}:\n  include: [{include_literal}]\n\n"


def _replace_top_level_section(lines, section_name, replacement_block, base_cfg_path):
    new_lines = []
    in_section = False
    replaced = False
    section_re = re.compile(rf'^{re.escape(section_name)}\s*:')

    for line in lines:
        if in_section:
            stripped = line.strip()
            if not stripped or line[0].isspace() or stripped.startswith('#'):
                continue
            in_section = False

        if not line[:1].isspace() and section_re.match(line):
            if replaced:
                fail(f"Multiple top-level '{section_name}:' sections found in '{base_cfg_path}'.")
            in_section = True
            replaced = True
            new_lines.append(replacement_block)
        else:
            new_lines.append(line)

    if not replaced:
        fail(f"No top-level '{section_name}:' section found in '{base_cfg_path}'.")

    return new_lines


def build_temp_config(base_cfg_path, plots_include, files_include=None, strip_systematics=False):
    with open(base_cfg_path, 'r') as f:
        lines = f.readlines()

    lines = _replace_top_level_section(lines, 'plots', _format_include_block('plots', plots_include), base_cfg_path)
    if files_include is not None:
        lines = _replace_top_level_section(lines, 'files', _format_include_block('files', files_include), base_cfg_path)
    if strip_systematics:
        lines = _replace_top_level_section(lines, 'systematics', "systematics: []\n\n", base_cfg_path)
    return ''.join(lines)


def find_split_data_files(input_dir, entries, channel):
    """Collect era-split data ROOT filenames (excluding merged files and syst variations)."""
    split_names = []
    for name in entries:
        if not name.endswith('.root') or name.startswith('.') or '__' in name:
            continue
        if not os.path.isfile(os.path.join(input_dir, name)):
            continue
        stem = name[:-5]
        if channel == 'muon':
            if name in ('hist_Muon.root', 'hist_SingleMuon.root'):
                continue
            if (stem.startswith('hist_Muon') and len(stem) > len('hist_Muon')) or (
                stem.startswith('hist_SingleMuon') and len(stem) > len('hist_SingleMuon')
            ):
                split_names.append(name)
        else:
            if name.lower() == 'hist_egamma.root':
                continue
            if stem.lower().startswith('hist_egamma') and len(stem) > len('hist_egamma'):
                split_names.append(name)
    return split_names


def ensure_merged_data(input_dir, channel, dry_run):
    """Ensure hist_Muon.root / hist_Egamma.root exists, auto-merging era-split data files via hadd if missing."""
    target_name = 'hist_Muon.root' if channel == 'muon' else 'hist_Egamma.root'
    target_path = os.path.join(input_dir, target_name)
    entries = sorted(os.listdir(input_dir))

    # Exact case check via os.listdir() so case-insensitive filesystems (e.g. macOS) are handled accurately
    if target_name in entries and os.path.isfile(target_path):
        return

    split_names = find_split_data_files(input_dir, entries, channel)
    if not split_names:
        fail(f"Merged data file '{target_path}' is missing and no era-split data ROOT files were found to hadd in '{input_dir}'.")

    split_paths = [os.path.join(input_dir, n) for n in split_names]
    if dry_run:
        preview = ', '.join(split_names[:3]) + (' ...' if len(split_names) > 3 else '')
        print(f"[DRY-RUN] Missing '{target_name}' in '{input_dir}' -> would run: hadd -f {target_path} <{len(split_names)} files: {preview}>")
        return

    hadd_bin = shutil.which('hadd')
    if not hadd_bin:
        fail(f"Merged data file '{target_path}' is missing (found {len(split_names)} split data files), but 'hadd' was not found in PATH.")

    print(f"[INFO] Missing '{target_name}' in '{input_dir}'; running hadd over {len(split_names)} split data files...")
    tmp_target = os.path.join(input_dir, f'.tmp_{os.getpid()}_{target_name}')
    try:
        proc = subprocess.run([hadd_bin, '-f', tmp_target] + split_paths, check=False)
        if proc.returncode != 0 or not os.path.isfile(tmp_target):
            fail(f"hadd failed (exit code {proc.returncode}) while creating '{target_path}'.", code=proc.returncode or 1)
        os.replace(tmp_target, target_path)
        print(f"[INFO] Created merged data file: {target_path}")
    finally:
        if os.path.exists(tmp_target):
            try:
                os.remove(tmp_target)
            except OSError:
                pass


def build_plotit_flags(args):
    flags = []
    if not args.no_yields:
        flags.append('-y')
    if not args.no_syst:
        flags.append('-s')
    if not args.no_qcd:
        flags.append('-q')
    if args.plotit_args:
        try:
            flags.extend(shlex.split(args.plotit_args))
        except ValueError as exc:
            fail(f"Invalid --plotit-args syntax: {exc}")
    return flags


def run_single(input_dir, output_dir, config_path, base_cfg_path, channel, era, plots_include, plotit_flags, args, exclude_signals=False, strip_systematics=False):
    if not os.path.isdir(input_dir):
        fail(f"Input directory does not exist: '{input_dir}'")

    ensure_merged_data(input_dir, channel, args.dry_run)

    orig_files = []
    filtered_files_content = ''
    removed_signals = 0
    if exclude_signals:
        orig_files = extract_files_include(base_cfg_path)
        chunks = []
        for fname in orig_files:
            fpath = os.path.join(config_path, fname)
            if not os.path.isfile(fpath):
                fail(f"Included files YAML '{fname}' (from '{base_cfg_path}') not found in '{config_path}'")
            chunk, n_removed = strip_signal_blocks(fpath)
            chunks.append(chunk)
            removed_signals += n_removed
        filtered_files_content = '\n'.join(chunks)

    if args.dry_run:
        dry_files_include = [f'.tmp_plot_era_files_{era}.yml'] if exclude_signals else None
        cfg_content = build_temp_config(
            base_cfg_path, plots_include, files_include=dry_files_include, strip_systematics=strip_systematics
        )
        cmd = [PLOTIT_BIN, '-i', input_dir, '-o', output_dir, os.path.join(config_path, f'.tmp_plot_era_{era}.yml')] + plotit_flags
        print(f"[DRY-RUN] Config base : {base_cfg_path}")
        print(f"[DRY-RUN] Input (-i)  : {input_dir}")
        print(f"[DRY-RUN] Output (-o) : {output_dir}")
        if exclude_signals:
            print(f"[DRY-RUN] files       : signal-excluded ({removed_signals} removed from {orig_files})")
        if strip_systematics:
            print("[DRY-RUN] systematics : cleared (luminosity-error preserved)")
        print(f"[DRY-RUN] plots: include: {plots_include}")
        print(f"[DRY-RUN] Command     : {' '.join(shlex.quote(c) for c in cmd)}")
        if args.keep_tmp:
            kept_files_include = None
            if exclude_signals:
                with tempfile.NamedTemporaryFile(
                    mode='w',
                    prefix=f'.tmp_plot_era_files_{era}_{os.getpid()}_',
                    suffix='.yml',
                    dir=config_path,
                    delete=False,
                ) as tmp_files_f:
                    tmp_files_f.write(filtered_files_content)
                kept_files_include = [os.path.basename(tmp_files_f.name)]
                print(f"[DRY-RUN] Kept temp files : {tmp_files_f.name}")
            cfg_content = build_temp_config(
                base_cfg_path, plots_include, files_include=kept_files_include, strip_systematics=strip_systematics
            )
            with tempfile.NamedTemporaryFile(
                mode='w',
                prefix=f'.tmp_plot_era_{era}_{os.getpid()}_',
                suffix='.yml',
                dir=config_path,
                delete=False,
            ) as tmp_f:
                tmp_f.write(cfg_content)
            print(f"[DRY-RUN] Kept temp config: {tmp_f.name}")
        return

    os.makedirs(output_dir, exist_ok=True)
    tmp_files_path = None
    tmp_cfg_path = None
    try:
        files_include = None
        if exclude_signals:
            with tempfile.NamedTemporaryFile(
                mode='w',
                prefix=f'.tmp_plot_era_files_{era}_{os.getpid()}_',
                suffix='.yml',
                dir=config_path,
                delete=False,
            ) as tmp_files_f:
                tmp_files_f.write(filtered_files_content)
                tmp_files_path = tmp_files_f.name
            files_include = [os.path.basename(tmp_files_path)]

        cfg_content = build_temp_config(
            base_cfg_path, plots_include, files_include=files_include, strip_systematics=strip_systematics
        )
        with tempfile.NamedTemporaryFile(
            mode='w',
            prefix=f'.tmp_plot_era_{era}_{os.getpid()}_',
            suffix='.yml',
            dir=config_path,
            delete=False,
        ) as tmp_f:
            tmp_f.write(cfg_content)
            tmp_cfg_path = tmp_f.name

        cmd = [PLOTIT_BIN, '-i', input_dir, '-o', output_dir, tmp_cfg_path] + plotit_flags
        print(f"[INFO] Running: {' '.join(shlex.quote(c) for c in cmd)}")
        proc = subprocess.run(cmd, check=False)
        if proc.returncode != 0:
            fail(f"plotIt failed (exit code {proc.returncode}) for era='{era}', input='{input_dir}'.", code=proc.returncode)
    finally:
        if not args.keep_tmp:
            for tmp_p in (tmp_cfg_path, tmp_files_path):
                if tmp_p and os.path.exists(tmp_p):
                    os.remove(tmp_p)


def parse_args():
    parser = argparse.ArgumentParser(description="Run plotIt for a single Run 3 era without editing config files")
    parser.add_argument('-I', '--input', dest='input', type=str, default='.', help="Input histogram directory (or version prefix with -f, e.g. process_1111_v1)")
    parser.add_argument('-O', '--output', dest='output', type=str, default='', help="Output figure directory (default: <parent>/figure_<input_basename><postfix>, or figure_<bare_era><postfix> in -f mode)")
    parser.add_argument('-Y', '--year', '--era', dest='year', type=str, default='', help="Single era: 2022, 2022EE, 2023, 2023BPix, or 2024 (auto-detected from input, or all 5 eras when omitted in -f version mode)")
    parser.add_argument('-C', '-ch', '--channel', dest='channel', type=str, default='', choices=['', 'electron', 'muon'], help="Channel: electron or muon (required with -f; auto-detected from input otherwise)")
    parser.add_argument('-D', '--dnn', dest='DNN', action='store_true', default=False, help="Run for DNN histograms (histos_dnn.yml)")
    parser.add_argument('-F', '--ff', dest='tauFF', action='store_true', default=False, help="Run for TauFF applied histograms (histos_FFapply.yml, histos_yield_S5.yml)")
    parser.add_argument('-f', '--calcff', dest='calcTauFF', action='store_true', default=False, help="Run for TauFF calculation histograms (histos_FFcalc.yml, excluding signals & systematics)")
    parser.add_argument('--regions', dest='regions', nargs='+', choices=FAKE_REGIONS, default=None, help="Subset of fake regions in -f version mode (default: lss los tss tos)")
    parser.add_argument('-y', '--yield-only', dest='yield_only', action='store_true', default=False, help="Use yield-only histogram list (histos_yield.yml or histos_yield_S5.yml with -F)")
    parser.add_argument('--plots', dest='plots', type=str, default='', help="Comma-separated histogram YAML files overriding mode selection")
    parser.add_argument('--no-yields', dest='no_yields', action='store_true', default=False, help="Omit -y when calling plotIt")
    parser.add_argument('--no-syst', dest='no_syst', action='store_true', default=False, help="Omit -s when calling plotIt")
    parser.add_argument('--no-qcd', dest='no_qcd', action='store_true', default=False, help="Omit -q when calling plotIt")
    parser.add_argument('--plotit-args', dest='plotit_args', type=str, default='', help="Additional flags passed verbatim to plotIt")
    parser.add_argument('--postfix', dest='postfix', type=str, default='', help="Postfix appended to default output directory name (and era subdir lookup in -f mode)")
    parser.add_argument('--dry-run', dest='dry_run', action='store_true', default=False, help="Validate config and print plotIt command without executing")
    parser.add_argument('--keep-tmp', dest='keep_tmp', action='store_true', default=False, help="Keep generated temporary YAML config inside plotIt config directory")
    return parser.parse_args()


def _validate_plotit_env(config_path, eras, plots_include, dry_run):
    if not os.path.isdir(config_path):
        fail(f"plotIt config directory does not exist: '{config_path}'")

    base_cfgs = {}
    for era in eras:
        base_cfg_path = os.path.join(config_path, f'config_{era}.yml')
        if not os.path.isfile(base_cfg_path):
            fail(f"Base era config does not exist: '{base_cfg_path}'")
        base_cfgs[era] = base_cfg_path

    for hist_yml in plots_include:
        hist_yml_path = os.path.join(config_path, hist_yml)
        if not os.path.isfile(hist_yml_path):
            fail(f"Histogram YAML file '{hist_yml}' not found in '{config_path}'")

    if not dry_run and not os.path.isfile(PLOTIT_BIN):
        fail(f"plotIt binary not found at '{PLOTIT_BIN}'")

    return base_cfgs


def main():
    args = parse_args()

    active_modes = sum([bool(args.DNN), bool(args.tauFF), bool(args.calcTauFF)])
    if active_modes > 1:
        fail("Mode flags -D (--dnn), -F (--ff), and -f (--calcff) are mutually exclusive.")
    if args.regions is not None and not args.calcTauFF:
        fail("--regions requires -f (--calcff).")

    plots_include = resolve_plots_include(args)
    plotit_flags = build_plotit_flags(args)
    abs_input = os.path.abspath(args.input.rstrip('/') or '.')

    if args.calcTauFF:
        if not args.channel:
            fail("-C/--channel ('muon' or 'electron') is required when using -f (--calcff).")
        channel = args.channel
        config_path = os.path.join(PLOTIT_ROOT, 'configs', f'Run3_{channel}')

        if FAKE_REGION_PATH_RE.search(abs_input):
            if args.regions is not None:
                fail("--regions can only be used in -f version-prefix mode (not when -I points to a single fake region directory).")
            if not os.path.isdir(abs_input):
                fail(f"Input directory does not exist: '{abs_input}'")
            era = resolve_era(args.year, abs_input)
            base_cfgs = _validate_plotit_env(config_path, [era], plots_include, args.dry_run)
            default_out = os.path.join(os.path.dirname(abs_input), f'figure_{era}{args.postfix}')
            abs_output = os.path.abspath(args.output) if args.output else default_out
            run_single(abs_input, abs_output, config_path, base_cfgs[era], channel, era, plots_include, plotit_flags, args, exclude_signals=True, strip_systematics=True)
            return

        if args.output:
            fail("-O/--output cannot be used with -f version mode; each region/era writes to <version>_fake_<region>/<channel>/figure_<era>.")
        regions = list(dict.fromkeys(args.regions)) if args.regions is not None else list(FAKE_REGIONS)
        targets = resolve_fake_targets(args.input, channel, args.year, regions, args.postfix)
        distinct_eras = list(dict.fromkeys(era for _, era, _, _ in targets))
        base_cfgs = _validate_plotit_env(config_path, distinct_eras, plots_include, args.dry_run)

        for region, era, era_dir, out_dir in targets:
            print(f"[INFO] [calcFF | {channel} | {era} | fake_{region}]")
            run_single(era_dir, out_dir, config_path, base_cfgs[era], channel, era, plots_include, plotit_flags, args, exclude_signals=True, strip_systematics=True)
        return

    if not os.path.isdir(abs_input):
        fail(f"Input directory does not exist: '{abs_input}'")

    channel = resolve_channel(args.channel, abs_input)
    era = resolve_era(args.year, abs_input)
    config_path = os.path.join(PLOTIT_ROOT, 'configs', f'Run3_{channel}')
    base_cfgs = _validate_plotit_env(config_path, [era], plots_include, args.dry_run)

    default_out = os.path.join(os.path.dirname(abs_input), f'figure_{os.path.basename(abs_input)}{args.postfix}')
    abs_output = os.path.abspath(args.output) if args.output else default_out
    run_single(abs_input, abs_output, config_path, base_cfgs[era], channel, era, plots_include, plotit_flags, args, exclude_signals=False)


if __name__ == '__main__':
    main()
