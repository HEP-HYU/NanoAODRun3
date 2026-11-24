#!/usr/bin/env python3
"""
Automate per-era plotIt execution by dynamically selecting histogram include files
based on workflow stage (process, TauFF calculation/application, DNN, yields) without
manually editing config_<era>.yml.
"""
# Usage examples:
#   python plot_era.py -I process_0715_v4/muon/v15_2024
#   python plot_era.py -I process_0715_v4/muon/v15_2024 -f --all-fake
#   python plot_era.py -I process_0715_v4/muon/v15_2024_FF -F
#   python plot_era.py -I dnn_out/muon/v15_2024_postprocess -D

import argparse
import os
import re
import shlex
import subprocess
import sys
import tempfile

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(SCRIPT_DIR)
PLOTIT_ROOT = os.path.join(REPO_ROOT, 'plotIt')
PLOTIT_BIN = os.path.join(PLOTIT_ROOT, 'plotIt')

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

FAKE_REGIONS = ['fake_los', 'fake_lss', 'fake_tos', 'fake_tss']
ERA_TOKEN_RE = re.compile(r'(2023_?BPix|2022_?EE|2024|2023|2022)', re.IGNORECASE)


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


def replace_plots_section(base_cfg_path, plots_include):
    with open(base_cfg_path, 'r') as f:
        lines = f.readlines()

    new_lines = []
    in_plots = False
    replaced = False
    include_literal = ', '.join(f"'{p}'" for p in plots_include)
    replacement_block = f"plots:\n  include: [{include_literal}]\n\n"

    for line in lines:
        if not in_plots:
            if not line[:1].isspace() and re.match(r'^plots\s*:', line):
                if replaced:
                    fail(f"Multiple top-level 'plots:' sections found in '{base_cfg_path}'.")
                in_plots = True
                replaced = True
                new_lines.append(replacement_block)
            else:
                new_lines.append(line)
        else:
            stripped = line.strip()
            if not stripped or line[0].isspace() or stripped.startswith('#'):
                continue
            in_plots = False
            new_lines.append(line)

    if not replaced:
        fail(f"No top-level 'plots:' section found in '{base_cfg_path}'.")

    return ''.join(new_lines)


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


def run_single(input_dir, output_dir, config_path, base_cfg_path, era, plots_include, plotit_flags, args):
    if not os.path.isdir(input_dir):
        fail(f"Input directory does not exist: '{input_dir}'")

    cfg_content = replace_plots_section(base_cfg_path, plots_include)

    if args.dry_run:
        cmd = [PLOTIT_BIN, '-i', input_dir, '-o', output_dir, os.path.join(config_path, f'.tmp_plot_era_{era}.yml')] + plotit_flags
        print(f"[DRY-RUN] Config base : {base_cfg_path}")
        print(f"[DRY-RUN] Input (-i)  : {input_dir}")
        print(f"[DRY-RUN] Output (-o) : {output_dir}")
        print(f"[DRY-RUN] plots: include: {plots_include}")
        print(f"[DRY-RUN] Command     : {' '.join(shlex.quote(c) for c in cmd)}")
        if args.keep_tmp:
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
    tmp_path = None
    try:
        with tempfile.NamedTemporaryFile(
            mode='w',
            prefix=f'.tmp_plot_era_{era}_{os.getpid()}_',
            suffix='.yml',
            dir=config_path,
            delete=False,
        ) as tmp_f:
            tmp_f.write(cfg_content)
            tmp_path = tmp_f.name

        cmd = [PLOTIT_BIN, '-i', input_dir, '-o', output_dir, tmp_path] + plotit_flags
        print(f"[INFO] Running: {' '.join(shlex.quote(c) for c in cmd)}")
        proc = subprocess.run(cmd, check=False)
        if proc.returncode != 0:
            fail(f"plotIt failed (exit code {proc.returncode}) for era='{era}', input='{input_dir}'.", code=proc.returncode)
    finally:
        if tmp_path and os.path.exists(tmp_path) and not args.keep_tmp:
            os.remove(tmp_path)


def parse_args():
    parser = argparse.ArgumentParser(description="Run plotIt for a single Run 3 era without editing config files")
    parser.add_argument('-I', '--input', dest='input', type=str, default='.', help="Input histogram directory (or parent directory with --all-fake)")
    parser.add_argument('-O', '--output', dest='output', type=str, default='', help="Output figure directory (default: <parent>/figure_<input_basename><postfix>)")
    parser.add_argument('-Y', '--year', '--era', dest='year', type=str, default='', help="Era: 2022, 2022EE, 2023, 2023BPix, or 2024 (auto-detected from input if empty)")
    parser.add_argument('-C', '-ch', '--channel', dest='channel', type=str, default='', choices=['', 'electron', 'muon'], help="Channel: electron or muon (auto-detected from input if empty)")
    parser.add_argument('-D', '--dnn', dest='DNN', action='store_true', default=False, help="Run for DNN histograms (histos_dnn.yml)")
    parser.add_argument('-F', '--ff', dest='tauFF', action='store_true', default=False, help="Run for TauFF applied histograms (histos_FFapply.yml, histos_yield_S5.yml)")
    parser.add_argument('-f', '--calcff', dest='calcTauFF', action='store_true', default=False, help="Run for TauFF calculation histograms (histos_FFcalc.yml)")
    parser.add_argument('-y', '--yield-only', dest='yield_only', action='store_true', default=False, help="Use yield-only histogram list (histos_yield.yml or histos_yield_S5.yml with -F)")
    parser.add_argument('--plots', dest='plots', type=str, default='', help="Comma-separated histogram YAML files overriding mode selection")
    parser.add_argument('--all-fake', dest='all_fake', action='store_true', default=False, help="Iterate over fake_los, fake_lss, fake_tos, fake_tss inside -I (requires -f)")
    parser.add_argument('--no-yields', dest='no_yields', action='store_true', default=False, help="Omit -y when calling plotIt")
    parser.add_argument('--no-syst', dest='no_syst', action='store_true', default=False, help="Omit -s when calling plotIt")
    parser.add_argument('--no-qcd', dest='no_qcd', action='store_true', default=False, help="Omit -q when calling plotIt")
    parser.add_argument('--plotit-args', dest='plotit_args', type=str, default='', help="Additional flags passed verbatim to plotIt")
    parser.add_argument('--postfix', dest='postfix', type=str, default='', help="Postfix appended to default output directory name")
    parser.add_argument('--dry-run', dest='dry_run', action='store_true', default=False, help="Validate config and print plotIt command without executing")
    parser.add_argument('--keep-tmp', dest='keep_tmp', action='store_true', default=False, help="Keep generated temporary YAML config inside plotIt config directory")
    return parser.parse_args()


def main():
    args = parse_args()

    active_modes = sum([bool(args.DNN), bool(args.tauFF), bool(args.calcTauFF)])
    if active_modes > 1:
        fail("Mode flags -D (--dnn), -F (--ff), and -f (--calcff) are mutually exclusive.")
    if args.all_fake and not args.calcTauFF:
        fail("--all-fake requires -f (--calcff).")

    abs_input = os.path.abspath(args.input.rstrip('/') or '.')
    if not os.path.isdir(abs_input):
        fail(f"Input directory does not exist: '{abs_input}'")

    channel = resolve_channel(args.channel, abs_input)
    era = resolve_era(args.year, abs_input)
    plots_include = resolve_plots_include(args)
    plotit_flags = build_plotit_flags(args)

    config_path = os.path.join(PLOTIT_ROOT, 'configs', f'Run3_{channel}')
    if not os.path.isdir(config_path):
        fail(f"plotIt config directory does not exist: '{config_path}'")

    base_cfg_path = os.path.join(config_path, f'config_{era}.yml')
    if not os.path.isfile(base_cfg_path):
        fail(f"Base era config does not exist: '{base_cfg_path}'")

    for hist_yml in plots_include:
        hist_yml_path = os.path.join(config_path, hist_yml)
        if not os.path.isfile(hist_yml_path):
            fail(f"Histogram YAML file '{hist_yml}' not found in '{config_path}'")

    if not args.dry_run and not os.path.isfile(PLOTIT_BIN):
        fail(f"plotIt binary not found at '{PLOTIT_BIN}'")

    if args.all_fake:
        explicit_out = os.path.abspath(args.output) if args.output else ''
        for region in FAKE_REGIONS:
            region_in = os.path.join(abs_input, region)
            if not os.path.isdir(region_in):
                fail(f"Expected fake region subdirectory missing: '{region_in}'")
        for region in FAKE_REGIONS:
            region_in = os.path.join(abs_input, region)
            region_out = os.path.join(explicit_out, region) if explicit_out else os.path.join(abs_input, f'figure_{region}{args.postfix}')
            run_single(region_in, region_out, config_path, base_cfg_path, era, plots_include, plotit_flags, args)
    else:
        default_out = os.path.join(os.path.dirname(abs_input), f'figure_{os.path.basename(abs_input)}{args.postfix}')
        abs_output = os.path.abspath(args.output) if args.output else default_out
        run_single(abs_input, abs_output, config_path, base_cfg_path, era, plots_include, plotit_flags, args)


if __name__ == '__main__':
    main()
