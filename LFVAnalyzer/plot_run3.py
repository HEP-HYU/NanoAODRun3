#!/usr/bin/env python3
"""
Generate combined Run 3 plotIt configurations (files_Run3.yml and config_Run3.yml)
across all 5 eras and run plotIt.
"""
# use (process plots): python plot_run3.py -I inputdir/channel
# use (DNN plots): python plot_run3.py -I inputdir/channel -D

import argparse
import os
import re
import subprocess
import sys

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(SCRIPT_DIR)
PLOTIT_ROOT = os.path.join(REPO_ROOT, 'plotIt')
PLOTIT_BIN = os.path.join(PLOTIT_ROOT, 'plotIt')

# [JME-03] Restored all JME systematics (11 reduced JES sources, 4 FlavorPure sources, jer, metUnclust) and tes.
# Retain only deprecated legacy aliases that do not exist in Run 3 histograms.
SYST_AVOID = ['muonhighscale']
SIGNAL_TYPE_RE = re.compile(r'^\s+type\s*:\s*[\'"]?signal[\'"]?(?:\s*(?:#.*)?)?$')
LOW_PT_QCD_FILES = {
    'hist_QCD_Pt10To30_EMEnriched.root',
    'hist_QCD_Pt15To20_EMEnriched.root',
    'hist_QCD_Pt20To30_EMEnriched.root',
    'hist_QCD_Pt15To20_MuEnriched.root',
    'hist_QCD_Pt20To30_MuEnriched.root',
}

# CMS_PdmVRun3Analysis.pdf (normtag_PHYSICS.json, Golden JSON, Era B excluded):
# 2022: 5.02+2.97=7.99, 2022EE: 5.81+17.78+3.09=26.68, 2023: 17.96, 2023BPix: 9.68, 2024: 109.82 fb^-1
ERAS = {
    '2022':     {'lumi': 7990.0,   'config_tag': '2022',     'dir_candidates': ['v12_2022', 'v2022', '2022']},
    '2022EE':   {'lumi': 26680.0,  'config_tag': '2022EE',   'dir_candidates': ['v12_2022EE', 'v2022EE', '2022EE']},
    '2023':     {'lumi': 17960.0,  'config_tag': '2023',     'dir_candidates': ['v12_2023', 'v2023', '2023']},
    '2023BPix': {'lumi': 9680.0,   'config_tag': '2023BPix', 'dir_candidates': ['v12_2023BPix', 'v12_2023_BPix', 'v2023_BPix', 'v2023BPix', '2023BPix', '2023_BPix']},
    '2024':     {'lumi': 109820.0, 'config_tag': '2024',     'dir_candidates': ['v15_2024', 'v2024', '2024']},
}
TOTAL_LUMI = sum(info['lumi'] for info in ERAS.values())  # 172130.0 pb^-1

# Unified inclusive signal cross sections (pb); TTto*ETau matched to NNLO muon values (hadronic / 0.67).
SIGNAL_SPECS = {
    'electron': [
        ('hist_TUETau-LFV-Scalar.root',     'LFVSTus', 0.1139,    'GLFVSTus', '1ST LFV \\cPqt\\cPqe\\tau Scalar',    1, None),
        ('hist_TUETau-LFV-Vector.root',     'LFVSTuv', 0.5267,    'GLFVSTuv', '1ST LFV \\cPqt\\cPqe\\tau Vector',    1, ('3ST LFV \\cPqt\\cPqu\\e\\tau Vector', 3)),
        ('hist_TUETau-LFV-Tensor.root',     'LFVSTut', 2.411,     'GLFVSTut', '2ST LFV \\cPqt\\cPqe\\tau Tensor',    1, None),
        ('hist_TCETau-LFV-Scalar.root',     'LFVSTcs', 0.007285,  'GLFVSTcs', '3ST LFV \\cPqt\\cPqc\\tau Scalar',    1, None),
        ('hist_TCETau-LFV-Vector.root',     'LFVSTcv', 0.03788,   'GLFVSTcv', '3ST LFV \\cPqt\\cPqc\\tau Vector',    3, None),
        ('hist_TCETau-LFV-Tensor.root',     'LFVSTct', 0.1874,    'GLFVSTct', '4ST LFV \\cPqt\\cPqc\\tau Tensor',    1, None),
        ('hist_TTtoUETau-LFV-Scalar.root',  'LFVTTus', 0.0044485, 'GLFVTTus', '5TT LFV \\cPqt\\cPqu\\e\\tau Scalar', 1, None),
        ('hist_TTtoUETau-LFV-Vector.root',  'LFVTTuv', 0.035555,  'GLFVTTuv', '5TT LFV \\cPqt\\cPqu\\e\\tau Vector', 2, ('7TT LFV \\cPqt\\cPqu\\e\\tau Vector', 4)),
        ('hist_TTtoUETau-LFV-Tensor.root',  'LFVTTut', 0.21333,   'GLFVTTut', '6TT LFV \\cPqt\\cPqu\\e\\tau Tensor', 1, None),
        ('hist_TTtoCETau-LFV-Scalar.root',  'LFVTTcs', 0.0044485, 'GLFVTTcs', '7TT LFV \\cPqt\\cPqc\\e\\tau Scalar', 1, None),
        ('hist_TTtoCETau-LFV-Vector.root',  'LFVTTcv', 0.035555,  'GLFVTTcv', '7TT LFV \\cPqt\\cPqc\\e\\tau Vector', 4, None),
        ('hist_TTtoCETau-LFV-Tensor.root',  'LFVTTct', 0.21333,   'GLFVTTct', '8TT LFV \\cPqt\\cPqc\\e\\tau Tensor', 1, None),
    ],
    'muon': [
        ('hist_TUMuTau-LFV-Scalar.root',    'LFVSTus', 0.1138,    'GLFVSTus', '1ST LFV \\cPqt\\cPqu\\mu\\tau Scalar',  1, None),
        ('hist_TUMuTau-LFV-Vector.root',    'LFVSTuv', 0.5266,    'GLFVSTuv', '1ST LFV \\cPqt\\cPqu\\mu\\tau Vector',  1, ('3ST LFV \\cPqt\\cPqu\\mu\\tau Vector', 3)),
        ('hist_TUMuTau-LFV-Tensor.root',    'LFVSTut', 2.401,     'GLFVSTut', '2ST LFV \\cPqt\\cPqu\\mu\\tau Tensor',  1, None),
        ('hist_TCMuTau-LFV-Scalar.root',    'LFVSTcs', 0.007288,  'GLFVSTcs', '3ST LFV \\cPqt\\cPqc\\mu\\tau Scalar',  1, None),
        ('hist_TCMuTau-LFV-Vector.root',    'LFVSTcv', 0.03753,   'GLFVSTcv', '3ST LFV \\cPqt\\cPqc\\mu\\tau Vector',  3, None),
        ('hist_TCMuTau-LFV-Tensor.root',    'LFVSTct', 0.1849,    'GLFVSTct', '4ST LFV \\cPqt\\cPqc\\mu\\tau Tensor',  1, None),
        ('hist_TTtoUMuTau-LFV-Scalar.root', 'LFVTTus', 0.0044485, 'GLFVTTus', '5TT LFV \\cPqt\\cPqu\\mu\\tau Scalar',  1, None),
        ('hist_TTtoUMuTau-LFV-Vector.root', 'LFVTTuv', 0.035555,  'GLFVTTuv', '5TT LFV \\cPqt\\cPqu\\mu\\tau Vector',  2, ('7TT LFV \\cPqt\\cPqu\\mu\\tau Vector', 4)),
        ('hist_TTtoUMuTau-LFV-Tensor.root', 'LFVTTut', 0.21333,   'GLFVTTut', '6TT LFV \\cPqt\\cPqu\\mu\\tau Tensor',  1, None),
        ('hist_TTtoCMuTau-LFV-Scalar.root', 'LFVTTcs', 0.0044485, 'GLFVTTcs', '7TT LFV \\cPqt\\cPqc\\mu\\tau Scalar',  1, None),
        ('hist_TTtoCMuTau-LFV-Vector.root', 'LFVTTcv', 0.035555,  'GLFVTTcv', '7TT LFV \\cPqt\\cPqc\\mu\\tau Vector',  4, None),
        ('hist_TTtoCMuTau-LFV-Tensor.root', 'LFVTTct', 0.21333,   'GLFVTTct', '8TT LFV \\cPqt\\cPqc\\mu\\tau Tensor',  1, None),
    ],
}

DATA_OR_SIGNAL_BASENAME_TAGS = ('LFV', 'SingleMuon', 'Egamma', 'EGamma', 'Muon')


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


def format_signal_blocks(channel, run3_dir_path, yield_only):
    blocks = ['']
    for fname, pretty, xsec, group, ygroup_all, order_all, plot_override in SIGNAL_SPECS[channel]:
        if yield_only:
            ygroup, order = ygroup_all, order_all
        elif plot_override is not None:
            ygroup, order = plot_override
        else:
            continue
        blocks.append(
            f"'{run3_dir_path}/{fname}':\n"
            f"  type: signal\n"
            f"  pretty-name: '{pretty}'\n"
            f"  cross-section: {xsec}\n"
            f"  generated-events: 1\n"
            f"  group: {group}\n"
            f"  yields-group: '{ygroup}'\n"
            f"  order: {order}\n"
        )
    return '\n'.join(blocks)


def collect_common_systs(config_path):
    common_syst_list = []
    for era_info in ERAS.values():
        cfg_file = os.path.join(config_path, f"config_{era_info['config_tag']}.yml")
        if not os.path.isfile(cfg_file):
            fail(f"Missing era config '{cfg_file}' needed to collect systematics.")
        with open(cfg_file, 'r') as f:
            in_syst = False
            for line in f:
                if not line[:1].isspace() and line.strip().startswith('systematics:'):
                    in_syst = True
                    continue
                if in_syst:
                    if line.strip() and not line[:1].isspace() and not line.strip().startswith('#'):
                        break
                    if any(s in line for s in SYST_AVOID):
                        continue
                    if line.startswith('  - '):
                        common_syst_list.append(line)
    if not common_syst_list:
        fail(f"Collected zero systematics from era configs in '{config_path}'.")
    ordered = list(dict.fromkeys(common_syst_list))
    return 'systematics:\n' + ''.join(ordered)


def resolve_era_dirs(dest_path, channel, postfix):
    resolved_eras = {}
    missing_eras = []
    for era_key, era_info in ERAS.items():
        matched_rel_path = None
        for cand in era_info['dir_candidates']:
            subpaths_to_test = [
                os.path.join(channel, f"{cand}_postprocess{postfix}"),
                os.path.join(channel, f"{cand}{postfix}"),
                os.path.join(channel, cand),
                f"{cand}_postprocess{postfix}",
                f"{cand}{postfix}",
                cand,
            ]
            for sub in subpaths_to_test:
                if os.path.isdir(os.path.join(dest_path, sub)):
                    matched_rel_path = sub
                    break
            if matched_rel_path:
                break
        if not matched_rel_path:
            missing_eras.append(era_key)
        else:
            resolved_eras[era_key] = matched_rel_path

    if missing_eras:
        fail(f"Missing histogram directories under '{dest_path}' for eras: {', '.join(missing_eras)}")
    return resolved_eras


def build_files_entries(files_yml_path, dest_path, era_rel_path, lumi, total_lumi, has_merged_run3, exclude_signals=False, skip_missing_low_pt_qcd=False):
    blocks = []
    current_header = None
    current_props = []

    with open(files_yml_path, 'r') as f:
        for lineno, raw_line in enumerate(f, 1):
            clean_line = raw_line.strip()
            if not clean_line or clean_line.startswith('#'):
                continue

            if not raw_line[0].isspace():
                if not (clean_line.startswith("'") and "':" in clean_line):
                    fail(f"Unexpected top-level YAML line in '{files_yml_path}' line {lineno}: {raw_line.rstrip()}")
                if current_header is not None:
                    blocks.append((current_header, current_props))
                current_header = raw_line
                current_props = []
            else:
                if current_header is None:
                    fail(f"Unexpected indented YAML line before any file entry in '{files_yml_path}' line {lineno}: {raw_line.rstrip()}")
                current_props.append(raw_line)

    if current_header is not None:
        blocks.append((current_header, current_props))

    out_blocks = []
    for header_line, prop_lines in blocks:
        clean_header = header_line.strip()
        is_signal = ('LFV' in clean_header) or any(SIGNAL_TYPE_RE.match(line) for line in prop_lines)
        if (has_merged_run3 or exclude_signals) and is_signal:
            continue

        fname = clean_header[1:clean_header.find("':")]
        basename = os.path.basename(fname)
        file_entry_path = os.path.normpath(os.path.join(dest_path, era_rel_path, basename))
        if skip_missing_low_pt_qcd and basename in LOW_PT_QCD_FILES and not os.path.isfile(file_entry_path):
            print(f"[WARNING] Missing low-pT QCD file (skipping): {file_entry_path}")
            continue
        block_lines = [f"'{file_entry_path}':\n"]
        if not any(tag in basename for tag in DATA_OR_SIGNAL_BASENAME_TAGS):
            block_lines.append(f"  scale: {float(lumi) / total_lumi}\n")
        for prop in prop_lines:
            block_lines.append(prop if prop.endswith('\n') else prop + '\n')
        out_blocks.append(''.join(block_lines))

    if not out_blocks:
        fail(f"No file entries remained after processing '{files_yml_path}'.")

    return '\n'.join(out_blocks) + '\n'


def resolve_plots_include(args):
    if args.tauFF:
        return ['histos_yield_S5.yml'] if args.yield_only else ['histos_FFapply.yml', 'histos_yield_S5.yml']
    if args.calcTauFF:
        return ['histos_FFcalc.yml']
    if args.DNN:
        return ['histos_dnn.yml']
    return ['histos_yield.yml'] if args.yield_only else ['histos_yield.yml', 'histos.yml', 'histos_reco.yml']


def build_config_run3_content(template_yml_path, common_syst, plots_include):
    out_lines = []
    found_root = False
    found_syst = False
    with open(template_yml_path, 'r') as f:
        for line in f:
            stripped = line.strip()
            if stripped.startswith('root:'):
                # Override empty root ('') with '#' so plotIt keeps absolute paths in files_Run3.yml intact
                found_root = True
                out_lines.append("  root: '#'\n")
            elif 'luminosity:' in line and not stripped.startswith('#'):
                out_lines.append(f"  luminosity: {TOTAL_LUMI:.1f}\n")
            elif stripped.startswith('systematics:'):
                found_syst = True
                break
            else:
                out_lines.append(line)
    if not found_root:
        fail(f"No 'root:' entry found in '{template_yml_path}'.")
    if not found_syst:
        fail(f"No 'systematics:' section found in '{template_yml_path}'.")
    out_lines.append(common_syst)
    include_literal = ', '.join(f"'{p}'" for p in plots_include)
    out_lines.append(f"\nplots:\n  include: [{include_literal}]\n")
    return ''.join(out_lines)


def parse_args():
    parser = argparse.ArgumentParser(description="Run plotIt for Run 3")
    parser.add_argument('-I', '--input', dest='input', type=str, default='process_0715_v4/electron', help="Input directory path")
    parser.add_argument('-C', '-ch', '--channel', dest='channel', type=str, default='', choices=['', 'electron', 'muon'], help="Channel: electron or muon (auto-detected from input if empty)")
    parser.add_argument('-D', dest='DNN', action='store_true', default=False, help="Run for DNN histograms")
    parser.add_argument('-F', '--ff', dest='tauFF', action='store_true', default=False, help="Run for tauFF histograms")
    parser.add_argument('-f', '--calcff', dest='calcTauFF', action='store_true', default=False, help="Run for calculate tauFF histograms (excluding signals & systematics)")
    parser.add_argument('-y', dest='yield_only', action='store_true', default=False, help="Run for yield histograms")
    parser.add_argument('-q', dest='qcd', action='store_true', default=False, help="Include QCD")
    parser.add_argument('--postfix', dest='postfix', type=str, default='', help="Add postfix to output here, to have rebinning for histograms")
    parser.add_argument('--dry-run', dest='dry_run', action='store_true', default=False, help="Validate inputs and print summary without modifying tracked configs or executing plotIt")
    return parser.parse_args()


def main():
    args = parse_args()

    if sum([bool(args.DNN), bool(args.tauFF), bool(args.calcTauFF)]) > 1:
        fail("Mode flags -D, -F, and -f are mutually exclusive.")

    dest_path = os.path.abspath(args.input.rstrip('/') or '.')
    if not os.path.isdir(dest_path):
        fail(f"Input directory does not exist: '{dest_path}'")

    channel = resolve_channel(args.channel, dest_path)
    config_path = os.path.join(PLOTIT_ROOT, 'configs', f'Run3_{channel}')
    if not os.path.isdir(config_path):
        fail(f"plotIt config directory does not exist: '{config_path}'")

    template_yml_path = os.path.join(config_path, 'template_Run3.yml')
    if not os.path.isfile(template_yml_path):
        fail(f"Template config does not exist: '{template_yml_path}'")

    plots_include = resolve_plots_include(args)
    for hist_yml in plots_include:
        if not os.path.isfile(os.path.join(config_path, hist_yml)):
            fail(f"Histogram YAML file '{hist_yml}' not found in '{config_path}'")

    postfix = args.postfix
    common_syst = "systematics: []\n" if args.calcTauFF else collect_common_systs(config_path)
    resolved_eras = resolve_era_dirs(dest_path, channel, postfix)

    run3_chan_dir = os.path.join(dest_path, channel, 'Run3' + postfix)
    run3_flat_dir = os.path.join(dest_path, 'Run3' + postfix)
    has_merged_run3 = os.path.isdir(run3_chan_dir) or os.path.isdir(run3_flat_dir)
    run3_dir_path = os.path.normpath(run3_chan_dir) if os.path.isdir(run3_chan_dir) else os.path.normpath(run3_flat_dir)

    string_for_files = ''
    for era_key, era_info in ERAS.items():
        files_yml_name = 'files24.yml' if '2024' in era_key else 'files.yml'
        files_yml_path = os.path.join(config_path, files_yml_name)
        if not os.path.isfile(files_yml_path):
            fail(f"Expected files YAML missing: '{files_yml_path}'")
        string_for_files += build_files_entries(
            files_yml_path=files_yml_path,
            dest_path=dest_path,
            era_rel_path=resolved_eras[era_key],
            lumi=era_info['lumi'],
            total_lumi=TOTAL_LUMI,
            has_merged_run3=has_merged_run3,
            exclude_signals=args.calcTauFF,
            skip_missing_low_pt_qcd=args.qcd,
        )

    files_run3_content = ''
    if has_merged_run3 and not args.calcTauFF:
        files_run3_content += format_signal_blocks(channel, run3_dir_path, args.yield_only)
    files_run3_content += string_for_files

    config_run3_content = build_config_run3_content(template_yml_path, common_syst, plots_include)

    out_fig_dir = os.path.join(dest_path, ('figure_run3_qcd' if args.qcd else 'figure_run3') + postfix)
    config_run3_yml = os.path.join(config_path, 'config_Run3.yml')
    files_run3_yml = os.path.join(config_path, 'files_Run3.yml')

    if args.yield_only and not args.qcd:
        plotit_flags = ['-y']
    elif args.yield_only and args.qcd:
        plotit_flags = ['-y', '-q']
    elif args.qcd:
        plotit_flags = ['-q', '-a']
    else:
        plotit_flags = ['-a']

    cmd = [PLOTIT_BIN, '-o', out_fig_dir, config_run3_yml] + plotit_flags

    if args.dry_run:
        if args.calcTauFF:
            sig_mode = "excluded (-f)"
        elif has_merged_run3:
            sig_mode = f"merged Run3 ({run3_dir_path})"
        else:
            sig_mode = "per-era"
        print(f"[DRY-RUN] Channel      : {channel}")
        print(f"[DRY-RUN] Total lumi   : {TOTAL_LUMI:.1f} pb^-1")
        print(f"[DRY-RUN] Merged Run3  : {has_merged_run3} ({run3_dir_path})")
        print(f"[DRY-RUN] Signals      : {sig_mode}")
        print(f"[DRY-RUN] Systematics  : {'cleared (-f, luminosity-error preserved)' if args.calcTauFF else 'combined (5 eras)'}")
        print(f"[DRY-RUN] plots        : {plots_include}")
        print(f"[DRY-RUN] files_Run3   : {len(files_run3_content.splitlines())} lines (not written in --dry-run)")
        print(f"[DRY-RUN] config_Run3  : {len(config_run3_content.splitlines())} lines (not written in --dry-run)")
        print(f"[DRY-RUN] Command      : {' '.join(cmd)}")
        return

    if not os.path.isfile(PLOTIT_BIN):
        fail(f"plotIt binary not found at '{PLOTIT_BIN}'")

    os.makedirs(out_fig_dir, exist_ok=True)
    with open(files_run3_yml, 'w') as fnew:
        fnew.write(files_run3_content)
    with open(config_run3_yml, 'w') as f1:
        f1.write(config_run3_content)

    print(f"[INFO] Running: {' '.join(cmd)}")
    proc = subprocess.run(cmd, check=False)
    if proc.returncode != 0:
        fail(f"plotIt failed (exit code {proc.returncode}) for Run 3 '{channel}'.", code=proc.returncode)


if __name__ == '__main__':
    main()
