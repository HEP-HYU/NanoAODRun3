#!/usr/bin/env python3
"""
validate_systematics.py
=======================
Systematic propagation validation framework for the LFV top analysis.

Checks:
  1. Branch existence  -- Central/Up/Down branches or files exist as expected.
  2. Object propagation -- Shape syst files have genuinely different object kinematics.
  3. Derived variable response -- Derived vars shift consistently with source objects.
  4. DNN input coverage -- All DNN input features are shifted in shape syst files.
  5. Weight syst isolation -- Weight branches vary; kinematic branches do NOT.
  6. Frozen variables -- Unaffected variables are bit-identical between files.
  7. Up/Down directionality -- Up and Down deviate in consistent opposite directions.

Usage
-----
  python validate_systematics.py \\
      --indir /path/to/process_output/ \\
      --sample TTto2L2Nu \\
      --year 2022 \\
      --ch muon \\
      [--nevents 5000] \\
      [--outfile syst_validation_report.json]

The --sample argument must be a substring of the ROOT filename
(e.g. "TTto2L2Nu", "ST_t-channel", "Tau-LFV").

Required Python packages: uproot, numpy
"""

import argparse
import json
import os
import sys
import warnings
from collections import defaultdict
from pathlib import Path

import numpy as np

try:
    import uproot
except ImportError:
    sys.exit("[ERROR] uproot not found.  Run:  pip install uproot")

# ─────────────────────────────────────────────────────────────
# Analysis-level knowledge  (mirrors TopLFVAnalyzer.cpp)
# ─────────────────────────────────────────────────────────────

# Shape systematics -> which files differ from nominal.
# 'objects'  : branches expected to change relative to nominal
# 'derived'  : downstream derived variables expected to change
# 'type'     : object | met  (used for frozen-var logic)
SHAPE_SYSTS = {
    # ── JES ───────────────────────────────────────────────────────────────────
    "jesAbsoluteup":
        {"type": "object",
         "objects": ["Jet1_pt","Jet1_mass","Jet2_pt","Jet2_mass","Jet3_pt","Jet3_mass","bJet1_pt"],
         "derived": ["chi2","chi2_SMW_mass","chi2_SMTop_mass","st_met"]},
    "jesAbsolutedown":
        {"type": "object",
         "objects": ["Jet1_pt","Jet1_mass","Jet2_pt","Jet2_mass","Jet3_pt","Jet3_mass","bJet1_pt"],
         "derived": ["chi2","chi2_SMW_mass","chi2_SMTop_mass","st_met"]},
    "jesBBEC1up":
        {"type": "object",
         "objects": ["Jet1_pt","Jet2_pt","Jet3_pt"],
         "derived": ["chi2","st_met"]},
    "jesBBEC1down":
        {"type": "object",
         "objects": ["Jet1_pt","Jet2_pt","Jet3_pt"],
         "derived": ["chi2","st_met"]},
    "jesRelativeBalup":
        {"type": "object",
         "objects": ["Jet1_pt","Jet2_pt","Jet3_pt"],
         "derived": ["chi2","st_met"]},
    "jesRelativeBaldown":
        {"type": "object",
         "objects": ["Jet1_pt","Jet2_pt","Jet3_pt"],
         "derived": ["chi2","st_met"]},
    "jesFlavorPureGluonup":
        {"type": "object",
         "objects": ["Jet1_pt","Jet2_pt","Jet3_pt"],
         "derived": ["chi2","st_met"]},
    "jesFlavorPureGluondown":
        {"type": "object",
         "objects": ["Jet1_pt","Jet2_pt","Jet3_pt"],
         "derived": ["chi2","st_met"]},
    "jesFlavorPureQuarkup":
        {"type": "object",
         "objects": ["Jet1_pt","Jet2_pt","Jet3_pt"],
         "derived": ["chi2","st_met"]},
    "jesFlavorPureQuarkdown":
        {"type": "object",
         "objects": ["Jet1_pt","Jet2_pt","Jet3_pt"],
         "derived": ["chi2","st_met"]},
    "jesFlavorPureCharmup":
        {"type": "object",
         "objects": ["Jet1_pt","Jet2_pt","Jet3_pt"],
         "derived": ["chi2","st_met"]},
    "jesFlavorPureCharmdown":
        {"type": "object",
         "objects": ["Jet1_pt","Jet2_pt","Jet3_pt"],
         "derived": ["chi2","st_met"]},
    "jesFlavorPureBottomup":
        {"type": "object",
         "objects": ["Jet1_pt","Jet2_pt","Jet3_pt"],
         "derived": ["chi2","st_met"]},
    "jesFlavorPureBottomdown":
        {"type": "object",
         "objects": ["Jet1_pt","Jet2_pt","Jet3_pt"],
         "derived": ["chi2","st_met"]},
    # ── JER ───────────────────────────────────────────────────────────────────
    "jerup":
        {"type": "object",
         "objects": ["Jet1_pt","Jet2_pt","Jet3_pt"],
         "derived": ["chi2","st_met"]},
    "jerdown":
        {"type": "object",
         "objects": ["Jet1_pt","Jet2_pt","Jet3_pt"],
         "derived": ["chi2","st_met"]},
    # ── TES ───────────────────────────────────────────────────────────────────
    "tesup":
        {"type": "object",
         "objects": ["Tau1_pt","Tau1_mass"],
         "derived": ["leptau_mass","leptau_dR","leptau_dEta","leptau_dPhi",
                     "chi2","st_met","lepMET_mt"]},
    "tesdown":
        {"type": "object",
         "objects": ["Tau1_pt","Tau1_mass"],
         "derived": ["leptau_mass","leptau_dR","leptau_dEta","leptau_dPhi",
                     "chi2","st_met","lepMET_mt"]},
    # ── MET Unclust ───────────────────────────────────────────────────────────
    "metUnclustup":
        {"type": "met",
         "objects": ["PuppiMET_pt","PuppiMET_phi"],
         "derived": ["lepMET_mt","st_met"]},
    "metUnclustdown":
        {"type": "met",
         "objects": ["PuppiMET_pt","PuppiMET_phi"],
         "derived": ["lepMET_mt","st_met"]},
    # ── Muon GE Scale ─────────────────────────────────────────────────────────
    "muonhighscaleup":
        {"type": "object",
         "objects": ["Muon1_pt","PuppiMET_pt"],
         "derived": ["leptau_mass","lepMET_mt","st_met"]},
    "muonhighscaledown":
        {"type": "object",
         "objects": ["Muon1_pt","PuppiMET_pt"],
         "derived": ["leptau_mass","lepMET_mt","st_met"]},
}

# Up/Down pairs for directionality check (base name without up/down suffix)
SHAPE_SYST_PAIRS = [
    "jesAbsolute", "jesBBEC1", "jesRelativeBal",
    "jesFlavorPureGluon", "jesFlavorPureQuark",
    "jesFlavorPureCharm", "jesFlavorPureBottom",
    "jer", "tes", "metUnclust", "muonhighscale",
]

# Weight systematics stored as branches inside the NOMINAL file
WEIGHT_SYST_BRANCHES = [
    "eventWeight__puup", "eventWeight__pudown",
    "eventWeight__topptup", "eventWeight__topptdown",
    "eventWeight__muidup", "eventWeight__muiddown",
    "eventWeight__muisoup", "eventWeight__muisodown",
    "eventWeight__mutrgup", "eventWeight__mutrgdown",
    "eventWeight__btaghfup", "eventWeight__btaghfdown",
    "eventWeight__btaglfup", "eventWeight__btaglfdown",
    "eventWeight__btaghfstats1up", "eventWeight__btaghfstats1down",
    "eventWeight__btaghfstats2up", "eventWeight__btaghfstats2down",
    "eventWeight__btaglfstats1up", "eventWeight__btaglfstats1down",
    "eventWeight__btaglfstats2up", "eventWeight__btaglfstats2down",
    "eventWeight__btagcferr1up",   "eventWeight__btagcferr1down",
    "eventWeight__btagcferr2up",   "eventWeight__btagcferr2down",
    "eventWeight__tauidjetUncert0up",    "eventWeight__tauidjetUncert0down",
    "eventWeight__tauidjetUncert1up",    "eventWeight__tauidjetUncert1down",
    "eventWeight__tauidjetSystallerasup","eventWeight__tauidjetSystallerasdown",
    "eventWeight__tauidelup",  "eventWeight__tauideldown",
    "eventWeight__tauidmuup",  "eventWeight__tauidmudown",
]

# Up/Down pairs for weight directionality check
WEIGHT_SYST_PAIRS = [
    ("eventWeight__puup",           "eventWeight__pudown"),
    ("eventWeight__topptup",        "eventWeight__topptdown"),
    ("eventWeight__muidup",         "eventWeight__muiddown"),
    ("eventWeight__muisoup",        "eventWeight__muisodown"),
    ("eventWeight__mutrgup",        "eventWeight__mutrgdown"),
    ("eventWeight__btaghfup",       "eventWeight__btaghfdown"),
    ("eventWeight__btaglfup",       "eventWeight__btaglfdown"),
    ("eventWeight__btaghfstats1up", "eventWeight__btaghfstats1down"),
    ("eventWeight__btaghfstats2up", "eventWeight__btaghfstats2down"),
    ("eventWeight__btaglfstats1up", "eventWeight__btaglfstats1down"),
    ("eventWeight__btaglfstats2up", "eventWeight__btaglfstats2down"),
    ("eventWeight__btagcferr1up",   "eventWeight__btagcferr1down"),
    ("eventWeight__btagcferr2up",   "eventWeight__btagcferr2down"),
]

# Kinematic branches that must NOT change for weight-only systematics
FROZEN_KINEMATIC_BRANCHES = [
    "Tau1_pt","Tau1_mass","Tau1_eta",
    "Jet1_pt","Jet1_mass","Jet1_eta",
    "Jet2_pt","Jet2_mass","Jet2_eta",
    "Jet3_pt","Jet3_mass","Jet3_eta",
    "bJet1_pt","bJet1_mass","bJet1_eta",
    "leptau_mass","leptau_dEta","leptau_dPhi","leptau_dR",
    "chi2","chi2_SMW_mass","chi2_SMTop_mass",
    "chi2_wqq_dEta","chi2_wqq_dPhi","chi2_wqq_dR",
    "PuppiMET_pt","PuppiMET_phi",
    "lepMET_mt","st_met",
]

# DNN input features per channel
DNN_FEATURES = {
    "muon": [
        "Muon1_pt","Muon1_eta",
        "Tau1_pt","Tau1_mass","Tau1_eta",
        "Jet1_pt","Jet1_mass","Jet1_eta","Jet1_btagPNetB",
        "Jet2_pt","Jet2_mass","Jet2_eta","Jet2_btagPNetB",
        "Jet3_pt","Jet3_mass","Jet3_eta","Jet3_btagPNetB",
        "chi2","chi2_SMW_mass","chi2_SMTop_mass",
        "chi2_wqq_dEta","chi2_wqq_dPhi","chi2_wqq_dR",
        "leptau_mass","leptau_dEta","leptau_dPhi","leptau_dR",
        "PuppiMET_pt",
    ],
    "electron": [
        "Electron1_pt","Electron1_eta",
        "Tau1_pt","Tau1_mass","Tau1_eta",
        "Jet1_pt","Jet1_mass","Jet1_eta","Jet1_btagPNetB",
        "Jet2_pt","Jet2_mass","Jet2_eta","Jet2_btagPNetB",
        "Jet3_pt","Jet3_mass","Jet3_eta","Jet3_btagPNetB",
        "chi2","chi2_SMW_mass","chi2_SMTop_mass",
        "chi2_wqq_dEta","chi2_wqq_dPhi","chi2_wqq_dR",
        "leptau_mass","leptau_dEta","leptau_dPhi","leptau_dR",
        "PuppiMET_pt",
    ],
}

TREENAME = "Events"


# ─────────────────────────────────────────────────────────────
# Data structures
# ─────────────────────────────────────────────────────────────

class CheckResult:
    PASS = "PASS"
    FAIL = "FAIL"
    WARN = "WARN"
    SKIP = "SKIP"

    def __init__(self, status, name, detail=""):
        self.status = status
        self.name   = name
        self.detail = detail

    def __repr__(self):
        return f"[{self.status}] {self.name}: {self.detail}"


# ─────────────────────────────────────────────────────────────
# I/O helpers
# ─────────────────────────────────────────────────────────────

def _load_arrays(fpath, branches, n=-1):
    """Read branches from TTree 'Events'. Returns {branch: np.array | None}."""
    out = {}
    if not fpath or not os.path.exists(fpath):
        return {}
    try:
        with uproot.open(fpath) as f:
            if TREENAME not in f:
                return {}
            tree = f[TREENAME]
            avail   = set(tree.keys())
            to_read = [b for b in branches if b in avail]
            missing = [b for b in branches if b not in avail]
            stop    = n if n > 0 else None
            if to_read:
                arrs = tree.arrays(to_read, entry_stop=stop, library="np")
                out.update(arrs)
            for b in missing:
                out[b] = None
    except Exception as e:
        warnings.warn(f"[uproot] {fpath}: {e}")
    return out


def _all_branches(fpath):
    try:
        with uproot.open(fpath) as f:
            return set(f[TREENAME].keys()) if TREENAME in f else set()
    except Exception:
        return set()


# ─────────────────────────────────────────────────────────────
# Statistical helpers
# ─────────────────────────────────────────────────────────────

def _mean_frac_diff(a, b):
    """Mean |a-b|/max(|a|,|b|,1e-6), nan-safe."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        af = a.astype(float)
        bf = b.astype(float)
        denom = np.maximum(np.maximum(np.abs(af), np.abs(bf)), 1e-6)
        return float(np.nanmean(np.abs(af - bf) / denom))


def _arrays_identical(a, b, rtol=1e-5):
    if a is None or b is None:
        return True
    n = min(len(a), len(b))
    if n == 0:
        return True
    return np.allclose(a[:n].astype(float), b[:n].astype(float),
                       rtol=rtol, atol=1e-9, equal_nan=True)


def _mean_shift(a, b):
    n = min(len(a), len(b))
    return float(np.nanmean(b[:n].astype(float) - a[:n].astype(float)))


# ─────────────────────────────────────────────────────────────
# File discovery
# ─────────────────────────────────────────────────────────────

def discover_sample_files(indir, sample):
    """
    Returns {syst_suffix_or_'nominal': filepath} for all files matching sample.
    Filename pattern from process.py:  hist_{dataset_name}{src}.root
    where src = "" for nominal, "__jesAbsoluteup" for shape systs.
    """
    result = {}
    for f in sorted(Path(indir).glob("*.root")):
        if sample not in f.stem:
            continue
        parts = f.stem.split("__")   # e.g. ["hist_TTto2L2Nu", "jesAbsoluteup"]
        if len(parts) == 1:
            result["nominal"] = str(f)
        else:
            result[parts[-1]] = str(f)
    return result


# ─────────────────────────────────────────────────────────────
# CHECK 1 – Branch / file existence
# ─────────────────────────────────────────────────────────────

def check_branch_existence(nominal_file, sample_files, ch):
    results = []

    # 1a. Nominal file accessible
    if not nominal_file or not os.path.exists(nominal_file):
        results.append(CheckResult(CheckResult.FAIL, "CHECK1_nominal_file",
                                   f"Not found: {nominal_file}"))
        return results
    results.append(CheckResult(CheckResult.PASS, "CHECK1_nominal_file"))

    nom_br = _all_branches(nominal_file)

    # 1b. Core kinematic branches
    for b in ["eventWeight","Tau1_pt","Jet1_pt","chi2","leptau_mass","PuppiMET_pt","st_met"]:
        st = CheckResult.PASS if b in nom_br else CheckResult.FAIL
        results.append(CheckResult(st, f"CHECK1_nominal_has_{b}",
                                   "" if st == CheckResult.PASS else "MISSING"))

    # 1c. Weight syst branches
    missing_wt = [b for b in WEIGHT_SYST_BRANCHES if b not in nom_br]
    if missing_wt:
        results.append(CheckResult(CheckResult.WARN, "CHECK1_weight_syst_branches",
                                   f"{len(missing_wt)} missing: {missing_wt[:6]}"
                                   f"{'...' if len(missing_wt) > 6 else ''}"))
    else:
        results.append(CheckResult(CheckResult.PASS, "CHECK1_weight_syst_branches",
                                   f"All {len(WEIGHT_SYST_BRANCHES)} weight syst branches present"))

    # 1d. Shape syst files
    found, miss = [], []
    for s in SHAPE_SYSTS:
        (found if s in sample_files else miss).append(s)
    if miss:
        results.append(CheckResult(CheckResult.WARN, "CHECK1_shape_syst_files",
                                   f"{len(miss)} files missing: {miss[:5]}"
                                   f"{'...' if len(miss) > 5 else ''}"))
    else:
        results.append(CheckResult(CheckResult.PASS, "CHECK1_shape_syst_files",
                                   f"All {len(SHAPE_SYSTS)} shape syst files found"))

    # 1e. DNN features in nominal
    dnn_feats = DNN_FEATURES.get(ch, DNN_FEATURES["muon"])
    miss_dnn  = [f for f in dnn_feats if f not in nom_br]
    if miss_dnn:
        results.append(CheckResult(CheckResult.FAIL, "CHECK1_dnn_features",
                                   f"Missing in nominal: {miss_dnn}"))
    else:
        results.append(CheckResult(CheckResult.PASS, "CHECK1_dnn_features",
                                   f"All {len(dnn_feats)} DNN features present"))

    return results


# ─────────────────────────────────────────────────────────────
# CHECK 2+3+4 – Shape syst propagation
# ─────────────────────────────────────────────────────────────

def check_shape_propagation(nominal_file, sample_files, ch, nevents):
    results = []

    # Collect all branches we might need in one read
    all_needed = set(FROZEN_KINEMATIC_BRANCHES + DNN_FEATURES.get(ch, DNN_FEATURES["muon"]))
    for info in SHAPE_SYSTS.values():
        all_needed.update(info["objects"] + info["derived"])
    all_needed = list(all_needed)

    nom = _load_arrays(nominal_file, all_needed, nevents)

    for syst_name, info in SHAPE_SYSTS.items():
        if syst_name not in sample_files:
            results.append(CheckResult(CheckResult.SKIP,
                                       f"CHECK2_object_{syst_name}", "file missing"))
            results.append(CheckResult(CheckResult.SKIP,
                                       f"CHECK3_derived_{syst_name}", "file missing"))
            results.append(CheckResult(CheckResult.SKIP,
                                       f"CHECK4_dnn_{syst_name}", "file missing"))
            continue

        syst = _load_arrays(sample_files[syst_name], all_needed, nevents)
        if not syst:
            for tag in ["CHECK2_object","CHECK3_derived","CHECK4_dnn"]:
                results.append(CheckResult(CheckResult.FAIL,
                                           f"{tag}_{syst_name}", "could not read file"))
            continue

        # ── CHECK 2: expected object branches are shifted ─────────────────────
        details2, all_shifted = [], True
        for br in info["objects"]:
            a, b = nom.get(br), syst.get(br)
            if a is None or b is None:
                details2.append(f"{br}=MISSING"); all_shifted = False; continue
            n = min(len(a), len(b))
            fd = _mean_frac_diff(a[:n], b[:n])
            if fd < 1e-5:
                details2.append(f"{br}=IDENTICAL({fd:.2e})"); all_shifted = False
            else:
                details2.append(f"{br}=OK({fd:.3f})")
        results.append(CheckResult(
            CheckResult.PASS if all_shifted else CheckResult.FAIL,
            f"CHECK2_object_{syst_name}", " | ".join(details2)))

        # ── CHECK 3: derived variables respond to the shift ───────────────────
        details3, all_derived = [], True
        for br in info["derived"]:
            a, b = nom.get(br), syst.get(br)
            if a is None or b is None:
                details3.append(f"{br}=MISSING"); continue
            n = min(len(a), len(b))
            fd = _mean_frac_diff(a[:n], b[:n])
            if fd < 1e-6:
                details3.append(f"{br}=FROZEN({fd:.2e})"); all_derived = False
            else:
                details3.append(f"{br}=OK({fd:.3f})")
        results.append(CheckResult(
            CheckResult.PASS if all_derived else CheckResult.FAIL,
            f"CHECK3_derived_{syst_name}", " | ".join(details3)))

        # ── CHECK 4: DNN features affected by this syst are shifted ───────────
        dnn_feats = DNN_FEATURES.get(ch, DNN_FEATURES["muon"])
        affected  = [f for f in dnn_feats
                     if f in info["objects"] or f in info["derived"]]
        if not affected:
            results.append(CheckResult(CheckResult.SKIP,
                                       f"CHECK4_dnn_{syst_name}", "no DNN features affected"))
            continue
        details4, dnn_ok = [], True
        for feat in affected:
            a, b = nom.get(feat), syst.get(feat)
            if a is None or b is None:
                details4.append(f"{feat}=MISSING"); dnn_ok = False; continue
            n = min(len(a), len(b))
            fd = _mean_frac_diff(a[:n], b[:n])
            if fd < 1e-6:
                details4.append(f"{feat}=FROZEN"); dnn_ok = False
            else:
                details4.append(f"{feat}=OK")
        results.append(CheckResult(
            CheckResult.PASS if dnn_ok else CheckResult.FAIL,
            f"CHECK4_dnn_{syst_name}", " | ".join(details4)))

    return results


# ─────────────────────────────────────────────────────────────
# CHECK 5 – Weight syst isolation
# ─────────────────────────────────────────────────────────────

def check_weight_isolation(nominal_file, sample_files, nevents):
    results = []

    # 5a. Weight branches in nominal differ from nominal eventWeight
    all_wt_br = ["eventWeight"] + WEIGHT_SYST_BRANCHES
    nom = _load_arrays(nominal_file,
                       all_wt_br + FROZEN_KINEMATIC_BRANCHES[:6], nevents)
    nom_ew = nom.get("eventWeight")
    if nom_ew is None:
        results.append(CheckResult(CheckResult.FAIL, "CHECK5_eventWeight_missing", ""))
        return results

    identical, differing = [], 0
    for wb in WEIGHT_SYST_BRANCHES:
        v = nom.get(wb)
        if v is None:
            continue
        n = min(len(nom_ew), len(v))
        if _arrays_identical(nom_ew[:n], v[:n]):
            identical.append(wb)
        else:
            differing += 1
    if identical:
        results.append(CheckResult(CheckResult.WARN,
                                   "CHECK5_weight_branches_differ",
                                   f"{len(identical)} weight branches identical to eventWeight: "
                                   f"{identical[:4]}{'...' if len(identical) > 4 else ''}"))
    else:
        results.append(CheckResult(CheckResult.PASS, "CHECK5_weight_branches_differ",
                                   f"All {differing} weight branches differ from nominal"))

    # 5b. btag ext_syst file: kinematics must be frozen
    btag_files = {k: v for k, v in sample_files.items() if "btag" in k.lower()}
    for bsyst, bfile in btag_files.items():
        s_arrs = _load_arrays(bfile, FROZEN_KINEMATIC_BRANCHES[:8], nevents)
        fail_br, ok_br = [], []
        for br in FROZEN_KINEMATIC_BRANCHES[:8]:
            a, b = nom.get(br), s_arrs.get(br)
            if a is None or b is None:
                continue
            n = min(len(a), len(b))
            if _arrays_identical(a[:n], b[:n]):
                ok_br.append(br)
            else:
                fail_br.append(br)
        if fail_br:
            results.append(CheckResult(CheckResult.FAIL,
                                       f"CHECK5_btag_kin_frozen_{bsyst}",
                                       f"Changed (should be frozen): {fail_br}"))
        elif ok_br:
            results.append(CheckResult(CheckResult.PASS,
                                       f"CHECK5_btag_kin_frozen_{bsyst}",
                                       f"All {len(ok_br)} kinematic branches frozen"))

    return results


# ─────────────────────────────────────────────────────────────
# CHECK 6 – Frozen variables for each shape syst
# ─────────────────────────────────────────────────────────────

def check_frozen_variables(nominal_file, sample_files, nevents):
    results = []

    # Branches expected to be FROZEN for each syst category
    frozen_map = {
        "tes":           ["Jet1_pt","Jet2_pt","Jet3_pt","Muon1_pt","Electron1_pt",
                          "PuppiMET_pt","PuppiMET_phi"],
        "jes":           ["Tau1_pt","Tau1_mass","Muon1_pt","Electron1_pt"],
        "jer":           ["Tau1_pt","Tau1_mass","Muon1_pt","Electron1_pt"],
        "metUnclust":    ["Tau1_pt","Tau1_mass","Jet1_pt","Jet2_pt","Jet3_pt",
                          "Muon1_pt","Electron1_pt","leptau_mass"],
        "muonhighscale": ["Tau1_pt","Tau1_mass","Jet1_pt","Jet2_pt","Jet3_pt"],
    }

    nom = _load_arrays(nominal_file, FROZEN_KINEMATIC_BRANCHES, nevents)

    for syst_name in SHAPE_SYSTS:
        if syst_name not in sample_files:
            continue

        # Match syst category
        cat = "jes"
        for key in frozen_map:
            if syst_name.startswith(key):
                cat = key
                break

        frozen_br = frozen_map[cat]
        s_arrs    = _load_arrays(sample_files[syst_name], frozen_br, nevents)

        fail, ok = [], []
        for br in frozen_br:
            a, b = nom.get(br), s_arrs.get(br)
            if a is None or b is None:
                continue
            n = min(len(a), len(b))
            if _arrays_identical(a[:n], b[:n]):
                ok.append(br)
            else:
                fd = _mean_frac_diff(a[:n], b[:n])
                fail.append(f"{br}(diff={fd:.3f})")

        if fail:
            results.append(CheckResult(CheckResult.FAIL,
                                       f"CHECK6_frozen_{syst_name}",
                                       f"Unexpectedly changed: {fail}"))
        elif ok:
            results.append(CheckResult(CheckResult.PASS,
                                       f"CHECK6_frozen_{syst_name}",
                                       f"All {len(ok)} unaffected branches frozen"))
        else:
            results.append(CheckResult(CheckResult.SKIP,
                                       f"CHECK6_frozen_{syst_name}", "no branches available"))

    return results


# ─────────────────────────────────────────────────────────────
# CHECK 7 – Up/Down directionality
# ─────────────────────────────────────────────────────────────

def check_directionality(nominal_file, sample_files, nevents):
    results = []

    # 7a. Shape syst pairs
    for base in SHAPE_SYST_PAIRS:
        # Find up/down keys (handle year-embedded variants like jesAbsolute_2022up)
        up_key   = next((k for k in sample_files if k.startswith(base) and k.endswith("up")),   None)
        down_key = next((k for k in sample_files if k.startswith(base) and k.endswith("down")), None)

        if up_key is None or down_key is None:
            results.append(CheckResult(CheckResult.SKIP,
                                       f"CHECK7_shape_{base}",
                                       f"up={up_key} down={down_key}"))
            continue

        # Primary branch for this syst type
        up_info = SHAPE_SYSTS.get(up_key) or SHAPE_SYSTS.get(base + "up")
        if up_info is None:
            continue
        primary = up_info["objects"][0]

        nom_v  = _load_arrays(nominal_file,            [primary], nevents).get(primary)
        up_v   = _load_arrays(sample_files[up_key],    [primary], nevents).get(primary)
        down_v = _load_arrays(sample_files[down_key],  [primary], nevents).get(primary)

        if any(v is None for v in [nom_v, up_v, down_v]):
            results.append(CheckResult(CheckResult.SKIP,
                                       f"CHECK7_shape_{base}",
                                       f"branch {primary} missing in one file"))
            continue

        n = min(len(nom_v), len(up_v), len(down_v))
        du = _mean_shift(nom_v[:n], up_v[:n])
        dd = _mean_shift(nom_v[:n], down_v[:n])
        fu = _mean_frac_diff(nom_v[:n], up_v[:n])
        fd_val = _mean_frac_diff(nom_v[:n], down_v[:n])

        if fu < 1e-6 and fd_val < 1e-6:
            results.append(CheckResult(CheckResult.FAIL,
                                       f"CHECK7_shape_{base}",
                                       f"Both Up & Down identical to nominal for {primary}"))
        elif du * dd < 0:
            results.append(CheckResult(CheckResult.PASS,
                                       f"CHECK7_shape_{base}",
                                       f"{primary}: DeltaUp={du:+.4f} DeltaDown={dd:+.4f} (opposite signs OK)"))
        else:
            results.append(CheckResult(CheckResult.WARN,
                                       f"CHECK7_shape_{base}",
                                       f"{primary}: DeltaUp={du:+.4f} DeltaDown={dd:+.4f} (same sign)"))

    # 7b. Weight syst pairs (all from nominal file)
    all_wt_br = ["eventWeight"] + [b for pair in WEIGHT_SYST_PAIRS for b in pair]
    nom = _load_arrays(nominal_file, all_wt_br, nevents)
    nom_ew = nom.get("eventWeight")
    if nom_ew is not None:
        for (up_br, dn_br) in WEIGHT_SYST_PAIRS:
            up_v = nom.get(up_br)
            dn_v = nom.get(dn_br)
            if up_v is None or dn_v is None:
                continue
            n = min(len(nom_ew), len(up_v), len(dn_v))
            du = _mean_shift(nom_ew[:n], up_v[:n])
            dd = _mean_shift(nom_ew[:n], dn_v[:n])
            label = up_br.replace("eventWeight__", "").rstrip("up")
            if du * dd < 0:
                results.append(CheckResult(CheckResult.PASS,
                                           f"CHECK7_weight_{label}",
                                           f"DeltaUp={du:+.4f} DeltaDown={dd:+.4f}"))
            else:
                results.append(CheckResult(CheckResult.WARN,
                                           f"CHECK7_weight_{label}",
                                           f"Same-sign: DeltaUp={du:+.4f} DeltaDown={dd:+.4f}"))

    return results


# ─────────────────────────────────────────────────────────────
# Report
# ─────────────────────────────────────────────────────────────

def render_report(all_results, outfile=None):
    from collections import Counter
    counts  = Counter(r.status for r in all_results)
    total   = len(all_results)
    sym_map = {CheckResult.PASS: "OK", CheckResult.FAIL: "FAIL",
               CheckResult.WARN: "WARN", CheckResult.SKIP: "SKIP"}

    sep = "=" * 80
    print(sep)
    print("  SYSTEMATIC PROPAGATION VALIDATION REPORT")
    print(sep)
    print(f"  Total checks : {total}")
    print(f"  PASS  : {counts[CheckResult.PASS]}")
    print(f"  WARN  : {counts[CheckResult.WARN]}")
    print(f"  FAIL  : {counts[CheckResult.FAIL]}")
    print(f"  SKIP  : {counts[CheckResult.SKIP]}")
    print(sep)

    for status in [CheckResult.FAIL, CheckResult.WARN, CheckResult.SKIP, CheckResult.PASS]:
        group = [r for r in all_results if r.status == status]
        if not group:
            continue
        print(f"\n[{status}]  ({len(group)} checks)")
        print("-" * 78)
        for r in group:
            tail = f"\n      -> {r.detail}" if r.detail else ""
            print(f"  [{sym_map[status]}]  {r.name}{tail}")

    if outfile:
        report = {
            "summary": dict(counts),
            "checks":  [{"status": r.status, "name": r.name, "detail": r.detail}
                        for r in all_results],
        }
        with open(outfile, "w") as fh:
            json.dump(report, fh, indent=2)
        print(f"\n[INFO] JSON report written to: {outfile}")


# ─────────────────────────────────────────────────────────────
# Gap discovery
# ─────────────────────────────────────────────────────────────

def discover_propagation_gaps(all_results):
    """
    Scan FAIL/WARN results and produce human-readable gap descriptions.
    """
    gaps = []
    fail_warn = [r for r in all_results if r.status in (CheckResult.FAIL, CheckResult.WARN)]

    by_check = defaultdict(list)
    for r in fail_warn:
        # Extract check number from name like CHECK2_object_jesAbsoluteup
        prefix = "_".join(r.name.split("_")[:2])
        by_check[prefix].append(r)

    desc = {
        "CHECK1": "[GAP-Existence]    Missing branches or files",
        "CHECK2": "[GAP-ObjectShift]  Object kinematics not shifted in syst file",
        "CHECK3": "[GAP-DerivedProp]  Derived variables frozen despite object shift",
        "CHECK4": "[GAP-DNNInput]     DNN input feature not shifted",
        "CHECK5": "[GAP-WeightIsolation] Weight syst leaks into kinematics",
        "CHECK6": "[GAP-FrozenVars]   Unexpected kinematic change in unrelated syst",
        "CHECK7": "[GAP-Direction]    Up/Down shift in same direction",
    }

    for prefix in sorted(by_check):
        res_list = by_check[prefix]
        check_id = "_".join(prefix.split("_")[:2])
        header   = desc.get(check_id, f"[GAP-{check_id}]")
        systs    = [r.name for r in res_list]
        gaps.append(f"{header}\n    Affected ({len(systs)}): {systs}")

    return gaps


# ─────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Systematic propagation validator for LFV top analysis",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("--indir",   required=True,
                        help="Directory with process output ROOT files")
    parser.add_argument("--sample",  required=True,
                        help="Dataset name substring (e.g. TTto2L2Nu)")
    parser.add_argument("--year",    required=True,
                        help="Year (e.g. 2022, 2022EE, 2023, 2024)")
    parser.add_argument("--ch",      default="muon",
                        choices=["muon","electron"])
    parser.add_argument("--nevents", type=int, default=10000,
                        help="Max events per file (0 = all)")
    parser.add_argument("--outfile", default="",
                        help="Optional path for JSON report")
    args = parser.parse_args()

    print(f"[INFO] Discovering files in: {args.indir}")
    print(f"[INFO] Sample: '{args.sample}'  Year: {args.year}  Channel: {args.ch}")

    sample_files = discover_sample_files(args.indir, args.sample)
    nominal_file = sample_files.get("nominal", "")

    print(f"[INFO] Files found: {len(sample_files)}  "
          f"(nominal: {'YES' if nominal_file else 'NO'})")
    print(f"[INFO] Syst files : {sorted(k for k in sample_files if k != 'nominal')}")

    if not nominal_file:
        sys.exit("[ERROR] No nominal file found. Check --indir and --sample.")

    nevents = args.nevents if args.nevents > 0 else None
    all_results = []

    print("\n[INFO] Running CHECK 1: Branch / file existence ...")
    all_results += check_branch_existence(nominal_file, sample_files, args.ch)

    print("[INFO] Running CHECK 2+3+4: Shape syst propagation ...")
    all_results += check_shape_propagation(nominal_file, sample_files, args.ch, nevents)

    print("[INFO] Running CHECK 5: Weight syst isolation ...")
    all_results += check_weight_isolation(nominal_file, sample_files, nevents)

    print("[INFO] Running CHECK 6: Frozen variables ...")
    all_results += check_frozen_variables(nominal_file, sample_files, nevents)

    print("[INFO] Running CHECK 7: Up/Down directionality ...")
    all_results += check_directionality(nominal_file, sample_files, nevents)

    print()
    render_report(all_results, args.outfile if args.outfile else None)

    gaps = discover_propagation_gaps(all_results)
    if gaps:
        print()
        print("=" * 80)
        print("  PROPAGATION GAP SUMMARY")
        print("=" * 80)
        for g in gaps:
            print(f"  {g}")
        print()

    n_fail = sum(1 for r in all_results if r.status == CheckResult.FAIL)
    sys.exit(1 if n_fail > 0 else 0)


if __name__ == "__main__":
    main()
