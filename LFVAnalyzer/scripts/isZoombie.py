import os
import argparse
import ROOT

# 필수로 존재해야 하는 object 이름들
REQUIRED_KEYS = [
    "Events",
    "hcounter",
    "hcounter_S1",
    "LHEPdfWeightSum",
    "PSWeightSum",
    "ScaleWeightSum"
]


def has_required_keys(tfile, is_skim):
    """필수 key 존재 여부 체크"""
    keys = [key.GetName() for key in tfile.GetListOfKeys()]

    if not is_skim:
        REQUIRED_KEYS = ["Events", "h_tau1_pt_S5"]

    for req in REQUIRED_KEYS:
        if req not in keys:
            return False
    return True


def find_problematic_root_files(base_dir=".", is_skim=False):
    bad_files = []

    for root, dirs, files in os.walk(base_dir):
        for file in files:
            if file.endswith(".root"):
                filepath = os.path.join(root, file)

                try:
                    f = ROOT.TFile.Open(filepath)

                    # 1. 파일 열기 실패
                    if not f:
                        bad_files.append(filepath)
                        continue

                    # 2. Zombie 파일
                    if f.IsZombie():
                        bad_files.append(filepath)
                        f.Close()
                        continue

                    # 3. 필수 key 누락
                    if not has_required_keys(f, is_skim):
                        bad_files.append(filepath)

                    f.Close()

                except Exception as e:
                    print(f"Error opening {filepath}: {e}")
                    bad_files.append(filepath)

    return bad_files


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Find problematic ROOT files")
    parser.add_argument("input_path", help="Base directory to scan")
    parser.add_argument("--isSkim", action="store_true",
                        help="Enable skim mode (check full key set)")

    args = parser.parse_args()

    output_file = "zombie.txt"

    bad_files = find_problematic_root_files(args.input_path, args.isSkim)

    with open(output_file, "w") as f:
        for bf in bad_files:
            f.write(bf + "\n")

    print(f"Saved {len(bad_files)} problematic files to {output_file}")
