import os
import argparse

parser = argparse.ArgumentParser(description="Find problematic ROOT files")
parser.add_argument("input_path", help="Base directory to scan")

args = parser.parse_args()

LOG_DIR = args.input_path+"/log"

ERROR_KEYWORDS = [
    "runtime_error",
    "fatal",
    "fault",
    "Traceback",
    "ERROR",
    "error"
]

SUCCESS_KEYWORD = "Job Done"
NOEVENT_KEYWORD = "There is NO EVENT to process, ending the processing!!"


def analyze_logs(base_dir):
    error_files = []
    error_but_done_files = []
    no_done_files = []
    count = 0

    for root, dirs, files in os.walk(base_dir):
        for file in files:
            if file.endswith(".log"):
                filepath = os.path.join(root, file)

                try:
                    with open(filepath, "r", errors="ignore") as f:
                        content = f.read()

                    # 에러 키워드 포함 여부
                    has_error = any(k in content for k in ERROR_KEYWORDS)

                    if has_error and SUCCESS_KEYWORD not in content:
                        error_files.append(filepath)

                    # Job Done 포함 여부
                    if has_error and SUCCESS_KEYWORD in content:
                        error_but_done_files.append(filepath)

                    if SUCCESS_KEYWORD not in content and NOEVENT_KEYWORD not in content:
                        no_done_files.append(filepath)

                    if not has_error and (SUCCESS_KEYWORD in content or NOEVENT_KEYWORD in content):
                        count = count + 1

                except Exception as e:
                    print(f"Error reading {filepath}: {e}")

    print ("done: ", count)
    return error_files, error_but_done_files, no_done_files


if __name__ == "__main__":
    error_files, error_but_done_files, no_done_files = analyze_logs(LOG_DIR)

    # 저장
    with open("error_logs.txt", "w") as f:
        for ef in error_files:
            f.write(ef + "\n")

    with open("error_but_done_logs.txt", "w") as f:
        for ef in error_but_done_files:
            f.write(ef + "\n")

    with open("not_done_logs.txt", "w") as f:
        for ef in no_done_files:
            f.write(ef + "\n")

    print(f"Total error logs: {len(error_files)}")
    print(f"Error logs with 'Job Done': {len(error_but_done_files)}")
    print(f"Job is not done: {len(no_done_files)}")
