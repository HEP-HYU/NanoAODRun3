#!/bin/bash

set -euo pipefail

if [ $# -ne 2 ]; then
    echo "Usage: $0 <check_dir> <version>"
    exit 1
fi

BASE_DIR="$1"
VERSION="$2"

DATASET_DIR="/home/itseyes/github/anti_NanoAODRun3/LFVAnalyzer/data/dataset/${VERSION}"

if [ ! -d "$DATASET_DIR" ]; then
    echo "ERROR: Dataset directory not found:"
    echo "  $DATASET_DIR"
    exit 1
fi

shopt -s nullglob

missing_file="missing_${VERSION}.txt"
> "$missing_file"

for dataset_file in "${DATASET_DIR}"/dataset_*.txt; do

    dataset_name=$(basename "$dataset_file")
    dataset_name=${dataset_name#dataset_}
    dataset_name=${dataset_name%.txt}

    echo "=================================================="
    echo "Checking dataset: ${dataset_name}"
    echo "Source file: ${dataset_file}"

    target_dirs=( "${BASE_DIR}"/*/"${VERSION}"/"${dataset_name}" )

    if [ ${#target_dirs[@]} -eq 0 ]; then
        echo "WARNING: No target directory found:"
        echo "  ${BASE_DIR}/*/${VERSION}/${dataset_name}"
        continue
    fi

    total_count=0
    missing_count=0

    while IFS= read -r root_path || [ -n "$root_path" ]; do

        [[ -z "$root_path" ]] && continue

        ((total_count+=1))

        subdir=$(basename "$(dirname "$root_path")")
        uuid_file=$(basename "$root_path")

        expected_file="${subdir}_${uuid_file}"

        found=false

        for target_dir in "${target_dirs[@]}"; do
            if [ -f "${target_dir}/${expected_file}" ]; then
                found=true
                break
            fi
        done

        if ! $found; then
            echo "$expected_file" >> "$missing_file"
            ((missing_count+=1))
        fi

    done < "$dataset_file"

    echo "Summary:"
    echo "  Total   : ${total_count}"
    echo "  Missing : ${missing_count}"
    echo "  Output  : ${missing_file}"
    echo

done
