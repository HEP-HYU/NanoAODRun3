#!/bin/bash

if $# -eq 0; then
    echo "use: $0 <year>"
    exit 1
fi

year=$1
JSON_FILE="data/dataset/${year}/dataset.json"

DATASETS=$(jq -c '.[] | {DataSetName, name}' "$JSON_FILE")

if ! command -v dasgoclient &> /dev/null; then
    source /cvmfs/cms.cern.ch/cmsset_default.sh
    source /cvmfs/cms.cern.ch/rucio/setup-py3.sh
fi

echo "$DATASETS" | while read -r line; do
    DATASET=$(echo "$line" | jq -r '.DataSetName')
    NAME=$(echo "$line" | jq -r '.name')

    OUTPUT_FILE="data/dataset/${year}/dataset_${NAME}.txt"

    echo "DAS query: $DATASET -> $OUTPUT_FILE"
    dasgoclient --query="file dataset=${DATASET}" > "$OUTPUT_FILE"
    sleep 1
done

