#!/bin/bash
set -e

cd ~/Projects/Cancer_pipeline/Scrna
mkdir -p logs

tail -n +2 metadata/sample_metadata.tsv | cut -f2 | while read SAMPLE
do
    echo "Running $SAMPLE"
    python scripts/run_one_sample.py --sample "$SAMPLE" \
        > "logs/${SAMPLE}_pipeline.log" 2>&1
    echo "Finished $SAMPLE"
done
