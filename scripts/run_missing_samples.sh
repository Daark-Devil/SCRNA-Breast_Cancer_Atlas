#!/bin/bash
set -e

cd ~/Projects/Cancer_pipeline/Scrna

find results -name "*_annotated.h5ad" | sed 's|results/||; s|/h5ad/.*||' | sort > completed_samples.txt
cut -f2 metadata/sample_metadata.tsv | tail -n +2 | sort > expected_samples.txt

comm -23 expected_samples.txt completed_samples.txt > missing_samples.txt

echo "Missing samples:"
cat missing_samples.txt

cat missing_samples.txt | while read SAMPLE
do
    echo "Running missing sample: $SAMPLE"
    python scripts/run_one_sample.py --sample "$SAMPLE" \
        > "logs/${SAMPLE}_pipeline.log" 2>&1
    echo "Finished missing sample: $SAMPLE"
done
