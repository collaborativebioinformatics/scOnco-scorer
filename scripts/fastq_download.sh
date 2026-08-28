#!/bin/bash

INPUT_FILE="filtered_relevant_runs.tsv"
LOG_DIR="slurm_logs"

mkdir -p "$LOG_DIR"

# Extract run_accession column (Column 1), skipping header
tail -n +2 "$INPUT_FILE" | cut -f1 | while read -r SRR; do
    # Skip empty lines
    [ -z "$SRR" ] && continue

    echo "Submitting Slurm job for accession: $SRR"

    sbatch <<EOF
#!/bin/bash
#SBATCH --job-name=sra_$SRR
#SBATCH --output=$LOG_DIR/${SRR}_%j.log
#SBATCH --error=$LOG_DIR/${SRR}_%j.err
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH -p general-cpu
#SBATCH --time=04:00:00

echo "Starting download for $SRR on \$(hostname) at \$(date)"

echo 'export PATH=$PATH:sratoolkit/bin' >> ~/.bashrc
source ~/.bashrc

# 1. Download .sra file archive via prefetch
echo "Running prefetch..."
prefetch "$SRR" --output-directory . --max-size 50g

# 2. Extract paired/single FASTQs with multi-threading
echo "Running fasterq-dump..."
fasterq-dump "$SRR" \
    --split-files \
    --threads \$SLURM_CPUS_PER_TASK \
    --mem 14G \
    --temp \${TMPDIR:-/tmp} \
    --outdir .

# 3. Clean up the intermediate .sra folder to free disk space
if [ -d "$SRR" ]; then
    rm -rf "$SRR"
fi

echo "Finished $SRR at \$(date)"
EOF

done