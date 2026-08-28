#!/bin/bash
#SBATCH --job-name=cellranger_count
#SBATCH --output=./logs/cellranger_%j.out
#SBATCH --error=./logs/cellranger_%j.err
#SBATCH --cpus-per-task=16
#SBATCH --mem=256G
#SBATCH -p general-cpu
#SBATCH --time=24:00:00

# Set sample ID parameter
ID="SRR9134717"

# Define directory paths
PROJECT_DIR="./Phuc_10x"
FASTQ_DIR="${PROJECT_DIR}/fastq"
REF_DIR="${PROJECT_DIR}/ref"
REF_TAR="${REF_DIR}/refdata-gex-GRCh38-2024-A.tar.gz"
EXTRACTED_REF="${REF_DIR}/refdata-gex-GRCh38-2024-A"

# 1. Unpack reference genome if not already extracted
if [ ! -d "${EXTRACTED_REF}" ]; then
  echo "Extracting reference genome..."
  tar -xvzf "${REF_TAR}" -C "${REF_DIR}"
fi

# 2. Setup directory structure and move/rename FASTQs
mkdir -p "${FASTQ_DIR}"

# 3. Run Cell Ranger inside container via Pyxis/Enroot
srun \
  --container-image='cumulusprod/cellranger:8.0.1' \
  --container-workdir="${PROJECT_DIR}" \
  --container-mounts="${PROJECT_DIR}:${PROJECT_DIR}" \
  bash -c "export PATH=/software/cellranger-8.0.1:\$PATH && cellranger count \
    --id=${ID}_run \
    --fastqs=${FASTQ_DIR} \
    --sample=${ID} \
    --create-bam=true \
    --transcriptome=${EXTRACTED_REF} \
    --localcores=${SLURM_CPUS_PER_TASK}"