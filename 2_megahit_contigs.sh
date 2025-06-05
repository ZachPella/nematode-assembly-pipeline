#!/bin/bash
#SBATCH --job-name=megahit_assembly
#SBATCH --output=megahit_assembly_%A_%a.out
#SBATCH --error=megahit_assembly_%A_%a.err
#SBATCH --array=1-9
#SBATCH --mem=15G
#SBATCH --time=8:00:00
#SBATCH --cpus-per-task=8

# Change to the working directory
cd /work/fauverlab/zachpella/braker_run/unmapped_reads/bams_surface_sterlizied_namericanus_l3s_for_ZP/fastq_og
module load megahit

# Create output directory for Megahit results
mkdir -p megahit_results

# Create an array of the BAM files (to use the same sample names)
BAM_FILES=($(ls *.sorted.refrename.bam))

# Get the current BAM file based on the array task ID
INDEX=$((SLURM_ARRAY_TASK_ID-1))
CURRENT_BAM=${BAM_FILES[$INDEX]}

# Extract the sample name without the extension
SAMPLE_NAME=$(basename "${CURRENT_BAM}" .sorted.refrename.bam)
echo "Assembling unmapped reads for ${SAMPLE_NAME}..."

# Run Megahit assembly using ALL unmapped reads instead of bacterial-only reads
megahit -1 fastq_for_assembly/${SAMPLE_NAME}.unmapped.R1.fastq \
        -2 fastq_for_assembly/${SAMPLE_NAME}.unmapped.R2.fastq \
        --min-contig-len 500 \
        --k-min 21 \
        --k-max 121 \
        --k-step 20 \
        --presets meta-sensitive \
        -m 14000000000 \
        -t $SLURM_CPUS_PER_TASK \
        -o megahit_results/${SAMPLE_NAME}

# Copy the final contigs to a central location for easier access
mkdir -p all_contigs
cp megahit_results/${SAMPLE_NAME}/final.contigs.fa all_contigs/${SAMPLE_NAME}.contigs.fa

echo "Assembly of unmapped reads completed for ${SAMPLE_NAME}"
