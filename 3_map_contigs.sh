#!/bin/bash
#SBATCH --job-name=contig_mapping
#SBATCH --output=contig_mapping_%A_%a.out
#SBATCH --error=contig_mapping_%A_%a.err
#SBATCH --array=1-9
#SBATCH --mem=16G
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=4

# Load necessary modules
module load bwa
module load samtools

# Change to the working directory
cd /work/fauverlab/zachpella/braker_run/unmapped_reads/bams_surface_sterlizied_namericanus_l3s_for_ZP/fastq_og

# Create output directory
MAPPING_DIR="contig_mapping"
mkdir -p ${MAPPING_DIR}

# Create an array of the sample directories in megahit_results
SAMPLES=($(ls -d megahit_results/Na*))

# Get the current sample
INDEX=$((SLURM_ARRAY_TASK_ID-1))
SAMPLE_PATH=${SAMPLES[$INDEX]}
SAMPLE_NAME=$(basename "${SAMPLE_PATH}")

echo "Processing sample: ${SAMPLE_NAME}"

# Path to assembled contigs
CONTIGS="${SAMPLE_PATH}/final.contigs.fa"

# Check if contigs file exists
if [ ! -f "$CONTIGS" ]; then
    echo "ERROR: Contigs file not found: ${CONTIGS}"
    echo "Has the Megahit assembly finished for this sample?"
    exit 1
fi

# Print contig file info
echo "Contig file details:"
ls -la ${CONTIGS}

# Directory for reference genomes - UPDATED PATH
REF_DIR="/work/fauverlab/zachpella/braker_run/unmapped_reads/bams_surface_sterlizied_namericanus_l3s_for_ZP/fastq_og/reference_genomes"

# List of reference genomes to map against
declare -a REFERENCES=(
  "S_pavanii.fna"
  "S_maltophilia.fna"
  "S_rhizophila.fna"
  "S_sp_610A2.fna"
  "S_acidaminiphila.fna"
  "P_rhodesiae.fna"
  "X_campestris.fna"
)

# Check if reference files exist and are indexed
for REF in "${REFERENCES[@]}"; do
  REF_PATH="${REF_DIR}/${REF}"
  if [ ! -f "${REF_PATH}" ]; then
    echo "ERROR: Reference file not found: ${REF_PATH}"
    exit 1
  fi
  if [ ! -f "${REF_PATH}.sa" ]; then
    echo "ERROR: BWA index not found for: ${REF_PATH}"
    echo "Please make sure the reference is indexed with BWA"
    exit 1
  fi
done

# Map contigs to each reference genome
for REF in "${REFERENCES[@]}"; do
  REF_NAME=$(basename ${REF} .fna)
  REF_PATH="${REF_DIR}/${REF}"

  echo "Mapping contigs from ${SAMPLE_NAME} against ${REF_NAME}..."
  echo "Using reference file: ${REF_PATH}"

  # Create reference-specific output directory
  OUT_DIR="${MAPPING_DIR}/${SAMPLE_NAME}_${REF_NAME}"
  mkdir -p ${OUT_DIR}

  # Map contigs using BWA with explicit paths
  echo "Running: bwa mem -t ${SLURM_CPUS_PER_TASK} ${REF_PATH} ${CONTIGS}"
  bwa mem -t ${SLURM_CPUS_PER_TASK} ${REF_PATH} ${CONTIGS} > ${OUT_DIR}/temp.sam

  # Check if SAM file was created successfully
  if [ ! -s "${OUT_DIR}/temp.sam" ]; then
    echo "ERROR: BWA mapping failed. SAM file is empty or not created."
    continue
  fi

  # Convert SAM to BAM, sort, and index
  echo "Converting SAM to BAM and sorting..."
  samtools view -bS ${OUT_DIR}/temp.sam | \
    samtools sort -o "${OUT_DIR}/${SAMPLE_NAME}_${REF_NAME}.sorted.bam" -

  # Remove temporary SAM file
  rm ${OUT_DIR}/temp.sam

  # Index BAM file
  echo "Indexing BAM file..."
  samtools index "${OUT_DIR}/${SAMPLE_NAME}_${REF_NAME}.sorted.bam"

  # Generate mapping statistics
  echo "Generating mapping statistics..."
  samtools flagstat "${OUT_DIR}/${SAMPLE_NAME}_${REF_NAME}.sorted.bam" > "${OUT_DIR}/${SAMPLE_NAME}_${REF_NAME}.flagstat.txt"

  # Calculate coverage
  echo "Calculating coverage..."
  samtools depth "${OUT_DIR}/${SAMPLE_NAME}_${REF_NAME}.sorted.bam" | \
    awk '{sum+=$3; bases++} END {if(bases>0) print "Average depth: "sum/bases; else print "Average depth: 0"}' > "${OUT_DIR}/${SAMPLE_NAME}_${REF_NAME}.depth.txt"

  # Calculate coverage breadth
  genome_size=$(grep -v "^>" ${REF_PATH} | tr -d '\n' | wc -c)
  covered_bases=$(samtools depth "${OUT_DIR}/${SAMPLE_NAME}_${REF_NAME}.sorted.bam" | wc -l)
  coverage_pct=$(echo "scale=2; 100 * ${covered_bases} / ${genome_size}" | bc)
  echo "Coverage breadth: ${coverage_pct}%" >> "${OUT_DIR}/${SAMPLE_NAME}_${REF_NAME}.depth.txt"

  echo "Completed mapping contigs from ${SAMPLE_NAME} against ${REF_NAME}"
done

echo "All mappings completed for sample ${SAMPLE_NAME}"
