#!/bin/bash
#SBATCH --job-name=comprehensive_stats
#SBATCH --output=comprehensive_stats.out
#SBATCH --error=comprehensive_stats.err
#SBATCH --mem=8G
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=1

# Load required modules
module load samtools

# Base directory
BASE_DIR="/work/fauverlab/zachpella/braker_run/unmapped_reads/bams_surface_sterlizied_namericanus_l3s_for_ZP/fastq_og"
cd ${BASE_DIR}

# Output files
ASSEMBLY_STATS="detailed_assembly_stats.tsv"
MAPPING_STATS="detailed_mapping_stats.tsv"

echo "=== COMPREHENSIVE ASSEMBLY AND MAPPING STATISTICS ==="

# Function to calculate N50
calculate_n50() {
    local fasta_file=$1
    grep -v "^>" ${fasta_file} | tr -d '\n' | fold -w1 | sort | uniq -c | \
    awk '{lengths[NR]=$2; total+=$2} END {
        # Sort lengths in descending order
        PROCINFO["sorted_in"] = "@val_num_desc"
        for (i in lengths) sorted_lengths[++j] = lengths[i]

        target = total * 0.5
        cumulative = 0
        for (i=1; i<=j; i++) {
            cumulative += sorted_lengths[i]
            if (cumulative >= target) {
                print sorted_lengths[i]
                break
            }
        }
    }'
}

# Function to get basic assembly stats
get_assembly_stats() {
    local fasta_file=$1
    local sample_name=$2

    if [ ! -f "$fasta_file" ]; then
        echo "WARNING: Assembly file not found: $fasta_file"
        return
    fi

    # Basic statistics
    total_contigs=$(grep -c "^>" ${fasta_file})
    total_size=$(grep -v "^>" ${fasta_file} | tr -d '\n' | wc -c)
    longest_contig=$(grep -v "^>" ${fasta_file} | awk '{print length}' | sort -nr | head -1)

    # Contigs by size threshold
    contigs_1kb=$(grep -v "^>" ${fasta_file} | awk '{if(length($0) >= 1000) count++} END {print count+0}')
    contigs_5kb=$(grep -v "^>" ${fasta_file} | awk '{if(length($0) >= 5000) count++} END {print count+0}')
    contigs_10kb=$(grep -v "^>" ${fasta_file} | awk '{if(length($0) >= 10000) count++} END {print count+0}')

    # Calculate N50 (simplified)
    n50=$(grep -v "^>" ${fasta_file} | awk '{print length}' | sort -nr | awk '
        {lengths[NR] = $1; total += $1}
        END {
            target = total * 0.5
            cumulative = 0
            for (i=1; i<=NR; i++) {
                cumulative += lengths[i]
                if (cumulative >= target) {
                    print lengths[i]
                    break
                }
            }
        }')

    # GC content
    gc_content=$(grep -v "^>" ${fasta_file} | tr -d '\n' | tr 'ATGCatgc' '0011' | \
                 awk '{gc=gsub(/1/,""); total=length($0); if(total>0) print gc/total*100; else print 0}')

    echo -e "${sample_name}\t${total_contigs}\t${total_size}\t${longest_contig}\t${n50}\t${contigs_1kb}\t${contigs_5kb}\t${contigs_10kb}\t${gc_content}"
}

# Create assembly statistics header
echo -e "Sample\tTotal_Contigs\tTotal_Size_bp\tLongest_Contig_bp\tN50_bp\tContigs_≥1kb\tContigs_≥5kb\tContigs_≥10kb\tGC_Content_%" > ${ASSEMBLY_STATS}

echo "Calculating assembly statistics..."

# Process each sample's assembly
for sample_dir in megahit_results/Na*; do
    if [ -d "$sample_dir" ]; then
        sample_name=$(basename "$sample_dir")
        contigs_file="${sample_dir}/final.contigs.fa"
        echo "Processing ${sample_name}..."
        get_assembly_stats "$contigs_file" "$sample_name" >> ${ASSEMBLY_STATS}
    fi
done

# Create mapping statistics header
echo -e "Sample\tReference\tTotal_Contigs_Mapped\tMapped_Contigs\tMapping_Rate_%\tTotal_Bases_Mapped\tCoverage_Breadth_%\tAverage_Depth\tGenome_Size_bp" > ${MAPPING_STATS}

echo "Calculating detailed mapping statistics..."

# Process mapping results
for mapping_dir in contig_mapping/Na*; do
    if [ -d "$mapping_dir" ]; then
        sample_name=$(basename "$mapping_dir" | cut -d'_' -f1-4)  # Extract sample name
        reference=$(basename "$mapping_dir" | cut -d'_' -f5-)     # Extract reference name

        bam_file="${mapping_dir}/${sample_name}_${reference}.sorted.bam"

        if [ -f "$bam_file" ]; then
            echo "Processing mapping: ${sample_name} -> ${reference}"

            # Get flagstat info
            flagstat_output=$(samtools flagstat "$bam_file")
            total_reads=$(echo "$flagstat_output" | head -1 | awk '{print $1}')
            mapped_reads=$(echo "$flagstat_output" | grep "mapped (" | head -1 | awk '{print $1}')

            # Calculate mapping rate
            if [ "$total_reads" -gt 0 ]; then
                mapping_rate=$(echo "scale=2; $mapped_reads * 100 / $total_reads" | bc)
            else
                mapping_rate=0
            fi

            # Get reference genome size
            ref_file="reference_genomes/${reference}.fna"
            if [ -f "$ref_file" ]; then
                genome_size=$(grep -v "^>" "$ref_file" | tr -d '\n' | wc -c)
            else
                genome_size="Unknown"
            fi

            # Calculate coverage breadth and depth
            if [ -f "${mapping_dir}/${sample_name}_${reference}.depth.txt" ]; then
                coverage_info=$(cat "${mapping_dir}/${sample_name}_${reference}.depth.txt")
                avg_depth=$(echo "$coverage_info" | grep "Average depth" | awk '{print $3}')
                coverage_breadth=$(echo "$coverage_info" | grep "Coverage breadth" | awk '{print $3}' | tr -d '%')
            else
                # Calculate on the fly
                covered_bases=$(samtools depth "$bam_file" | wc -l)
                if [ "$genome_size" != "Unknown" ] && [ "$genome_size" -gt 0 ]; then
                    coverage_breadth=$(echo "scale=2; $covered_bases * 100 / $genome_size" | bc)
                else
                    coverage_breadth=0
                fi

                avg_depth=$(samtools depth "$bam_file" | awk '{sum+=$3; count++} END {if(count>0) print sum/count; else print 0}')
            fi

            # Count mapped contigs
            mapped_contigs=$(samtools view -F 4 "$bam_file" | cut -f1 | sort -u | wc -l)

            # Total bases mapped
            total_bases_mapped=$(samtools view -F 4 "$bam_file" | awk '{sum+=length($10)} END {print sum+0}')

            echo -e "${sample_name}\t${reference}\t${total_reads}\t${mapped_contigs}\t${mapping_rate}\t${total_bases_mapped}\t${coverage_breadth}\t${avg_depth}\t${genome_size}" >> ${MAPPING_STATS}
        fi
    fi
done

# Summary report
echo ""
echo "=== SUMMARY REPORT ==="
echo "Assembly statistics saved to: ${ASSEMBLY_STATS}"
echo "Mapping statistics saved to: ${MAPPING_STATS}"
echo ""

echo "Quick assembly summary:"
echo "Sample\tTotal_Size_bp\tContigs\tN50\tLongest"
tail -n +2 ${ASSEMBLY_STATS} | awk -F'\t' '{printf "%-20s\t%8s\t%6s\t%6s\t%8s\n", $1, $3, $2, $5, $4}'

echo ""
echo "Best mapping hits (>10% coverage breadth):"
echo "Sample\tReference\tCoverage_Breadth_%\tAvg_Depth"
tail -n +2 ${MAPPING_STATS} | awk -F'\t' '{if($7 > 10) printf "%-20s\t%-15s\t%8s\t%8s\n", $1, $2, $7, $8}'

echo ""
echo "Analysis complete! Check the detailed TSV files for full statistics."
