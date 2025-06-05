# Nematode Assembly Pipeline

A comprehensive metagenome assembly pipeline for *de novo* assembly and quality assessment of nematode-associated microbial communities from unmapped sequencing reads.

## Overview

This pipeline extracts unmapped reads from host genomic data, performs metagenome assembly using MEGAHIT, maps assemblies to reference genomes, and provides comprehensive statistical analysis of assembly quality and reference genome coverage. The workflow is optimized for characterizing microbial communities associated with nematodes.

## Pipeline Steps

### 1. Extract Unmapped Reads (`1_run_get_unmapped.sh`)
- Extracts paired-end reads where both mates are unmapped from BAM files
- Converts to FASTQ format for metagenome assembly
- **Input**: `*.sorted.refrename.bam`
- **Output**: `fastq_for_assembly/*.unmapped.R{1,2}.fastq`

### 2. Metagenome Assembly (`2_megahit_contigs.sh`)
- *De novo* assembly using MEGAHIT optimized for metagenomes
- Handles variable coverage depths typical of microbial communities
- **Output**: `megahit_results/*/final.contigs.fa`

### 3. Reference Genome Mapping (`3_map_contigs.sh`)
- Maps assembled contigs to reference bacterial genomes
- Identifies potential matches to known bacterial species
- **Output**: `contig_mapping/*/sample_reference.sorted.bam`

### 4. Comprehensive Statistics (`calc_stats.sh`)
- Calculates detailed assembly quality metrics
- Analyzes mapping coverage to reference genomes
- Generates summary reports and visualizations
- **Output**: 
  - `detailed_assembly_stats.tsv`
  - `detailed_mapping_stats.tsv`
  - Summary reports

## Requirements

### Software Dependencies
- **SLURM** (job scheduler)
- **samtools** (≥1.9)
- **MEGAHIT** (≥1.2.9)
- **BWA-MEM** or **Bowtie2**
- **bc** (basic calculator for statistics)

### System Requirements
- HPC cluster with SLURM
- Memory: 8-32GB depending on step
- CPUs: 1-16 cores per job
- Storage: ~20GB per sample for assemblies and mapping

## Usage

### Sequential Workflow
```bash
# 1. Extract unmapped reads (array job for multiple samples)
sbatch 1_run_get_unmapped.sh

# 2. Assemble contigs with MEGAHIT
sbatch 2_megahit_contigs.sh

# 3. Map contigs to reference genomes
sbatch 3_map_contigs.sh

# 4. Calculate comprehensive statistics
sbatch calc_stats.sh
```

### Input Data Structure
```
project_directory/
├── *.sorted.refrename.bam    # Input BAM files
├── reference_genomes/        # Reference bacterial genomes
│   └── *.fna                # FASTA reference files
├── 1_run_get_unmapped.sh     # Pipeline scripts
├── 2_megahit_contigs.sh
├── 3_map_contigs.sh
└── calc_stats.sh
```

### Output Structure
```
project_directory/
├── fastq_for_assembly/          # Unmapped FASTQ files
├── megahit_results/             # Assembly outputs
│   └── Sample_Name/
│       ├── final.contigs.fa     # Assembled contigs
│       └── log                  # Assembly logs
├── contig_mapping/              # Mapping results
│   └── Sample_Reference/
│       ├── *.sorted.bam         # Mapped reads
│       └── *.depth.txt          # Coverage depth
├── detailed_assembly_stats.tsv  # Assembly quality metrics
├── detailed_mapping_stats.tsv   # Reference mapping analysis
└── comprehensive_stats.out      # Summary report
```

## Key Metrics

### Assembly Quality Assessment
- **Total Contigs**: Number of assembled sequences
- **Total Assembly Size**: Combined length of all contigs
- **N50**: Median contig length (quality indicator)
- **Longest Contig**: Maximum single contig length
- **Size Distribution**: Contigs ≥1kb, ≥5kb, ≥10kb
- **GC Content**: Overall nucleotide composition

### Reference Mapping Analysis
- **Mapping Rate**: Percentage of reads mapping to references
- **Coverage Breadth**: Percentage of reference genome covered
- **Average Depth**: Mean sequencing depth across covered regions
- **Genome Size**: Reference genome statistics
- **Best Matches**: Identification of closest bacterial relatives

## Expected Results

### Assembly Metrics
- **N50 Range**: 1-20kb typical for nematode microbiomes
- **Contig Count**: 100-10,000 depending on diversity
- **Total Size**: 1-100Mb assembly size range
- **Quality Threshold**: >1000 contigs ≥1kb indicates good assembly

### Reference Mapping
- **Coverage Identification**: >10% breadth indicates potential matches
- **Depth Analysis**: >5× average depth suggests reliable detection
- **Taxonomic Insights**: Identifies dominant bacterial associates

## Quality Control

### Assembly Quality Indicators
- **High Quality**: N50 >5kb, >50% contigs ≥1kb
- **Medium Quality**: N50 1-5kb, >25% contigs ≥1kb  
- **Low Quality**: N50 <1kb, <25% contigs ≥1kb

### Mapping Quality Thresholds
- **Strong Match**: >50% coverage breadth, >10× depth
- **Moderate Match**: 10-50% coverage breadth, >5× depth
- **Weak Match**: <10% coverage breadth or <5× depth

## Customization

### Assembly Parameters
- **K-mer Sizes**: Modify MEGAHIT k-mer range for different data types
- **Memory Allocation**: Adjust based on sample complexity
- **Minimum Length**: Change contig length thresholds

### Reference Database
- **Custom References**: Add specific bacterial genomes of interest
- **Database Updates**: Include latest genomic sequences
- **Filtering**: Focus on particular taxonomic groups

### Statistical Analysis
- **Coverage Thresholds**: Adjust minimum coverage requirements
- **Quality Metrics**: Modify assembly quality standards
- **Reporting**: Customize summary output formats

## Troubleshooting

### Common Issues
1. **Assembly Failures**: Increase memory allocation or reduce k-mer size
2. **Low Mapping Rates**: Check reference genome quality and relevance
3. **Poor Assembly Quality**: Verify input read quality and depth
4. **Memory Errors**: Use high-memory nodes for large datasets

### Performance Optimization
- **Parallel Processing**: Run multiple samples simultaneously
- **Resource Allocation**: Match CPU/memory to data complexity
- **Storage Management**: Clean intermediate files to save space

## Interpretation Guide

### Assembly Success Indicators
- **N50 >5kb**: Excellent assembly quality
- **Many long contigs**: Good coverage and complexity resolution
- **Reasonable total size**: Matches expected microbial diversity

### Reference Mapping Insights
- **High coverage matches**: Likely presence of related bacteria
- **Multiple weak matches**: Diverse microbial community
- **Novel sequences**: Potential undescribed bacterial associates

## Citation

If you use this pipeline, please cite:
- **MEGAHIT**: Li et al. (2015) Bioinformatics 31(10):1674-1676
- **BWA**: Li & Durbin (2009) Bioinformatics 25(14):1754-1760
- **SAMtools**: Li et al. (2009) Bioinformatics 25(16):2078-2079

## Best Practices

### Data Preparation
- **Quality Control**: Ensure high-quality unmapped reads
- **Read Depth**: >1M read pairs recommended per sample
- **Reference Selection**: Use high-quality, complete genomes

### Analysis Strategy
- **Multiple References**: Include diverse bacterial representatives
- **Validation**: Cross-reference assembly and mapping results
- **Documentation**: Record assembly parameters and versions

## License

MIT License - free to use and modify for research purposes.

---

**Pipeline for nematode-associated microbial assembly and analysis**  
*From unmapped reads → Assembled contigs → Reference mapping → Quality metrics*
