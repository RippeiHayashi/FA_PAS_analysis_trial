# FA_PAS_analysis

A lightweight pipeline for calling polyadenylation sites (PAS) from **Oxford Nanopore direct RNA-seq** data and identifying **formaldehyde-enriched PAS** at the site level.

This workflow was developed for **nuclear direct RNA sequencing of untreated and formaldehyde-treated HeLa cells**. It extracts read 3′ ends from genome-aligned BAM files, builds a stringent PAS catalog, assigns reads back to PAS, generates PAS count matrices, performs **PAS-level differential abundance analysis with edgeR**, and produces browser tracks for visualization.

## Overview

The pipeline is designed for Nanopore direct RNA sequencing data after:

1. **Dorado basecalling** with poly(A) estimation enabled
2. **minimap2 genome alignment**
3. filtering to **primary mapped alignments**

The core PAS workflow then consists of:

1. extracting one genomic 3′ end per aligned read
2. clustering read ends into candidate PAS
3. assigning reads back to PAS
4. generating PAS count matrices
5. testing PAS-level differential abundance with **edgeR**
6. exporting browser tracks and BED files for visualization

## Included scripts

### Python
- `py/extract_read_ends.py`  
  Extract read 3′ end coordinates from primary BAM alignments.

- `py/cluster_pas_strict.py`  
  Build a stringent PAS catalog from pooled read-end coordinates.

- `py/assign_reads_to_pas.py`  
  Assign reads back to the nearest PAS summit.

- `py/make_pas_matrices.py`  
  Build PAS count matrices for downstream analysis.

### R
- `r/run_edger_pas.R`  
  Perform PAS-level differential abundance analysis using edgeR.

### Config
- `config/samples.tsv`  
  Sample sheet describing sample IDs and conditions.

## Input requirements

### Required inputs
- one **coordinate-sorted BAM** per biological replicate
- BAM index (`.bai`)
- a sample sheet with:
  - `sample_id`
  - `condition`
- reference genome FASTA
- chromosome sizes file for bigWig generation

### Assumptions
This pipeline assumes that BAM files:
- were generated from **Nanopore direct RNA-seq**
- contain **primary genome alignments**
- were filtered to remove unmapped, secondary, and supplementary/chimeric alignments

## Software requirements

### Basecalling and alignment
- Dorado v1.0.0
- minimap2 v2.28
- samtools
- bedtools
- UCSC `bedGraphToBigWig`

### Python
Recommended Python packages:
- `pysam`
- standard library modules only for the remaining scripts

### R
Required R packages:
- `edgeR`
- `readr`
- `dplyr`
- `tibble`

## Example workflow

## 1. Dorado basecalling with poly(A) estimation

```bash
export PATH="/path/to/dorado/bin:${PATH}"

dorado basecaller sup,m5C_2OmeC,inosine_m6A_2OmeA,pseU_2OmeU,2OmeG   ${PBS_JOBFS}/pod5_dir/   --estimate-poly-a   -r -b 416 -c 9216   --models-directory /path/to/dorado/bin   > ${PBS_JOBFS}/bam_dir/${pod5_name}.bam

cp ${PBS_JOBFS}/bam_dir/${pod5_name}.bam ${output_directory}/
```

Extract FASTQ while preserving Dorado tags:

```bash
TAGS=$(samtools view ${output_directory}/${pod5_name}.bam | cut -f 12- | tr '\t' '\n' | cut -d ':' -f 1 | sort | uniq | tr '\n' ',' | sed 's/,$//')

samtools fastq -T $TAGS ${output_directory}/${pod5_name}.bam   > ${output_directory}/${pod5_name}_basecall_dorado-v1.0.0.fastq
```

## 2. Genome alignment with minimap2

```bash
module load minimap2/2.28

genome="/path/to/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa"

cat ${output_directory}/*_basecall_dorado-v1.0.0.fastq   > ${output_directory}/${library_name}_merged.fastq

merged_fastq="${output_directory}/${library_name}_merged.fastq"

minimap2 -a -x splice -y -t 24 "${genome}" "${merged_fastq}"   | samtools view -@ 24 -b   > ${output_directory}/${library_name}_merged_total_genome.bam

samtools view -F 2308 -b -@ 24   "${output_directory}/${library_name}_merged_total_genome.bam"   | samtools sort -@ 24   > "${output_directory}/${library_name}_merged_primary_genome.bam"

samtools index "${output_directory}/${library_name}_merged_primary_genome.bam"
```

## 3. PAS calling and differential abundance analysis

### Example sample setup

```bash
analysis="/path/to/FA_PAS_analysis"
scripts="/path/to/FA_PAS_analysis"

mkdir -p ${analysis}/results/read_ends
mkdir -p ${analysis}/results/pas_catalog
mkdir -p ${analysis}/results/counts
mkdir -p ${analysis}/results/stats
```

### 3.1 Extract read ends

```bash
python3 ${scripts}/py/extract_read_ends.py   --bam /path/to/ctrl_1_merged_primary_genome.bam   --sample-id ctrl_1   --condition control   --out-tsv ${analysis}/results/read_ends/ctrl_1.read_ends.tsv

python3 ${scripts}/py/extract_read_ends.py   --bam /path/to/ctrl_2_merged_primary_genome.bam   --sample-id ctrl_2   --condition control   --out-tsv ${analysis}/results/read_ends/ctrl_2.read_ends.tsv

python3 ${scripts}/py/extract_read_ends.py   --bam /path/to/fa_1_merged_primary_genome.bam   --sample-id fa_1   --condition formaldehyde   --out-tsv ${analysis}/results/read_ends/fa_1.read_ends.tsv

python3 ${scripts}/py/extract_read_ends.py   --bam /path/to/fa_2_merged_primary_genome.bam   --sample-id fa_2   --condition formaldehyde   --out-tsv ${analysis}/results/read_ends/fa_2.read_ends.tsv
```

### 3.2 Build PAS catalog

```bash
python3 ${scripts}/py/cluster_pas_strict.py   --inputs ${analysis}/results/read_ends/*.read_ends.tsv   --max-dist 24   --min-pos-reads 2   --min-pos-samples 2   --peak-half-window 12   --min-peak-reads 3   --min-peak-fraction 0.10   --min-final-reads 10   --min-final-samples 2   --out-tsv ${analysis}/results/pas_catalog/pas_catalog.strict.tsv
```

The `cluster_pas_strict.py` step builds a shared PAS catalog from pooled read 3′ ends across all samples. It first filters weak exact positions, then groups nearby positions into coarse clusters, identifies local peaks within those clusters, and finally removes weak PAS with poor pooled support.

Key parameters:

- `--max-dist 24`: group nearby exact positions into the same **coarse PAS neighborhood** if they are within 24 nt
- `--min-pos-reads 2`: keep exact positions with at least 2 pooled reads before clustering
- `--min-pos-samples 2`: or keep positions seen in at least 2 samples, even if total reads are low
- `--peak-half-window 12`: within each coarse cluster, detect local maxima in a ±12 nt neighborhood; very nearby competing peaks are collapsed, while sufficiently separated strong peaks can be retained as separate summits
- `--min-peak-reads 3`: require at least 3 reads at a summit position for it to be considered a local peak
- `--min-peak-fraction 0.10`: summit must represent at least 10% of reads in the coarse cluster
- `--min-final-reads 10`: keep final PAS with at least 10 pooled reads after peak splitting
- `--min-final-samples 2`: keep final PAS supported by at least 2 samples

### 3.3 Assign reads back to PAS

```bash
python3 ${scripts}/py/assign_reads_to_pas.py   --read-ends ${analysis}/results/read_ends/ctrl_1.read_ends.tsv   --pas-catalog ${analysis}/results/pas_catalog/pas_catalog.strict.tsv   --out-tsv ${analysis}/results/counts/ctrl_1.assign.tsv

python3 ${scripts}/py/assign_reads_to_pas.py   --read-ends ${analysis}/results/read_ends/ctrl_2.read_ends.tsv   --pas-catalog ${analysis}/results/pas_catalog/pas_catalog.strict.tsv   --out-tsv ${analysis}/results/counts/ctrl_2.assign.tsv

python3 ${scripts}/py/assign_reads_to_pas.py   --read-ends ${analysis}/results/read_ends/fa_1.read_ends.tsv   --pas-catalog ${analysis}/results/pas_catalog/pas_catalog.strict.tsv   --out-tsv ${analysis}/results/counts/fa_1.assign.tsv

python3 ${scripts}/py/assign_reads_to_pas.py   --read-ends ${analysis}/results/read_ends/fa_2.read_ends.tsv   --pas-catalog ${analysis}/results/pas_catalog/pas_catalog.strict.tsv   --out-tsv ${analysis}/results/counts/fa_2.assign.tsv
```

### 3.4 Build PAS count matrices

```bash
python3 ${scripts}/py/make_pas_matrices.py   --assignments ${analysis}/results/counts/*.assign.tsv   --out-counts ${analysis}/results/counts/pas_counts.tsv   --out-long ${analysis}/results/counts/pas_counts.long.tsv
```

### 3.5 Run PAS-level differential abundance with edgeR

```bash
Rscript ${scripts}/r/run_edger_pas.R   ${analysis}/results/counts/pas_counts.tsv   ${scripts}/config/samples.tsv   ${analysis}/results/stats/pas_edger.tsv
```

## Visualization tracks

### Generate bigWig files of read ends

```bash
for lib in ctrl_1 ctrl_2 fa_1 fa_2; do

  awk 'BEGIN{OFS="\t"} NR>1 {print $4,$5,$6,$1,1,$7}'     ${analysis}/results/read_ends/${lib}.read_ends.tsv     | sort -k1,1 -k2,2n     > ${analysis}/results/read_ends/${lib}.read_ends.bed

  READS=$(wc -l ${analysis}/results/read_ends/${lib}.read_ends.bed | awk '{print $1}')

  bedtools genomecov -strand + -split -bg     -i ${analysis}/results/read_ends/${lib}.read_ends.bed     -g /path/to/chrNameLength.txt     | awk -v READS=${READS} '{print $1,$2,$3,$4/READS*1000000}'     | tr ' ' '\t'     > ${analysis}/results/read_ends/${lib}.read_ends.plus.bg

  bedtools genomecov -strand - -split -bg     -i ${analysis}/results/read_ends/${lib}.read_ends.bed     -g /path/to/chrNameLength.txt     | awk -v READS=${READS} '{print $1,$2,$3,-$4/READS*1000000}'     | tr ' ' '\t'     > ${analysis}/results/read_ends/${lib}.read_ends.minus.bg

  for strand in plus minus; do
    bedGraphToBigWig       ${analysis}/results/read_ends/${lib}.read_ends.${strand}.bg       /path/to/chrNameLength.txt       ${analysis}/results/read_ends/${lib}.read_ends.${strand}.bw
  done
done
```

### Generate PAS BED files

```bash
awk 'BEGIN{OFS="\t"} NR>1 {print "chr"$2,$7,$8,$1,$10,$5}'   ${analysis}/results/pas_catalog/pas_catalog.strict.tsv   > ${analysis}/results/pas_catalog/pas_catalog.strict.bed
```

Extract PAS significantly increased in formaldehyde:

```bash
awk 'BEGIN{FS=OFS="\t"}
NR==FNR {
  if (FNR>1 && $2>0 && $6<0.05) keep[$1]=1
  next
}
($4 in keep)' ${analysis}/results/stats/pas_edger.tsv ${analysis}/results/pas_catalog/pas_catalog.strict.bed > ${analysis}/results/pas_catalog/pas_formaldehyde_up.FDR0.05.bed
```

## Output files

### Main outputs
- `results/read_ends/*.read_ends.tsv`  
  per-read 3′ end coordinates

- `results/pas_catalog/pas_catalog.strict.tsv`  
  stringent PAS catalog

- `results/counts/*.assign.tsv`  
  read-to-PAS assignments

- `results/counts/pas_counts.tsv`  
  PAS-by-sample count matrix

- `results/stats/pas_edger.tsv`  
  edgeR PAS-level differential abundance results

### Browser outputs
- `results/read_ends/*.bw`  
  strand-specific read-end bigWig files

- `results/pas_catalog/pas_catalog.strict.bed`  
  PAS coordinates for genome browser visualization

- `results/pas_catalog/pas_formaldehyde_up.FDR0.05.bed`  
  PAS significantly increased in formaldehyde

## Notes and caveats

- This repository currently implements **PAS calling and PAS-level differential abundance**.
- PAS-level edgeR analysis tests **absolute PAS abundance changes**, not PAS usage within genes.
- Replicates should remain separate throughout the workflow.
- One shallow library can usually be retained if it behaves concordantly with other replicates.
- Chromosome naming must be consistent between BAM, BED, bigWig, and reference files.
- The current public workflow does not yet include:
  - gene annotation of PAS
  - DRIMSeq-based differential PAS usage within genes
  - poly(A) tail length comparison
  - m6A analysis from Dorado modification tags

## Future extensions

Planned additions may include:
- gene-level PAS annotation
- intronic/internal PAS classification
- DRIMSeq-based differential PAS usage
- poly(A)-tail length analysis using Dorado `pt` tags
- m6A analysis on PAS-stratified read sets

## Citation

If you use this workflow, please cite the relevant software:
- Dorado
- minimap2
- edgeR
- bedtools
- samtools
