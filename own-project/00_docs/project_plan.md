# Project Plan

## Title
Integrated CLIP-seq and Ribosome Profiling Analysis of LIN28A-Related Translational Regulation

## Biological Question
Which functional gene groups are preferentially LIN28A-bound and show translational changes after `Lin28a` knockdown?

## Hypothesis
Genes with stronger LIN28A-associated CLIP signal are enriched in specific biological processes and show directional translation changes under knockdown.

## Inputs
- `01_data/raw/CLIP-35L33G.bam`
- `01_data/raw/RPF-siLuc.bam`
- `01_data/raw/RPF-siLin28a.bam`
- `01_data/raw/RNA-siLuc.bam`
- `01_data/raw/RNA-siLin28a.bam`
- `01_data/raw/gencode.gtf`
- Reused class-generated `featureCounts` table: `../../binfo1-work/read-counts.txt`

## Workflow
1. Build annotation and gene regions (already completed).
2. Build per-assay count tables from the integrated count matrix.
3. Construct gene-level matrix with CLIP/RPF/RNA counts.
4. Compute CLIP enrichment and TE change metrics.
5. Define quadrant groups from thresholds (`x`, `y`).
6. Run enrichment summary and assemble bubble input table.
7. Export final bubble plot figure (`png`, `pdf`).

## Key Definitions
- `x = clip_enrichment_log2`: RNA-normalized CLIP enrichment on log2 scale.
- `y = te_change_log2 = log2((RPF_siLin28a/RNA_siLin28a)/(RPF_siLuc/RNA_siLuc))`.
- Quadrants:
  - Q1: `x > 0.5`, `y > 0.3`
  - Q2: `x <= -0.5`, `y > 0.3`
  - Q3: `x <= -0.5`, `y <= -0.3`
  - Q4: `x > 0.5`, `y <= -0.3`

## Outputs
- `03_results/counts/*.txt`
- `03_results/metrics/gene_matrix.tsv`
- `03_results/metrics/gene_metrics.tsv`
- `03_results/enrichment/go_enrichment.tsv`
- `03_results/enrichment/go_bubble_input.tsv`
- `03_results/figures/figure_go_clip_ribo_bubble.png`
- `03_results/figures/figure_go_clip_ribo_bubble.pdf`

## Reproducibility
Run scripts in order:
- `02_scripts/02_count_reads.sh`
- `02_scripts/03_build_gene_matrix.py`
- `02_scripts/04_compute_metrics.py`
- `02_scripts/05_go_enrichment.py`
- `02_scripts/06_plot_bubble.py`
