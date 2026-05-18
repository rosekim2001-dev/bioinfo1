# Own Bioinformatics Project Skeleton

## 1) Project Title
**Integrated CLIP-seq and Ribosome Profiling Analysis of LIN28A-Related Translational Regulation**

## 2) Core Goal
Reproduce a figure style similar to the reference bubble-scatter plot by integrating:
- CLIP enrichment per gene (x-axis)
- Ribosome density change after `Lin28a` knockdown (y-axis)
- GO term enrichment significance and term size (bubble color/size)

---

## 3) Biological Question and Hypothesis
### Question
Which functional gene groups are preferentially LIN28A-bound and show translational changes after `Lin28a` knockdown?

### Hypothesis
Genes with stronger LIN28A-associated CLIP signal are enriched in specific biological processes, and these processes show directional translation changes (ribosome density up/down) under knockdown.

---

## 4) Data and Paths
Use your existing class dataset.

### Required files
- `binfo1-datapack1/CLIP-35L33G.bam`
- `binfo1-datapack1/RPF-siLuc.bam`
- `binfo1-datapack1/RPF-siLin28a.bam`
- `binfo1-datapack1/RNA-siLuc.bam`
- `binfo1-datapack1/RNA-siLin28a.bam`
- `binfo1-datapack1/gencode.gtf`

### Working directory
- `binfo1-work/own-project/`

---

## 5) Project Folder Structure
```text
own-project/
├── 00_docs/
│   ├── project_plan.md
│   └── figure_notes.md
├── 01_data/
│   ├── raw/                  # symlink or copied input BAM/GTF
│   ├── ref/
│   └── metadata/
├── 02_scripts/
│   ├── 01_prepare_annotation.sh
│   ├── 02_count_reads.sh
│   ├── 03_build_gene_matrix.py
│   ├── 04_compute_metrics.py
│   ├── 05_go_enrichment.py
│   └── 06_plot_bubble.py
├── 03_results/
│   ├── counts/
│   ├── metrics/
│   ├── enrichment/
│   ├── figures/
│   └── tables/
├── 04_notebooks/
│   └── own_analysis.ipynb
├── environment.yml
└── README.md
```

---

## 6) Environment Setup
Create `environment.yml`:
```yaml
name: own-bioinfo
channels:
  - conda-forge
  - bioconda
  - defaults
dependencies:
  - python=3.11
  - samtools
  - bedtools
  - subread
  - pandas
  - numpy
  - scipy
  - statsmodels
  - matplotlib
  - seaborn
  - jupyter
  - pip
  - pip:
      - gseapy
```

Install:
```bash
conda env create -f environment.yml
conda activate own-bioinfo
```

---

## 7) Analysis Strategy

### Step A. Prepare gene annotation
1. Extract `gene_id`, `gene_name`, `gene_type`, genomic intervals from GTF.
2. Build gene-level BED (merged exons or gene span; pick one and keep consistent).

Output:
- `03_results/tables/gene_annotation.tsv`
- `03_results/tables/gene_regions.bed`

### Step B. Count reads per gene for each assay
Count mapped reads overlapping gene regions for:
- CLIP
- RPF (`siLuc`, `siLin28a`)
- RNA (`siLuc`, `siLin28a`)

Recommended command template (example with `featureCounts`):
```bash
featureCounts -a 01_data/raw/gencode.gtf -o 03_results/counts/clip_counts.txt \
  -T 4 -t exon -g gene_id 01_data/raw/CLIP-35L33G.bam
```
Run similarly for each BAM.

Outputs:
- `03_results/counts/clip_counts.txt`
- `03_results/counts/rpf_siLuc_counts.txt`
- `03_results/counts/rpf_siLin28a_counts.txt`
- `03_results/counts/rna_siLuc_counts.txt`
- `03_results/counts/rna_siLin28a_counts.txt`

### Step C. Build integrated gene matrix
Merge counts by `gene_id`:
- `clip_count`
- `rpf_siLuc`, `rpf_siLin28a`
- `rna_siLuc`, `rna_siLin28a`

Add pseudocount (e.g., `+1`) and filter low-expression genes.

Output:
- `03_results/metrics/gene_matrix.tsv`

### Step D. Compute key metrics for figure axes
Define:
- **CLIP enrichment (x):**
  - Option 1 (simple): `log2((clip_count + 1) / median(clip_count + 1))`
  - Option 2 (better): normalize by RNA abundance, then log2 transform.
- **Ribosome density change (y):**
  - `TE_siLuc = (rpf_siLuc + 1) / (rna_siLuc + 1)`
  - `TE_siLin28a = (rpf_siLin28a + 1) / (rna_siLin28a + 1)`
  - `y = log2(TE_siLin28a / TE_siLuc)`

Output:
- `03_results/metrics/gene_metrics.tsv`

### Step E. Split gene sets for enrichment
Define groups for GO analysis, for example:
- `CLIP_enriched`: x > 0.5
- `CLIP_depleted`: x < -0.5
- `TE_up`: y > 0.3
- `TE_down`: y < -0.3
- Combined quadrants:
  - Q1: x high, y high
  - Q2: x low, y high
  - Q3: x low, y low
  - Q4: x high, y low

Run GO enrichment for selected groups.

Output:
- `03_results/enrichment/*.tsv`

### Step F. Build bubble-scatter summary table (GO-level)
For each significant GO term, compute:
- `go_term`
- `n_genes_in_term` (bubble size)
- `x_go`: average/weighted average x across term genes
- `y_go`: average/weighted average y across term genes
- `fdr` (for bubble color)
- optional label score for auto-annotation

Output:
- `03_results/enrichment/go_bubble_input.tsv`

### Step G. Plot final figure
Plot style target:
- Scatter-bubble chart
- X = CLIP enrichment (`log2`)
- Y = ribosome density change (`log2`)
- Bubble size = term size
- Bubble color = `-log10(FDR)` or log-scaled FDR colorbar
- Add labeled top terms with arrows
- Horizontal line at y=0 and vertical line at x=0

Output:
- `03_results/figures/figure_go_clip_ribo_bubble.png`
- `03_results/figures/figure_go_clip_ribo_bubble.pdf`

---

## 8) Notebook Blueprint (`04_notebooks/own_analysis.ipynb`)

### Section 1. Setup and imports
- environment check
- package versions

### Section 2. Load annotation and count files
- parse `featureCounts` outputs
- build merged gene table

### Section 3. Compute x/y metrics
- TE calculation
- CLIP enrichment normalization choice
- save metrics table

### Section 4. Define gene groups
- thresholds
- group size sanity checks

### Section 5. GO enrichment
- run gseapy or your preferred tool
- collect significant terms

### Section 6. Bubble plot assembly
- build GO-level coordinates
- map size and color
- label selected terms

### Section 7. Export
- figures
- final TSV for report

### Section 8. Interpretation notes
- 3–5 key biological takeaways
- method limitations

---

## 9) Script Skeleton Details

### `02_scripts/01_prepare_annotation.sh`
- Input: GTF
- Output: cleaned annotation TSV + BED

### `02_scripts/02_count_reads.sh`
- Input: BAMs + GTF
- Output: count files for each assay
- Include basic QC summary (total mapped, assigned reads)

### `02_scripts/03_build_gene_matrix.py`
- Merge count tables
- Join gene names/types
- Save unified matrix

### `02_scripts/04_compute_metrics.py`
- Compute x/y metrics with pseudocount/filter
- Save gene metrics

### `02_scripts/05_go_enrichment.py`
- Run GO enrichment for selected groups
- Save all and significant-only tables

### `02_scripts/06_plot_bubble.py`
- Generate publication-style bubble plot
- Optional auto-label top N terms by significance

---

## 10) Quality Control Checklist
- BAM index exists for all BAM files (`.bai`)
- Gene counts are non-zero for expected genes
- Replicate absence acknowledged (single sample per condition)
- Threshold sensitivity checked (at least 2 threshold sets)
- Plot axis labels and colorbar units clearly shown

---

## 11) What to Submit
1. `04_notebooks/own_analysis.ipynb`
2. `02_scripts/*.sh` and `02_scripts/*.py`
3. Final figure (`.png` + optional `.pdf`)
4. `03_results/enrichment/go_bubble_input.tsv`
5. Short report (1–2 pages) including:
   - hypothesis
   - methods
   - main figure
   - interpretation
   - limitations and future work

---

## 12) Suggested Timeline (Fast Execution)
1. Day 1: setup, counts, matrix build
2. Day 2: metric definition, GO enrichment
3. Day 3: final plot polishing, write-up

---

## 13) Notes on Similarity vs Originality
This project is aligned with your class workflow and reference style, but is still your own work because:
- You integrate CLIP + RPF + RNA at gene/GO level.
- You define your own metric thresholds and summarization.
- You produce a custom figure-level synthesis rather than reproducing one locus-level exercise.
