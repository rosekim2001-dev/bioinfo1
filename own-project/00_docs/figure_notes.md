# Figure Notes

## Figure Target
Bubble scatter summarizing GO-like term behavior in CLIP enrichment vs translation change space.

## Axes
- X-axis: mean `clip_enrichment_log2` across genes in the term.
- Y-axis: mean `te_change_log2` across genes in the term.

## Visual Encodings
- Bubble size: `n_genes_in_term`.
- Bubble color: `-log10(FDR)`.
- Reference lines: `x=0` and `y=0` dashed lines.

## Labels
- Label top 12 terms by highest `-log10(FDR)`.
- Use short term text (`go_term`) and small offset annotations.

## Current Implementation Notes
- Plot script: `02_scripts/06_plot_bubble.py`.
- Exported files:
  - `03_results/figures/figure_go_clip_ribo_bubble.png`
  - `03_results/figures/figure_go_clip_ribo_bubble.pdf`

## Interpretation Guide
- Right side (`x > 0`): relatively stronger CLIP enrichment.
- Upper side (`y > 0`): relatively increased TE in `siLin28a` vs `siLuc`.
- Large, darker bubbles: larger terms with stronger enrichment significance.

## Limitations
- Current enrichment is a local proxy based on `gene_type` categories.
- For publication-grade biology interpretation, replace with true GO database enrichment and preserve the same bubble schema.
