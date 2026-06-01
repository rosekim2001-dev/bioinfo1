#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SRC_COUNTS="$ROOT_DIR/../../binfo1-work/read-counts.txt"
OUT_DIR="$ROOT_DIR/03_results/counts"

mkdir -p "$OUT_DIR"

# Reuse previously generated featureCounts output from class work.
cp "$SRC_COUNTS" "$OUT_DIR/all_assays_counts.txt"

# Split one-column count files to match the expected skeleton outputs.
python3 - "$OUT_DIR/all_assays_counts.txt" "$OUT_DIR" <<'PY'
import sys
import pandas as pd

src, out_dir = sys.argv[1], sys.argv[2]
df = pd.read_csv(src, sep="\t", comment="#")

mapping = {
    "CLIP-35L33G.bam": "clip_counts.txt",
    "RPF-siLuc.bam": "rpf_siLuc_counts.txt",
    "RPF-siLin28a.bam": "rpf_siLin28a_counts.txt",
    "RNA-siLuc.bam": "rna_siLuc_counts.txt",
    "RNA-siLin28a.bam": "rna_siLin28a_counts.txt",
}

for col, out_name in mapping.items():
    out = df[["Geneid", col]].rename(columns={"Geneid": "gene_id", col: "count"})
    out.to_csv(f"{out_dir}/{out_name}", sep="\t", index=False)
PY

echo "Count outputs written to: $OUT_DIR"
