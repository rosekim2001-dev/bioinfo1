#!/usr/bin/env python3
import argparse
import pandas as pd


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--counts", required=True)
    p.add_argument("--annotation", required=True)
    p.add_argument("--out", required=True)
    args = p.parse_args()

    raw = pd.read_csv(args.counts, sep="\t", comment="#")
    ann = pd.read_csv(args.annotation, sep="\t")

    keep = [
        "Geneid",
        "CLIP-35L33G.bam",
        "RPF-siLuc.bam",
        "RPF-siLin28a.bam",
        "RNA-siLuc.bam",
        "RNA-siLin28a.bam",
    ]
    m = raw[keep].rename(
        columns={
            "Geneid": "gene_id",
            "CLIP-35L33G.bam": "clip_count",
            "RPF-siLuc.bam": "rpf_siLuc",
            "RPF-siLin28a.bam": "rpf_siLin28a",
            "RNA-siLuc.bam": "rna_siLuc",
            "RNA-siLin28a.bam": "rna_siLin28a",
        }
    )

    df = ann.merge(m, on="gene_id", how="inner")

    # Keep genes with evidence in either RNA or RPF conditions.
    expr = (df["rna_siLuc"] + df["rna_siLin28a"] + df["rpf_siLuc"] + df["rpf_siLin28a"]) >= 10
    df = df.loc[expr].copy()

    df.to_csv(args.out, sep="\t", index=False)
    print(f"Saved {len(df)} genes to {args.out}")


if __name__ == "__main__":
    main()
