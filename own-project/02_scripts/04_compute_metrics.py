#!/usr/bin/env python3
import argparse
import numpy as np
import pandas as pd


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--matrix", required=True)
    p.add_argument("--out", required=True)
    args = p.parse_args()

    df = pd.read_csv(args.matrix, sep="\t")

    # RNA-normalized CLIP enrichment.
    clip_per_rna = (df["clip_count"] + 1.0) / (df[["rna_siLuc", "rna_siLin28a"]].mean(axis=1) + 1.0)
    df["clip_enrichment_log2"] = np.log2(clip_per_rna / np.median(clip_per_rna))

    te_siLuc = (df["rpf_siLuc"] + 1.0) / (df["rna_siLuc"] + 1.0)
    te_siLin28a = (df["rpf_siLin28a"] + 1.0) / (df["rna_siLin28a"] + 1.0)
    df["te_change_log2"] = np.log2(te_siLin28a / te_siLuc)

    x, y = df["clip_enrichment_log2"], df["te_change_log2"]
    df["group"] = "other"
    df.loc[(x > 0.5) & (y > 0.3), "group"] = "Q1_clip_high_te_up"
    df.loc[(x <= -0.5) & (y > 0.3), "group"] = "Q2_clip_low_te_up"
    df.loc[(x <= -0.5) & (y <= -0.3), "group"] = "Q3_clip_low_te_down"
    df.loc[(x > 0.5) & (y <= -0.3), "group"] = "Q4_clip_high_te_down"

    df.to_csv(args.out, sep="\t", index=False)
    print(f"Saved metrics for {len(df)} genes to {args.out}")


if __name__ == "__main__":
    main()
