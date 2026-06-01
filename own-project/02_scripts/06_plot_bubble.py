#!/usr/bin/env python3
import argparse
import pandas as pd
import matplotlib.pyplot as plt


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--bubble", required=True)
    p.add_argument("--png", required=True)
    p.add_argument("--pdf", required=True)
    args = p.parse_args()

    df = pd.read_csv(args.bubble, sep="\t")

    plt.figure(figsize=(11, 8))
    sizes = 40 + (df["n_genes_in_term"] / max(df["n_genes_in_term"].max(), 1)) * 900
    sc = plt.scatter(
        df["x_go"],
        df["y_go"],
        s=sizes,
        c=df["neglog10_fdr"],
        cmap="viridis",
        alpha=0.8,
        edgecolor="black",
        linewidth=0.4,
    )

    top = df.sort_values("neglog10_fdr", ascending=False).head(12)
    for _, r in top.iterrows():
        plt.annotate(
            r["go_term"],
            (r["x_go"], r["y_go"]),
            xytext=(5, 5),
            textcoords="offset points",
            fontsize=8,
        )

    plt.axhline(0, color="gray", linestyle="--", linewidth=1)
    plt.axvline(0, color="gray", linestyle="--", linewidth=1)
    plt.xlabel("Mean CLIP enrichment (log2)")
    plt.ylabel("Mean TE change siLin28a/siLuc (log2)")
    plt.title("GO-like Bubble Summary: CLIP Enrichment vs Translation Change")
    cbar = plt.colorbar(sc)
    cbar.set_label("-log10(FDR)")
    plt.tight_layout()
    plt.savefig(args.png, dpi=300)
    plt.savefig(args.pdf)
    print(f"Saved: {args.png}")
    print(f"Saved: {args.pdf}")


if __name__ == "__main__":
    main()
