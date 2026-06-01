#!/usr/bin/env python3
import argparse
import math
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests


def fdr_table(metrics_df, universe_col="gene_id"):
    terms = sorted(metrics_df["gene_type"].dropna().unique())
    total = len(metrics_df)

    rows = []
    for grp, gdf in metrics_df.groupby("group"):
        if grp == "other":
            continue
        gset = set(gdf[universe_col])
        gsize = len(gset)
        for t in terms:
            term_set = set(metrics_df.loc[metrics_df["gene_type"] == t, universe_col])
            a = len(gset & term_set)
            b = gsize - a
            c = len(term_set) - a
            d = total - (a + b + c)
            _, p = fisher_exact([[a, b], [c, d]], alternative="greater")
            rows.append({
                "group": grp,
                "go_term": f"gene_type:{t}",
                "n_genes_in_term": len(term_set),
                "overlap": a,
                "group_size": gsize,
                "pvalue": p,
            })

    out = pd.DataFrame(rows)
    out["fdr"] = multipletests(out["pvalue"], method="fdr_bh")[1]
    return out.sort_values(["fdr", "pvalue"]) 


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--metrics", required=True)
    p.add_argument("--enrich_out", required=True)
    p.add_argument("--bubble_out", required=True)
    args = p.parse_args()

    df = pd.read_csv(args.metrics, sep="\t")
    enr = fdr_table(df)
    enr.to_csv(args.enrich_out, sep="\t", index=False)

    sig = enr.loc[(enr["fdr"] < 0.2) & (enr["overlap"] >= 20)].copy()
    if sig.empty:
        sig = enr.head(20).copy()

    bubble_rows = []
    for _, r in sig.iterrows():
        gt = r["go_term"].split(":", 1)[1]
        term_genes = df.loc[df["gene_type"] == gt]
        bubble_rows.append({
            "group": r["group"],
            "go_term": r["go_term"],
            "n_genes_in_term": int(r["n_genes_in_term"]),
            "x_go": float(term_genes["clip_enrichment_log2"].mean()),
            "y_go": float(term_genes["te_change_log2"].mean()),
            "fdr": float(r["fdr"]),
            "neglog10_fdr": -math.log10(max(float(r["fdr"]), 1e-300)),
            "overlap": int(r["overlap"]),
        })

    pd.DataFrame(bubble_rows).to_csv(args.bubble_out, sep="\t", index=False)
    print(f"Saved enrichment: {args.enrich_out}")
    print(f"Saved bubble input: {args.bubble_out}")


if __name__ == "__main__":
    main()
