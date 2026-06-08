#!/usr/bin/env python3
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


GROUP_ORDER = [
    "Q1_clip_high_te_up",
    "Q2_clip_low_te_up",
    "Q3_clip_low_te_down",
    "Q4_clip_high_te_down",
    "other",
]

GROUP_LABELS = {
    "Q1_clip_high_te_up": "CLIP high / TE up",
    "Q2_clip_low_te_up": "CLIP low / TE up",
    "Q3_clip_low_te_down": "CLIP low / TE down",
    "Q4_clip_high_te_down": "CLIP high / TE down",
    "other": "Other",
}

PALETTE = {
    "Q1_clip_high_te_up": "#2f7d32",
    "Q2_clip_low_te_up": "#4f6fb3",
    "Q3_clip_low_te_down": "#c24c3a",
    "Q4_clip_high_te_down": "#b88326",
    "other": "#9a9a9a",
}


def save_current(fig, outdir, stem):
    for ext in ("png", "pdf"):
        fig.savefig(outdir / f"{stem}.{ext}", dpi=300, bbox_inches="tight")
    plt.close(fig)


def clipped_metrics(df):
    out = df.copy()
    out["clip_plot"] = out["clip_enrichment_log2"].clip(-5, 5)
    out["te_plot"] = out["te_change_log2"].clip(-5, 5)
    return out


def plot_sample_counts(df, outdir):
    cols = ["clip_count", "rpf_siLuc", "rpf_siLin28a", "rna_siLuc", "rna_siLin28a"]
    labels = ["CLIP", "RPF siLuc", "RPF siLin28a", "RNA siLuc", "RNA siLin28a"]
    total = df[cols].sum()
    detected = (df[cols] > 0).sum()
    plot_df = pd.DataFrame({"sample": labels, "total_reads": total.values, "detected_genes": detected.values})

    fig, ax1 = plt.subplots(figsize=(8.8, 5.2))
    sns.barplot(data=plot_df, x="sample", y="total_reads", color="#567c8d", ax=ax1)
    ax1.set_ylabel("Total assigned reads")
    ax1.set_xlabel("")
    ax1.tick_params(axis="x", rotation=20)
    ax1.ticklabel_format(style="sci", axis="y", scilimits=(0, 0))

    ax2 = ax1.twinx()
    ax2.plot(plot_df["sample"], plot_df["detected_genes"], color="#c44e52", marker="o", linewidth=2)
    ax2.set_ylabel("Detected genes")
    ax2.set_ylim(0, max(plot_df["detected_genes"]) * 1.2)
    ax1.set_title("Read Assignment and Gene Detection by Assay")
    save_current(fig, outdir, "figure_01_sample_count_qc")


def plot_metric_distributions(df, outdir):
    p = clipped_metrics(df)
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.6))
    sns.histplot(p["clip_plot"], bins=70, kde=True, color="#567c8d", ax=axes[0])
    axes[0].axvline(0, color="black", linestyle="--", linewidth=1)
    axes[0].set_xlabel("CLIP enrichment log2, clipped to [-5, 5]")
    axes[0].set_ylabel("Genes")
    axes[0].set_title("CLIP Enrichment Distribution")

    sns.histplot(p["te_plot"], bins=70, kde=True, color="#b88326", ax=axes[1])
    axes[1].axvline(0, color="black", linestyle="--", linewidth=1)
    axes[1].set_xlabel("TE change log2, clipped to [-5, 5]")
    axes[1].set_ylabel("Genes")
    axes[1].set_title("Translation Efficiency Change Distribution")
    save_current(fig, outdir, "figure_02_metric_distributions")


def plot_gene_scatter(df, outdir):
    p = clipped_metrics(df)
    fig, ax = plt.subplots(figsize=(8.4, 7.2))
    for grp in GROUP_ORDER:
        sub = p.loc[p["group"] == grp]
        if sub.empty:
            continue
        alpha = 0.18 if grp == "other" else 0.6
        size = 9 if grp == "other" else 18
        ax.scatter(
            sub["clip_plot"],
            sub["te_plot"],
            s=size,
            alpha=alpha,
            c=PALETTE[grp],
            label=GROUP_LABELS[grp],
            linewidths=0,
        )
    ax.axhline(0, color="black", linestyle="--", linewidth=0.9)
    ax.axvline(0, color="black", linestyle="--", linewidth=0.9)
    ax.axhline(0.3, color="#777777", linestyle=":", linewidth=0.9)
    ax.axhline(-0.3, color="#777777", linestyle=":", linewidth=0.9)
    ax.axvline(0.5, color="#777777", linestyle=":", linewidth=0.9)
    ax.axvline(-0.5, color="#777777", linestyle=":", linewidth=0.9)
    ax.set_xlabel("CLIP enrichment log2, clipped to [-5, 5]")
    ax.set_ylabel("TE change siLin28a/siLuc log2, clipped to [-5, 5]")
    ax.set_title("Gene-Level CLIP Enrichment and TE Response")
    ax.legend(loc="upper right", frameon=False, fontsize=8)
    save_current(fig, outdir, "figure_03_gene_level_scatter")


def plot_group_counts(df, outdir):
    counts = df["group"].value_counts().reindex(GROUP_ORDER).fillna(0).astype(int).reset_index()
    counts.columns = ["group", "n_genes"]
    counts["label"] = counts["group"].map(GROUP_LABELS)

    fig, ax = plt.subplots(figsize=(8.4, 5))
    sns.barplot(
        data=counts,
        y="label",
        x="n_genes",
        hue="label",
        palette=dict(zip(counts["label"], [PALETTE[g] for g in counts["group"]])),
        legend=False,
        ax=ax,
    )
    ax.set_xlabel("Number of genes")
    ax.set_ylabel("")
    ax.set_title("Gene Counts in CLIP/TE Quadrants")
    for i, value in enumerate(counts["n_genes"]):
        ax.text(value + max(counts["n_genes"]) * 0.01, i, f"{value:,}", va="center", fontsize=9)
    save_current(fig, outdir, "figure_04_group_counts")


def plot_gene_type_composition(df, outdir):
    main_types = df["gene_type"].value_counts().head(8).index
    comp = df.assign(gene_type_plot=np.where(df["gene_type"].isin(main_types), df["gene_type"], "other_gene_type"))
    tab = pd.crosstab(comp["group"], comp["gene_type_plot"])
    tab = tab.reindex(GROUP_ORDER).fillna(0)
    tab_frac = tab.div(tab.sum(axis=1).replace(0, np.nan), axis=0).fillna(0)

    fig, ax = plt.subplots(figsize=(10.5, 5.8))
    bottom = np.zeros(len(tab_frac))
    colors = sns.color_palette("tab10", n_colors=tab_frac.shape[1])
    for idx, col in enumerate(tab_frac.columns):
        ax.barh([GROUP_LABELS[g] for g in tab_frac.index], tab_frac[col], left=bottom, color=colors[idx], label=col)
        bottom += tab_frac[col].to_numpy()
    ax.set_xlim(0, 1)
    ax.set_xlabel("Fraction of genes")
    ax.set_ylabel("")
    ax.set_title("Gene-Type Composition by Quadrant")
    ax.legend(loc="center left", bbox_to_anchor=(1.01, 0.5), frameon=False, fontsize=8)
    save_current(fig, outdir, "figure_05_gene_type_composition")


def plot_enrichment_heatmap(enrichment, outdir):
    enr = enrichment.copy()
    enr["score"] = -np.log10(enr["fdr"].clip(lower=1e-300))
    top_terms = enr.sort_values("score", ascending=False)["go_term"].drop_duplicates().head(14)
    heat = enr.loc[enr["go_term"].isin(top_terms)].pivot_table(
        index="go_term", columns="group", values="score", aggfunc="max", fill_value=0
    )
    heat = heat.reindex(top_terms)
    heat = heat[[g for g in GROUP_ORDER if g in heat.columns]]
    heat.columns = [GROUP_LABELS[g] for g in heat.columns]

    fig, ax = plt.subplots(figsize=(9.6, 6.8))
    sns.heatmap(heat, cmap="mako", linewidths=0.4, linecolor="white", cbar_kws={"label": "-log10(FDR)"}, ax=ax)
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_title("Top Category Enrichment Across Quadrants")
    save_current(fig, outdir, "figure_06_enrichment_heatmap")


def write_summary(df, enrichment, outdir):
    cols = ["clip_count", "rpf_siLuc", "rpf_siLin28a", "rna_siLuc", "rna_siLin28a"]
    sample_summary = pd.DataFrame({
        "sample": cols,
        "total_assigned_reads": df[cols].sum().astype(int).values,
        "detected_genes": (df[cols] > 0).sum().astype(int).values,
    })
    sample_summary.to_csv(outdir / "sample_qc_summary.tsv", sep="\t", index=False)

    group_summary = (
        df.groupby("group", dropna=False)
        .agg(
            n_genes=("gene_id", "size"),
            mean_clip_enrichment_log2=("clip_enrichment_log2", "mean"),
            median_clip_enrichment_log2=("clip_enrichment_log2", "median"),
            mean_te_change_log2=("te_change_log2", "mean"),
            median_te_change_log2=("te_change_log2", "median"),
        )
        .reindex(GROUP_ORDER)
        .reset_index()
    )
    group_summary.to_csv(outdir / "group_metric_summary.tsv", sep="\t", index=False)

    top = enrichment.sort_values(["fdr", "pvalue"]).head(12)
    top.to_csv(outdir / "top_enrichment_summary.tsv", sep="\t", index=False)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--metrics", required=True)
    p.add_argument("--enrichment", required=True)
    p.add_argument("--outdir", required=True)
    args = p.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    sns.set_theme(style="whitegrid", context="notebook")

    df = pd.read_csv(args.metrics, sep="\t")
    enrichment = pd.read_csv(args.enrichment, sep="\t")

    plot_sample_counts(df, outdir)
    plot_metric_distributions(df, outdir)
    plot_gene_scatter(df, outdir)
    plot_group_counts(df, outdir)
    plot_gene_type_composition(df, outdir)
    plot_enrichment_heatmap(enrichment, outdir)
    write_summary(df, enrichment, outdir)
    print(f"Saved extended figures and summaries to {outdir}")


if __name__ == "__main__":
    main()
