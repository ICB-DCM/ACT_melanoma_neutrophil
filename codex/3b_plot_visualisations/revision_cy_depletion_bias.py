#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
# '''=================================================
# @Project -> File :ACT_melanoma_neutrophil -> revision_cy_depletion_bias.py
# @Author : Gemma van der Voort
# @Time   : 3/30/26 9:30 PM
# @Desc   : Investigation into Ki67 expression pre- and post-Cy treatment, as a proxy for proliferation and potential Cy
# sensitivity. Also includes visualisation of the distribution of Ki67 values, and the mean Ki67 per mouse,
# split by treatment and organ.
# '''=================================================
# %% imports
import sys
from pathlib import Path
import numpy as np
import anndata as ad
import matplotlib.pyplot as plt
import pandas as pd
import os
import seaborn as sns
from pathlib import Path
from typing import Optional, Tuple, Dict

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import sparse
from scipy.stats import mannwhitneyu

PROJECT_ROOT = Path(__file__).resolve()
while PROJECT_ROOT.name != "ACT_melanoma_neutrophil" and PROJECT_ROOT.parent != PROJECT_ROOT:
    PROJECT_ROOT = PROJECT_ROOT.parent
if PROJECT_ROOT.name != "ACT_melanoma_neutrophil":
    raise RuntimeError("ACT_melanoma_neutrophil root not found")
sys.path.append(str(PROJECT_ROOT / "helper_files"))

from paths_parameters import (fig_dir, modifier_base, parent_dir, processed_spatial_data_dir, x_axis_tick_labels_codex,
                              bar_colours_hex_codex, updated_ln_label_dict, marker_per_organ_dict)
from functions_vis import boxplot_unpaired, single_paired_plot, boxplot_paired


# load style
plt.style.use(str(PROJECT_ROOT / "helper_files" / "journal_style_old.mplstyle"))

#%% preamble
cell_type_col = 'Cell type'
cell_type_col_file_name = 'cell_type'  # for use in file names, no spaces or caps. also poss shorter
version_key = 'v4_8'  # version of the metacluster key, large scale classifier.
metacluster_key_orig = f'Metacluster {version_key}'
cn_col = 'CN_k50_n20'  # column name for the CN in adata
# relative = True  # whether to plot relative or absolute abundances
flavour = "tb_cyclo"
# flavour = 'all_conditions'
majority_threshold = 0.5  # threshold for filtering metaclusters
metacluster_col_name = f"{metacluster_key_orig}_filtered_{majority_threshold}_{version_key}"
img_modifier = f'v4_corrected'  # for the cn clustering
fig_dir = f'{fig_dir}/manuscipt_plots/revision'
os.makedirs(fig_dir, exist_ok=True)

adata = ad.read_h5ad(f"{processed_spatial_data_dir}/sa_{modifier_base}_inclregions_{img_modifier}_{version_key}.h5ad")

#%% plot ki67 histogram
ki67_values = adata[:, 'Ki67'].X
nonzero = ki67_values[ki67_values > 0]
thr = np.percentile(nonzero, 90) # set threshold at 90th percentile of nonzero values
fig, ax = plt.subplots()
sns.histplot(ki67_values, bins=100, ax=ax, legend=False)
ax.axvline(thr)
# print the percentage of cells above the threshold
percent_above_thr = (ki67_values > thr).mean() * 100
print(f"Percentage of cells above threshold: {percent_above_thr:.2f}%")
# %% add percentage to plot as text
ax.text(thr + 0.1 , ax.get_ylim()[1] * 0.9, f"{percent_above_thr:.2f}%", ha='left', va='top')
ax.set_yscale('log')
sns.despine()
fig.tight_layout()
fig.savefig(f"{fig_dir}/ki67_histogram_for_revision.pdf", dpi=300, bbox_inches="tight")

#%%
df_images = adata.obs[["dataset_name", "Treatment", "Organ"]].drop_duplicates()
counts = df_images.value_counts(["Treatment", "Organ"])

#%% continious & thresholded ki67 to obs for easiery handling:
ki67 = adata[:, "Ki67"].X
if sparse.issparse(ki67):
    ki67 = ki67.toarray()
ki67 = np.asarray(ki67).flatten()

adata.obs["ki67_pos"] = ki67 > thr
adata.obs["ki67"] = ki67

#%% does Cy cause biased cell death, meaning A: less cells express Ki67

if flavour == "tb_cyclo":
    df = adata.obs[
        adata.obs["Treatment"].isin(["tumour_bearing", "Cyclo_d1"])
    ].copy()
    x_axis_tick_labels_codex = x_axis_tick_labels_codex[:2]
else:
    df = adata.obs.copy()
df["ki67_log"] = np.log1p(df["ki67"])

per_mouse = (
    df.groupby(["Treatment", "mouse_id"])
    .agg(
        ki67_pos_frac=("ki67_pos", "mean"),
        ki67_mean=("ki67_log", "mean"),
        ki67_median=("ki67_log", "median"),  # optional, often more robust
    )
    .reset_index()
)
for test in ["ki67_pos_frac", "ki67_mean", "ki67_median"]:
    tb = per_mouse.loc[per_mouse["Treatment"] == "tumour_bearing", test]
    cyclo = per_mouse.loc[per_mouse["Treatment"] == "Cyclo_d1", test]
    tb = tb.dropna()  # groupby created nans where treatment/mice combo's didn't exist, drop them.
    cyclo = cyclo.dropna()
    print(f"Overall: {len(tb)} tumour_bearing mice, {len(cyclo)} cyclo mice")
    result_over_treatment = mannwhitneyu(tb, cyclo)
    print(f"{test}: {result_over_treatment}")

#%% 2: Does Cy effect differ between LN types? “Is bias uniform or LN-specific?”
per_mouse_ln = (
    df.groupby(["Treatment", "mouse_id", "Organ"])
    .agg(
        ki67_pos_frac=("ki67_pos", "mean"),
        ki67_mean=("ki67_log", "mean"),
        ki67_median=("ki67_log", "median"),
    )
    .reset_index()
)
for test in ["ki67_pos_frac", "ki67_mean", "ki67_median"]:
    result_per_organ = {}
    for organ in per_mouse_ln["Organ"].unique():
        sub = per_mouse_ln[per_mouse_ln["Organ"] == organ]
        tb = sub.loc[sub["Treatment"] == "tumour_bearing", test]
        cyclo = sub.loc[sub["Treatment"] == "Cyclo_d1", test]
        tb = tb.dropna()
        cyclo = cyclo.dropna()
        print(f"{organ}: {len(tb)} tumour_bearing, {len(cyclo)} cyclo")

        if len(tb) > 0 and len(cyclo) > 0:
            result_per_organ[organ] = mannwhitneyu(tb, cyclo)
        else:
            result_per_organ[organ] = "Not enough data"
        print(f"{organ} {test}: {result_per_organ[organ]}")

#%% so not significant, but trends are there. sorting shows clear separation of treatment groups.
per_mouse_sorted = per_mouse.sort_values("ki67_mean").dropna(subset=["ki67_mean"], inplace=False)
per_mouse_ln_sorted = per_mouse_ln.sort_values("ki67_mean").dropna(subset=["ki67_mean"], inplace=False)
#%% visual checks:
if flavour == "tb_cyclo":
    condition_order = ["tumour_bearing", "Cyclo_d1"]
else:
    condition_order = ["tumour_bearing", "Cyclo_d1", "ACT_d3", "ACT_d7", "ACT_d14"]
condition_palette = dict(zip(condition_order, bar_colours_hex_codex[:len(condition_order)]))
condition_label_dict = dict(zip(condition_order, x_axis_tick_labels_codex[:len(condition_order)]))
ln_order = ["inLNr", "inLNl", "brLNr"]
ln_label_dict = {organ: updated_ln_label_dict.get(organ, organ) for organ in ln_order}

#%% plot the data without stats (since non-significant) in manuscript plotting style.
#A: without LN-split
group_col = 'Treatment'
value_col = 'ki67_mean'
value_col_label = 'Mean log(Ki67) per mouse'
# drop unused categories from per_mouse_sorted.Treatment
per_mouse_sorted.Treatment = per_mouse_sorted.Treatment.cat.remove_unused_categories()
# set category order & order the df by condition order (since palette gets confused)
per_mouse_sorted.Treatment = per_mouse_sorted.Treatment.cat.set_categories(condition_order, ordered=True)
per_mouse_ordered = per_mouse_sorted.sort_values("Treatment")
n_categories = len(per_mouse_sorted[group_col].unique())
fig, ax = plt.subplots(1, 1, figsize=(n_categories + 2, 6))
boxplot_unpaired(per_mouse_ordered, group_col, value_col, stats_df=None, ax=ax, p_val_col_name=None,
                 strip_kwargs={'palette': condition_palette, 'marker': '8'})
ax.set_title('')  # remove title
ax.set_ylabel(value_col_label)
ax.set_xlabel("")
ax.set_xticklabels(x_axis_tick_labels_codex)
ax.set_xlabel('')  # remove xlabel
plt.tight_layout()
plt.savefig(f"{fig_dir}/ki67_mean_per_mouse_{flavour}_v2.pdf", dpi=300, bbox_inches="tight")
plt.close()
# %% B: with LN-split
ncols = 1
organ_col = 'Organ'
markers_LN = [marker_per_organ_dict[ln] for ln in ln_order]
# drop unused categories
per_mouse_ln_sorted.Treatment = per_mouse_ln_sorted.Treatment.cat.remove_unused_categories()
per_mouse_ln_sorted.Treatment = per_mouse_ln_sorted.Treatment.cat.set_categories(condition_order, ordered=True)
per_mouse_ln_ordered = per_mouse_ln_sorted.sort_values("Treatment")
fig, ax = plt.subplots(1, ncols , figsize=(((n_categories*1.5)+1), 6))
color_list =list(condition_palette.values()) # dict doesn't work for this function (would need to generalise, low priority).
boxplot_paired(per_mouse_ln_ordered, 'Treatment', value_col, x_group_col='Organ', stats_df=None, ax=ax,
               x_group_order=ln_order, markers=markers_LN, condition_colours=color_list,
               omit_ns=True
               )
ax.set_ylabel(value_col_label)
ax.set_xlabel('')
ax.set_xticklabels(x_axis_tick_labels_codex)
plt.tight_layout()
plt.savefig(f"{fig_dir}/ki67_mean_per_organ_{flavour}.pdf", dpi=300, bbox_inches="tight")
plt.close()
#%% add distribution (KDE) plot
def plot_distribution(
    df,
    figsize=(6, 5),
    xlabel="log(Ki67)",
    ylabel="Density",
    kde_kwargs=None,
):
    if kde_kwargs is None:
        kde_kwargs = {}

    fig, ax = plt.subplots(figsize=figsize)

    sns.kdeplot(
        data=df,
        x="ki67_log",
        hue="Treatment",
        hue_order=condition_order,
        palette=condition_palette,
        common_norm=False,
        ax=ax,
        **kde_kwargs
    )

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    legend = ax.get_legend()
    if legend is not None:
        legend.set_title("Treatment")
        for text, treatment in zip(legend.texts, condition_order):
            text.set_text(condition_label_dict[treatment])
        # increase font size of legend
        for text in legend.texts:
            text.set_fontsize(18)
        # increase font size of legend title
        legend.get_title().set_fontsize(20)

    sns.despine()
    fig.tight_layout()
    return fig

fig_dist = plot_distribution(df, figsize=(9, 6))
fig_dist.savefig(f"{fig_dir}/ki67_distribution_{flavour}.pdf", dpi=300, bbox_inches="tight")
plt.close(fig_dist)
