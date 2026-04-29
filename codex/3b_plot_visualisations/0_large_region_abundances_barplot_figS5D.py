#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
# '''=================================================
# @Project: ACT_melanoma_neutrophil (data: act_neutrophil_data_repository)
# @Author : Gemma van der Voort
# @Desc   : fig 5 panel D: region abundances (stacked barplot for now)
# @Desc updated: CODEx region abundance stacked barplots.
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

PROJECT_ROOT = Path(__file__).resolve()
while PROJECT_ROOT.name != "ACT_melanoma_neutrophil" and PROJECT_ROOT.parent != PROJECT_ROOT:
    PROJECT_ROOT = PROJECT_ROOT.parent
if PROJECT_ROOT.name != "ACT_melanoma_neutrophil":
    raise RuntimeError("ACT_melanoma_neutrophil root not found")
sys.path.append(str(PROJECT_ROOT / "helper_files"))

from paths_parameters import fig_dir, modifier_base, parent_dir, processed_spatial_data_dir, cell_colours_hex

# load style
plt.style.use(str(PROJECT_ROOT / "helper_files" / "journal_style.mplstyle"))

#%% preamble
cell_type_col = 'Cell type'
cell_type_col_file_name = 'cell_type'  # for use in file names, no spaces or caps. also poss shorter
version_key = 'v4_8'  # version of the metacluster key, large scale classifier.
metacluster_key_orig = f'Metacluster {version_key}'
cn_col = 'CN_k50_n20'  # column name for the CN in adata
relative = False  # whether to plot relative or absolute abundances

majority_threshold = 0.5  # threshold for filtering metaclusters
metacluster_key = 'Metacluster v4_8_filtered_0.5_maxit_25_v4_8' # f"{metacluster_key_orig}_filtered_{majority_threshold}_{version_key}"
img_modifier = f'v4_corrected'  # for the cn clustering
fig_dir = f'{fig_dir}/manuscipt_plots/abundances'
os.makedirs(fig_dir, exist_ok=True)
adata = ad.read_h5ad(f"{processed_spatial_data_dir}/sa_{modifier_base}_inclregions_{img_modifier}_{metacluster_key}.h5ad")

# %% make a stacked barplot of the three regions. one bar per treatment, one colour per region.
# get the regions and their colours
region_order = adata.uns[f'{metacluster_key_orig}_order']
region_colours = adata.uns[f'{metacluster_key_orig}_colors']
conditions_codex = ['tumour_bearing', 'Cyclo_d1', 'ACT_d3', 'ACT_d7', 'ACT_d14']
x_axis_tick_labels_codex = ["Tumor\nbearing", "Cy", "ACT\nday 3", "ACT\nday 7", "ACT\nday 14"]
# Get region and treatment data
df = adata.obs[[metacluster_key, 'Treatment']].copy()

# Count cells per region and treatment
counts = df.groupby(['Treatment', metacluster_key]).size().unstack(fill_value=0)

# Reorder rows and columns
counts = counts.loc[:, region_order]
# Normalize to get proportions (optional)
counts_norm = counts.div(counts.sum(axis=1), axis=0)
# Plot
fig, ax = plt.subplots(figsize=(8, 6))
bottom = pd.Series([0] * len(counts), index=counts.index)

for region, color in zip(region_order, region_colours):
    ax.bar(counts.index, counts[region], bottom=bottom, label=region, color=color)
    bottom += counts[region]

ax.set_ylabel("normalised cell count")
ax.set_xlabel("treatment")

# Fix tick labels
ax.set_xticks(range(len(counts.index)))
ax.set_xticklabels(x_axis_tick_labels_codex)

sns.despine()

# Save
plt.tight_layout()
plt.savefig(f"{fig_dir}/stacked_barplot_regions_nolegend_{img_modifier}_{metacluster_key}.pdf", dpi=300)
plt.close()

#%% per organ per treatment stacked barplot
organs = ['inLNr', 'inLNl', 'brLNr']
organ_labels = ['tdLN', 'clLN', 'intLN']
# Get counts
df = adata.obs[[metacluster_key, 'Treatment', 'Organ']].copy()
counts = df.groupby(['Treatment', 'Organ', metacluster_key]).size().unstack(fill_value=0)
counts = counts.loc[(conditions_codex, organs), region_order]

# Normalize within bar
counts = counts.div(counts.sum(axis=1), axis=0)

# Bar positions
n_organs = len(organs)
x_positions = []
tick_labels = []
for i, treatment in enumerate(conditions_codex):
    for j, organ in enumerate(organ_labels):
        x_positions.append(i * (n_organs + 1) + j)
        tick_labels.append(organ)

# Plot
fig, ax = plt.subplots(figsize=(16, 6))
bottom = pd.Series([0] * len(counts), index=counts.index)

for region, color in zip(region_order, region_colours):
    ax.bar(x_positions, counts[region].values, bottom=bottom.values, label=region, color=color, width=0.8)
    bottom += counts[region].values
ax.set_ylabel("Cell count")
ax.set_xlabel("")
# Set organ tick labels
ax.set_xticks(x_positions)
ax.set_xticklabels(tick_labels, rotation=0, fontsize=14)

# Add treatment labels below groups
group_centers = [i * (n_organs + 1) + (n_organs - 1) / 2 for i in range(len(conditions_codex))]
y_min = -0.08  # adjust if needed
for center, label in zip(group_centers, x_axis_tick_labels_codex):
    ax.text(center, y_min, label, ha='center', va='top', fontsize=16, transform=ax.get_xaxis_transform())

ax.legend(title="region", bbox_to_anchor=(1.05, 1), loc='upper left')
sns.despine()

# Save
plt.tight_layout()
plt.savefig(f"{fig_dir}/stacked_barplot_by_organ_by_treatment_{img_modifier}_{metacluster_key}.pdf", dpi=300)
plt.close()

#%% incl errorbars
# Count per sample
df = adata.obs[[metacluster_key, 'Treatment', 'Organ', 'mouse_id']].copy()
grouped = df.groupby(['Treatment', 'Organ', 'mouse_id', metacluster_key]).size().unstack(fill_value=0)

# Reorder columns
grouped = grouped[region_order]

# drop rows that are zero for the region columns (since we currently copy mice to conditions they were not in)
grouped = grouped[(grouped[region_order].sum(axis=1) > 0)]
# save this df to file for use in blood project (for consistent visualisation style)
grouped.to_csv(f'{processed_spatial_data_dir}/region_counts_per_LN_{img_modifier}_{metacluster_key}.csv')
if relative:
    # Normalize per sample to get region proportions
    grouped_proportion = grouped.div(grouped.sum(axis=1), axis=0)
    grouped = grouped_proportion
# Compute mean and SEM across samples
means = grouped.groupby(['Treatment', 'Organ']).mean()
sems = grouped.groupby(['Treatment', 'Organ']).sem()

# Set up bar positions
n_organs = len(organs)
x_positions = []
tick_labels = []
for i, treatment in enumerate(conditions_codex):
    for j, organ in enumerate(organ_labels):
        x_positions.append(i * (n_organs + 1) + j)
        tick_labels.append(organ)

# Plot
fig, ax = plt.subplots(figsize=(16, 6))
bar_width = 0.8
bottom = np.zeros(len(x_positions))

for k, (region, color) in enumerate(zip(region_order, region_colours)):
    # Bar heights and error for current region
    region_means = []
    region_sems = []
    for treatment in conditions_codex:
        for organ in organs:
            region_means.append(means.loc[(treatment, organ), region] if (treatment, organ) in means.index else 0)
            region_sems.append(sems.loc[(treatment, organ), region] if (treatment, organ) in sems.index else 0)

    ax.bar(x_positions, region_means, bottom=bottom, yerr=region_sems, label=region,
           color=color, width=bar_width, capsize=3)
    bottom += region_means
if relative:
    ax.set_ylabel("Normalized cell count")
    # remove y tick labels greater than 1
    yticks = ax.get_yticks().tolist()
    yticks = [tick for tick in yticks if tick <= 1.0]
    ax.set_yticks(yticks)
else:
    ax.set_ylabel("Cell count")
ax.set_xticks(x_positions)
ax.set_xticklabels(tick_labels, rotation=0, fontsize=14)

# Treatment labels under organ bars
group_centers = [i * (n_organs + 1) + (n_organs - 1) / 2 for i in range(len(conditions_codex))]
y_min = -0.08
for center, label in zip(group_centers, x_axis_tick_labels_codex):
    ax.text(center, y_min, label, ha='center', va='top', fontsize=16, transform=ax.get_xaxis_transform())

ax.legend(title="region", bbox_to_anchor=(1.05, 1), loc='upper left')
sns.despine()

plt.tight_layout()
if relative:
    plt.savefig(f"{fig_dir}/stacked_barplot_by_organ_with_error_relative_{img_modifier}_{metacluster_key}.pdf", dpi=300)
else:
    plt.savefig(f"{fig_dir}/stacked_barplot_by_organ_with_error_absolute_{img_modifier}_{metacluster_key}.pdf", dpi=300)
plt.close()
