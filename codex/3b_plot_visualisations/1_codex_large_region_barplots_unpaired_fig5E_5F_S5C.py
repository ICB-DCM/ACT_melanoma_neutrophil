#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
# '''=================================================
# @Project: ACT_melanoma_neutrophil (data: act_neutrophil_data_repository)
# @Author : Gemma van der Voort
# @Desc   : revisit fig 5E with less crowded plot in the style of the other plots in the manuscript.
# @Desc updated: Unpaired large-region barplots for CODEx.
# unpaired barplots of the 3 LNs over all conditions, one per large region (T, B, Med). Plot both absolute and
# relative abundances. main: T zone, inLNr. supplement: B, Med, also inLNr. relative, but check absolute as well.
# update jan 2026: larger font sizes for plotting, here manually adjusted to be a bit smaller to fit.
# '''=================================================
import sys
from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt
from pandas.api.types import CategoricalDtype
import math
import os

PROJECT_ROOT = Path(__file__).resolve()
while PROJECT_ROOT.name != "ACT_melanoma_neutrophil" and PROJECT_ROOT.parent != PROJECT_ROOT:
    PROJECT_ROOT = PROJECT_ROOT.parent
if PROJECT_ROOT.name != "ACT_melanoma_neutrophil":
    raise RuntimeError("ACT_melanoma_neutrophil root not found")
sys.path.append(str(PROJECT_ROOT / "helper_files"))

from paths_parameters import x_axis_tick_labels_codex, bar_colours_hex_codex, organ_names_short, \
    marker_per_organ_dict, data_repo_path
from functions_stat_test import unpaired_multiple_conditions  # unpaired_test
from functions_vis import boxplot_unpaired

plt.style.use(str(PROJECT_ROOT / "helper_files" / "journal_style_old.mplstyle"))

# paths etc
processed_spatial_data_dir = f'{data_repo_path}/processed/codex'
stats_dir = f'{processed_spatial_data_dir}/statistics'
os.makedirs(stats_dir, exist_ok=True)
modifier_base = f'20241021_all_ln_tree6_rep2'
file_modifier = 'v4_corrected'
version_key = 'v4_8'
metacluster_key = 'Metacluster v4_8_filtered_0.5_maxit_25_v4_8'
df_grouped_regions = pd.read_csv(f'{processed_spatial_data_dir}/region_counts_per_LN_{file_modifier}_'
                                 f'{metacluster_key}.csv')

fig_dir = f'{data_repo_path}/figures/codex/manuscipt_plots/large_regions'
os.makedirs(fig_dir, exist_ok=True)

test_function_used = 'unpaired_multiple_conditions'
mtc = True
non_parametric = False  # decided on parametric after testing
force_parametric = True # set after initial normalicy testing
print(f'Using {test_function_used} with multiple testing correction: {mtc}')
modifier_fig = f'{file_modifier}_{metacluster_key}_force_parametric_{force_parametric}_mtc_{mtc}_nonparam_{non_parametric}'
regions= ['B-cell follicle', 'T-cell zone', 'Medulla-Interfollicular zone-SCS']
conditions_ordered = ['tumour_bearing', 'Cyclo_d1', 'ACT_d3', 'ACT_d7', 'ACT_d14']
# rename treatment to condition for consistency
df_grouped_regions = df_grouped_regions.rename(columns={'Treatment': 'Condition'})
# reorder according to conditions_ordered
df_grouped_regions['Condition'] = pd.Categorical(df_grouped_regions['Condition'], categories=conditions_ordered, ordered=True)
df_grouped_regions = df_grouped_regions.sort_values('Condition') # actually reorder the df
#%% add the proportion columns to the df_grouped_all_cells
df_grouped_regions['Total cells in organ'] = df_grouped_regions[regions].sum(axis=1)
# for each region, make a new column with the proportion of cells in that region over total cells in organ
for region in regions:
    df_grouped_regions[f'{region}_proportion'] = df_grouped_regions[region] / df_grouped_regions['Total cells in organ']
# %%version 1: absolute counts
# value_cols = regions
# value_cols_y_labels = [f'Total cells\nin {region}' for region in regions]
# modifier_fig = f'{modifier_fig}_absolute_counts'
# version 2: relative counts
value_cols = [f'{region}_proportion' for region in regions]
value_cols_y_labels = [f'Proportion of cells\nin {region}' for region in regions]
modifier_fig = f'{modifier_fig}_relative_counts'
# general
ln_order = ['inLNr', 'inLNl', 'brLNr']
organ_col = 'Organ'
condition_col = 'Condition'
conditions_to_compare = [(conditions_ordered[i], conditions_ordered[i + 1]) for i in
                         range(len(conditions_ordered) - 1)]
# %% df with the results of the statistical tests.
stats_df_columns = ['Cell type', organ_col, 'Condition 1', 'Condition 2', 'p-value', 'test_used', 'mtc_res']
stats_df = pd.DataFrame(columns=stats_df_columns)
for value_col in value_cols:
    for i_LN, LN in enumerate(ln_order):
        df_1ln = df_grouped_regions[df_grouped_regions[organ_col] == LN]
        for condition1, condition2 in conditions_to_compare:
            t_stat, p_val, test_used, mtc_res = unpaired_multiple_conditions(df_1ln, value_col, condition_col,
                                                                             (condition1, condition2),
                                                                             non_parametric=non_parametric,
                                                                             multiple_testing_correction=mtc,
                                                                             force_parametric=force_parametric
                                                                             )
            stats_df = pd.concat([stats_df, pd.DataFrame([[value_col, LN, condition1, condition2, p_val,
                                                           test_used, mtc_res]],
                                                         columns=stats_df_columns)])
#%% write to file
stats_df.to_csv(f'{stats_dir}/large_regions_unpaired_stats_{modifier_fig}.csv', index=False)
#%% plot the data
n_categories = len(df_grouped_regions['Condition'].unique())
ncols = len(ln_order)
fig_width = (n_categories-1)*ncols
markers_LN = [marker_per_organ_dict[ln] for ln in ln_order]

for cell_type_i, cell_type in enumerate(value_cols):
    fig, axes = plt.subplots(1, ncols, figsize=(fig_width, 6), sharey=True)
    axes = axes.flatten()
    # get the highest value for y-axis limit
    max_y = math.ceil(df_grouped_regions[cell_type].max())
    # for each lymph node in LN_type, plot the cell type count
    for i_LN, LN in enumerate(ln_order):
        ax = axes[i_LN]
        data = df_grouped_regions[(df_grouped_regions[organ_col] == LN)]  # [cell_type]
        # get stats df for cell type and LN
        stats_df_LN = stats_df[(stats_df['Cell type'] == cell_type) & (stats_df[organ_col] == LN)]
        stats_df_LN = stats_df_LN.drop(['Cell type', organ_col], axis=1)  # drop the LN_type column for the plot
        if mtc:
            p_val_col_name = 'mtc_res'
        else:
            p_val_col_name = 'p-value'
        boxplot_unpaired(data, condition_col, cell_type, stats_df=None, ax=ax, p_val_col_name=p_val_col_name,
                         strip_kwargs = {'palette':bar_colours_hex_codex, 'marker':markers_LN[i_LN]}, y_lim=max_y)
        # ax.set_title(f'{lymph_nodes_label[i_LN]}', y=1.05)
        ax.set_title('')  # remove title
        ax.set_ylabel(f'{value_cols_y_labels[cell_type_i]}', size=21)
        # ax.set_xlabel(organ_col)
        ax.set_xlabel('')  # remove xlabel
        ax.set_xticklabels(x_axis_tick_labels_codex, size=15)
    plt.tight_layout()
    plt.savefig(f'{fig_dir}/{cell_type}_flow_boxplot_large_regions_{modifier_fig}_no_sign.pdf')
    plt.close()
