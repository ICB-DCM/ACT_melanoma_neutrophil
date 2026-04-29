#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
# '''=================================================
# @Project: ACT_melanoma_neutrophil (data: act_neutrophil_data_repository)
# @Author : Gemma van der Voort
# @Desc   : vis of the flow data, unpaired, no CpG. generates a plot per cell type consisting of 3 subplots,
# @Desc updated: Unpaired pan-immune flow visualizations for main conditions.
# one for each LN type. title is the LN type, remove if needed (also visible in marker shape)
# update 13.03: mtc revisited, as in blood_new_vis on feb 20th.
# update 8.8.2025: removed titles and x-labels, revisited mtc and test, fixed pmel abundances
# also updated to us english
# '''=================================================
import os
import sys
from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt
import math

PROJECT_ROOT = Path(__file__).resolve()
while PROJECT_ROOT.name != "ACT_melanoma_neutrophil" and PROJECT_ROOT.parent != PROJECT_ROOT:
    PROJECT_ROOT = PROJECT_ROOT.parent
if PROJECT_ROOT.name != "ACT_melanoma_neutrophil":
    raise RuntimeError("ACT_melanoma_neutrophil root not found")
sys.path.append(str(PROJECT_ROOT / "helper_files"))

from paths_parameters import x_axis_tick_labels_flow, bar_colours_hex_blood_flow, lymph_nodes_label, marker_per_organ_dict, data_repo_path
from functions_stat_test import unpaired_multiple_conditions  # unpaired_test
from functions_vis import boxplot_unpaired

plt.style.use(str(PROJECT_ROOT / "helper_files" / "journal_style.mplstyle"))

# paths etc
processed_data_path = f'{data_repo_path}/processed/flow_cytometry'
stats_dir = f'{processed_data_path}/statistics'
os.makedirs(stats_dir, exist_ok=True)
fig_path = f'{data_repo_path}/figures/flow_cytometry/main_conditions'
ln_df = pd.read_csv(f'{processed_data_path}/ln_pan_immune_normalised_v2_cd45_corrected.csv')
test_function_used = 'unpaired_multiple_conditions'
modifier = 'cd45_corrected_v2'
mtc = True
non_parametric = False
force_parametric = True # set after initial normalicy testing
print(f'Using {test_function_used} with multiple testing correction: {mtc}')
# make dir and subdir single_paired if it doesn't exist
os.makedirs(fig_path, exist_ok=True)
os.makedirs(f'{fig_path}/unpaired', exist_ok=True)
#%%
# drop condition if it contains the substring 'no CpG'
ln_df = ln_df[~ln_df['Condition'].str.contains('no CpG')]
df_condition_names = ln_df['Condition'].unique()
df_condition_names_pmel = df_condition_names[-4:]  # last 4 conditions are the pmel conditions, so only compare those for pmel T cells
value_cols = ['DCs percent', 'CD11b+ Dcs percent', 'Xcr1+ Dcs percent',
              'CD8 T percent', 'Endo CD8 T percent', 'Pmel T percent',
              'CD4 Tconv percent', 'CD25+ CD4+ T percent',
              'Neutrophils percent', 'Nk cells percent', 'B cells percent']
value_cols_y_label = ['DCs\n% of CD45.2$^+$ cells', 'CD11b$^+$ DCs\n% of CD45.2$^+$ cells', 'Xcr1$^+$ DCs\n% of CD45.2$^+$ cells',
                        'CD8 T\n% of CD45.2$^+$ cells', 'Endogenous CD8 T\n% of CD45.2$^+$ cells', 'Pmel-1 CD8 T\n% of CD45.2$^+$ cells',
                        'Conventional CD4 T\n% of CD45.2$^+$ cells', 'CD25$^+$ CD4 T\n% of CD45.2$^+$ cells',
                        'Neutrophils\n% of CD45.2$^+$ cells', 'NK cells\n% of CD45.2$^+$ cells', 'B cells\n% of CD45.2$^+$ cells'
                      ]  # plural except when ending in T (is short for T cells, but that takes too much space)

group_col = 'Condition'
conditions_to_compare_general = [(df_condition_names[i], df_condition_names[i + 1]) for i in
                         range(len(df_condition_names) - 1)]
conditions_to_compare_pmel = conditions_to_compare_general[-3:]  # only compare the last 3 conditions for Pmel T cells
ln_order = ['inLNr', 'inLNl', 'brLNr']
# df with the results of the statistical tests.
stats_df_columns = ['Cell type', 'LN_type', 'Condition 1', 'Condition 2', 'p-value', 'test_used', 'mtc_res']
stats_df = pd.DataFrame(columns=stats_df_columns)
for value_col in value_cols:
    for i_LN, LN in enumerate(ln_order):
        # for each condition, compare it to the next one
        if value_col == 'Pmel T percent':
            conditions_to_compare = conditions_to_compare_pmel
            df_1ln = ln_df[ln_df['LN_type'] == LN]
            df_1ln = df_1ln[df_1ln['Condition'].isin(df_condition_names_pmel)]  # only take the last 4 conditions for Pmel T cells
        else:
            conditions_to_compare = conditions_to_compare_general
            # continue
            df_1ln = ln_df[ln_df['LN_type'] == LN]

        for condition1, condition2 in conditions_to_compare:
            t_stat, p_val, test_used, mtc_res = unpaired_multiple_conditions(df_1ln, value_col, group_col,
                                                     (condition1, condition2), non_parametric=non_parametric,
                                                                    multiple_testing_correction=mtc,
                                                                             force_parametric=force_parametric
                                                                             )
            stats_df = pd.concat([stats_df, pd.DataFrame([[value_col, LN, condition1, condition2, p_val,
                                                           test_used, mtc_res]],
                                                         columns=stats_df_columns)])
# write to file
stats_df.to_csv(f'{stats_dir}/ln_pan_immune_stats_no_CpG_{test_function_used}_mtc{mtc}_np_{non_parametric}_fp_{force_parametric}_fdr_by_{modifier}.csv')
#%% plot the data
n_categories = len(ln_df['Condition'].unique())
ncols = len(ln_order)
fig_width = (n_categories-1)*ncols
markers_LN = [marker_per_organ_dict[ln] for ln in ln_order]

for cell_type_i, cell_type in enumerate(value_cols):
    if cell_type == 'Pmel T percent':
        continue  # skip pmel t cells and do outside of loop, since they need different fig dimensions and colours
    fig, axes = plt.subplots(1, ncols, figsize=(fig_width, 6), sharey=True)
    axes = axes.flatten()
    # get the highest value for y-axis limit
    max_y = math.ceil(ln_df[cell_type].max())
    # for each lymph node in LN_type, plot the cell type count
    for i_LN, LN in enumerate(ln_order):
        ax = axes[i_LN]
        data = ln_df[(ln_df['LN_type'] == LN)]  # [cell_type]
        # get stats df for cell type and LN
        stats_df_LN = stats_df[(stats_df['Cell type'] == cell_type) & (stats_df['LN_type'] == LN)]
        stats_df_LN = stats_df_LN.drop(['Cell type', 'LN_type'], axis=1)  # drop the LN_type column for the plot
        if mtc:
            p_val_col_name = 'mtc_res'
        else:
            p_val_col_name = 'p-value'
        boxplot_unpaired(data, group_col, cell_type, stats_df=stats_df_LN, ax=ax, p_val_col_name=p_val_col_name,
                         strip_kwargs = {'palette':bar_colours_hex_blood_flow, 'marker':markers_LN[i_LN]}, y_lim=max_y)
        # ax.set_title(f'{lymph_nodes_label[i_LN]}', y=1.05)
        ax.set_title('')  # remove title
        ax.set_ylabel(f'{value_cols_y_label[cell_type_i]}')
        # ax.set_xlabel(group_col)
        ax.set_xlabel('')  # remove xlabel
        ax.set_xticklabels(x_axis_tick_labels_flow)
    plt.tight_layout()
    plt.savefig(f'{fig_path}/unpaired/{cell_type}_flow_boxplot_no_CpG_{test_used}_mtc{mtc}_np_{non_parametric}_fp_{force_parametric}_{modifier}.pdf')
    plt.close()

# %%now for pmels
pmel_cell_type = 'Pmel T percent'
pmel_conditions = ln_df['Condition'].unique()[-4:]  # get the final 4 conditions
pmel_palette = bar_colours_hex_blood_flow[-4:]  # remove the first 3 colours, since they are not relevant for pmel
pmel_x_labels = x_axis_tick_labels_flow[-4:]  # remove the first 3 labels, since they are not relevant for pmel

fig_width_pmel = len(pmel_conditions) * len(ln_order)
fig, axes = plt.subplots(1, ncols, figsize=(fig_width_pmel, 6), sharey=True)
axes = axes.flatten()
max_y = math.ceil(ln_df[pmel_cell_type].max())

for i_LN, LN in enumerate(ln_order):
    ax = axes[i_LN]
    data = ln_df[(ln_df['LN_type'] == LN) & (ln_df['Condition'].isin(pmel_conditions))]
    stats_df_LN = stats_df[(stats_df['Cell type'] == pmel_cell_type) &
                           (stats_df['LN_type'] == LN) &
                           (stats_df['Condition 1'].isin(pmel_conditions)) &
                           (stats_df['Condition 2'].isin(pmel_conditions))]
    stats_df_LN = stats_df_LN.drop(['Cell type', 'LN_type'], axis=1)

    if mtc:
        p_val_col_name = 'mtc_res'
    else:
        p_val_col_name = 'p-value'

    boxplot_unpaired(data, group_col, pmel_cell_type, stats_df=stats_df_LN, ax=ax,
                     p_val_col_name=p_val_col_name,
                     strip_kwargs={'palette': pmel_palette, 'marker': markers_LN[i_LN]},
                     y_lim=max_y)
    ax.set_ylabel(f'{value_cols_y_label[value_cols.index(pmel_cell_type)]}')
    ax.set_xticklabels(pmel_x_labels)
    ax.set_title('')  # remove title
    ax.set_xlabel('')  # remove xlabel

plt.tight_layout()
plt.savefig(
    f'{fig_path}/unpaired/{pmel_cell_type}_flow_boxplot_final4_no_CpG_{test_used}_mtc{mtc}_np_{non_parametric}_fp_{force_parametric}_{modifier}.pdf')
plt.close()
