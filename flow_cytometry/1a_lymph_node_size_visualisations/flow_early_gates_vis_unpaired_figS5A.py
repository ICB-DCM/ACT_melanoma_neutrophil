#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
# '''=================================================
# @Project: ACT_melanoma_neutrophil (data: act_neutrophil_data_repository)
# @Author : Gemma van der Voort
# @Desc   : use the preprocessed, summed data for fig 5D. do unpaired tests.
# @Desc updated: Unpaired early-gate lymph node size visualizations.
# update: us english, remove x-axis title and title
# update sep 2025: add noCPG conditions
# '''=================================================
import math
import os
import sys
from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt

import numpy as np

PROJECT_ROOT = Path(__file__).resolve()
while PROJECT_ROOT.name != "ACT_melanoma_neutrophil" and PROJECT_ROOT.parent != PROJECT_ROOT:
    PROJECT_ROOT = PROJECT_ROOT.parent
if PROJECT_ROOT.name != "ACT_melanoma_neutrophil":
    raise RuntimeError("ACT_melanoma_neutrophil root not found")
sys.path.append(str(PROJECT_ROOT / "helper_files"))

from paths_parameters import x_axis_tick_labels_flow, bar_colours_hex_blood_flow, lymph_nodes_label, marker_per_organ_dict, data_repo_path
from functions_stat_test import ratio_paired_test, unpaired_multiple_conditions
from functions_vis import boxplot_paired, single_paired_plot, boxplot_unpaired

plt.style.use(str(PROJECT_ROOT / "helper_files" / "journal_style.mplstyle"))

processed_data_path = f'{data_repo_path}/processed/flow_cytometry'
stats_dir = f'{processed_data_path}/statistics'
os.makedirs(stats_dir, exist_ok=True)
fig_path = f'{data_repo_path}/figures/flow_cytometry/lymph_node_sizes_all_conditions/unpaired'
os.makedirs(fig_path, exist_ok=True)
mtc = True
non_parametric = False
force_parametric = True # set after initial normalicy testing
print(f'Using parametric unpaired with multiple testing correction: {mtc}')

df_total_counts = pd.read_csv(f'{processed_data_path}/early_gates_slo_counts_summed_over_panels.csv')
plot_CpG_conditions = False  # if False, drop the no CpG conditions. if true, plot only those conditions.
if plot_CpG_conditions:
    modifier = 'CpG_conditions'
    conditions_flow = ['ACT-d7_no CpG', 'ACT-d14_no CpG']
    df_total_counts = df_total_counts[df_total_counts['Condition'].str.contains('no CpG')]
    bar_colours_hex_blood_flow = bar_colours_hex_blood_flow[-3:-1] # only the last two colours
    x_axis_tick_labels_flow = ['ACT-d7\nno CpG', 'ACT-d14\nno CpG']
else:
    modifier = 'no_CpG_conditions'
    # drop condition if it contains the substring 'no CpG'
    df_total_counts = df_total_counts[~df_total_counts['Condition'].str.contains('no CpG')]
    conditions_flow = ['Naive', 'Tumour-bearing 3-5 mm', 'Cyclo only', 'ACT-d3', 'ACT-d7',
           'ACT-d14', 'ACT_Relapse']  # , 'ACT-d14_no CpG', 'ACT-d7_no CpG']
modifier = f'{modifier}_force_parametric_{force_parametric}_mtc_{mtc}_nonparam_{non_parametric}'
# reorder df_total_counts['Condition'] according to conditions_flow
df_total_counts['Condition'] = pd.Categorical(df_total_counts['Condition'], categories=conditions_flow, ordered=True)
# actually reorder the df
df_total_counts = df_total_counts.sort_values('Condition')

# lymph_nodes = ['inLNr', 'inLNl', 'brLNr']
value_cols = ['all_count_norm']
value_cols_y_label = ['Total cells\nper lymph node']
organ_col = 'Organ'
condition_col = 'Condition'
df_condition_names = df_total_counts[condition_col].unique()
conditions_to_compare = [(df_condition_names[i], df_condition_names[i + 1]) for i in
                         range(len(df_condition_names) - 1)]
ln_order = ['inLNr', 'inLNl', 'brLNr']
# %% df with the results of the statistical tests.
# stats_df = pd.DataFrame(columns=['Cell type', 'Condition', 'LN 1', 'LN 2', 'p-value', 'test_used'])
stats_df_columns = ['Cell type', organ_col, 'Condition 1', 'Condition 2', 'p-value', 'test_used', 'mtc_res']
stats_df = pd.DataFrame(columns=stats_df_columns)
for value_col in value_cols:
    for i_LN, LN in enumerate(ln_order):
        df_1ln = df_total_counts[df_total_counts[organ_col] == LN]
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
stats_df.to_csv(f'{stats_dir}/early_gates_stats_unpaired_stats_{modifier}.csv', index=False)
#%% plot the data
n_categories = len(df_total_counts['Condition'].unique())
ncols = len(ln_order)
fig_width = (n_categories-1)*ncols
markers_LN = [marker_per_organ_dict[ln] for ln in ln_order]

for cell_type_i, cell_type in enumerate(value_cols):
    fig, axes = plt.subplots(1, ncols, figsize=(fig_width, 6), sharey=True)
    axes = axes.flatten()
    # get the highest value for y-axis limit
    max_y = math.ceil(df_total_counts[cell_type].max())
    # for each lymph node in LN_type, plot the cell type count
    for i_LN, LN in enumerate(ln_order):
        ax = axes[i_LN]
        data = df_total_counts[(df_total_counts[organ_col] == LN)]  # [cell_type]
        # get stats df for cell type and LN
        stats_df_LN = stats_df[(stats_df['Cell type'] == cell_type) & (stats_df[organ_col] == LN)]
        stats_df_LN = stats_df_LN.drop(['Cell type', organ_col], axis=1)  # drop the LN_type column for the plot
        if mtc:
            p_val_col_name = 'mtc_res'
        else:
            p_val_col_name = 'p-value'
        boxplot_unpaired(data, condition_col, cell_type, stats_df=stats_df_LN, ax=ax, p_val_col_name=p_val_col_name,
                         strip_kwargs = {'palette':bar_colours_hex_blood_flow, 'marker':markers_LN[i_LN]}, y_lim=max_y)
        # ax.set_title(f'{lymph_nodes_label[i_LN]}', y=1.05)
        ax.set_title('')  # remove title
        ax.set_ylabel(f'{value_cols_y_label[cell_type_i]}')
        # ax.set_xlabel(organ_col)
        ax.set_xlabel('')  # remove xlabel
        ax.set_xticklabels(x_axis_tick_labels_flow)
        # set y to log
        ax.set_yscale('log')
        # set scale from 10^5 to 10^9
        ax.set_ylim(10**5, 10**9)
    plt.tight_layout()
    plt.savefig(f'{fig_path}/{cell_type}_flow_boxplot_early_gates_{modifier}_log.pdf')
    plt.close()
