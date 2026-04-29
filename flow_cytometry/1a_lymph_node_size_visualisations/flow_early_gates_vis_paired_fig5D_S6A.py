#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
# '''=================================================
# @Project: ACT_melanoma_neutrophil (data: act_neutrophil_data_repository)
# @Author : Gemma van der Voort
# @Desc   : use the preprocessed, summed data for fig 5D. do paired tests.
# @Desc updated: Paired early-gate lymph node size visualizations.
# update: us english, remove x-axis title and title
# update sep 2025: add noCPG conditions
# update 20260210: includes both options for no CpG and CpG conditions.
# '''=================================================
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

from paths_parameters import (x_axis_tick_labels_flow, bar_colours_hex_blood_flow, lymph_nodes_label,
                              marker_per_organ_dict, data_repo_path)
from functions_stat_test import ratio_paired_test
from functions_vis import boxplot_paired, single_paired_plot

plt.style.use(str(PROJECT_ROOT / "helper_files" / "journal_style.mplstyle"))

processed_data_path = f'{data_repo_path}/processed/flow_cytometry'
stats_dir = f'{processed_data_path}/statistics'
os.makedirs(stats_dir, exist_ok=True)
fig_path = f'{data_repo_path}/figures/flow_cytometry/lymph_node_sizes_all_conditions/paired'
os.makedirs(fig_path, exist_ok=True)
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

# reorder df_total_counts['Condition'] according to conditions_flow
df_total_counts['Condition'] = pd.Categorical(df_total_counts['Condition'], categories=conditions_flow, ordered=True)
# actually reorder the df
df_total_counts = df_total_counts.sort_values('Condition')

lymph_nodes = ['inLNr', 'inLNl', 'brLNr']
value_cols = ['all_count_norm']
value_cols_y_label = ['Total cells per lymph node']
group_col = 'Organ'
condition_col = 'Condition'
ln_order = ['inLNr', 'inLNl', 'brLNr']
LN_to_compare = [('inLNr', 'inLNl')]
# %% df with the results of the statistical tests.
stats_df = pd.DataFrame(columns=['Cell type', 'Condition', 'LN 1', 'LN 2', 'p-value', 'test_used'])
for value_col in value_cols:
    for condition in conditions_flow:
        for LN1, LN2 in LN_to_compare:
            # for some mice not both lymph nodes were processed. drop these mice
            ln_df_LN1 = df_total_counts[(df_total_counts[group_col] == LN1) & (df_total_counts[condition_col] == condition)]
            ln_df_LN2 = df_total_counts[(df_total_counts[group_col] == LN2) & (df_total_counts[condition_col] == condition)]
            # check if the mouse_ids are the same and in the same order
            if not  ln_df_LN1['Mouse_ID'].tolist() == ln_df_LN2['Mouse_ID'].tolist():
                # drop the mouse ID that is not in both organs
                only_ln1 = set(ln_df_LN1['Mouse_ID'].tolist()) - set(ln_df_LN2['Mouse_ID'].tolist())
                only_ln2 = set(ln_df_LN2['Mouse_ID'].tolist()) - set(ln_df_LN1['Mouse_ID'].tolist())
                print(f'Only in {LN1}: {only_ln1}, only in {LN2}: {only_ln2}. Dropped single mice')
                ln_df_LN1 = ln_df_LN1[~ln_df_LN1['Mouse_ID'].isin(only_ln1)]
                ln_df_LN2 = ln_df_LN2[~ln_df_LN2['Mouse_ID'].isin(only_ln2)]
            else:
                # print(f'Mice in {LN1} and {LN2} are the same')
                pass
            t_stat, p_val, test_used = ratio_paired_test(ln_df_LN1, ln_df_LN2, value_col, non_parametric=False)
            stats_df = pd.concat([stats_df, pd.DataFrame([[value_col, condition, LN1, LN2, p_val,
                                                           test_used]],
                                                         columns=['Cell type', 'Condition', 'LN 1', 'LN 2', 'p-value',
                                                                  'test_used'])])
#%% write to file
stats_df.to_csv(f'{stats_dir}/early_gates_stats_paired_stats_{modifier}.csv', index=False)
#%% plot the data
n_categories = len(conditions_flow)
ncols = 1
if plot_CpG_conditions:
    fig_width = ((n_categories*1.5)+1) * ncols
else:
    fig_width = ((n_categories*1.5)-1)*ncols  # a bit broader because we have 3 entries per condition
markers_LN = [marker_per_organ_dict[ln] for ln in ln_order]
for cell_type_i, cell_type in enumerate(value_cols):  # make single plots for each cell type
    print(f'Working on {cell_type}')
    for LN1, LN2 in LN_to_compare:  # plot for each LN combination, too crowed otherwise.
        fig, ax = plt.subplots(1, ncols , figsize=(fig_width, 6))
        data = df_total_counts  # [cell_type]
        # get stats df for cell type and LN combination
        stats_df_cell = stats_df[(stats_df['Cell type'] == cell_type) & (stats_df['LN 1'] == LN1) & (stats_df['LN 2']
                                                                                                     == LN2)]
        stats_df_cell = stats_df_cell.drop(['Cell type', 'Condition'], axis=1)  # drop the LN_type column for the plot
        boxplot_paired(data, 'Condition', cell_type, x_group_col=group_col, stats_df=stats_df_cell, ax=ax,
                       x_group_order=ln_order, markers=markers_LN, condition_colours=bar_colours_hex_blood_flow,
                       omit_ns=True
                       )
        # supported by stripplot?
        # add title with the LN combination
        # ax.set_title(f'Lymph nodes, paired test {LN1} and {LN2}', y=1.05)
        ax.set_ylabel(f'{value_cols_y_label[cell_type_i]}')
        ax.set_xlabel('')
        # y log scale
        ax.set_yscale('log')
        ax.set_xticklabels(x_axis_tick_labels_flow)
        plt.tight_layout()
        plt.savefig(f'{fig_path}/{cell_type}_{LN1}_{LN2}_ratio_paired_boxplot_log_{modifier}.pdf')
        plt.close()
