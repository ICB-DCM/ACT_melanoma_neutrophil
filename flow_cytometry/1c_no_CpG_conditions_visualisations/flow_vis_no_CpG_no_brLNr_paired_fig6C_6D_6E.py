#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
# '''=================================================
# @Project: ACT_melanoma_neutrophil (data: act_neutrophil_data_repository)
# @Author : Gemma van der Voort
# @Desc   : cd45 corrected data
# @Desc updated: Paired no-CpG analysis without brLNr for pan-immune and CD8 panels.
# update: now without d7. note: possibly a bit broad now, but I don't set any dimensions and it looks spacious, so
# keep for now.
# '''=================================================
import os
import sys
from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt
from pandas.api.types import CategoricalDtype
import numpy as np

PROJECT_ROOT = Path(__file__).resolve()
while PROJECT_ROOT.name != "ACT_melanoma_neutrophil" and PROJECT_ROOT.parent != PROJECT_ROOT:
    PROJECT_ROOT = PROJECT_ROOT.parent
if PROJECT_ROOT.name != "ACT_melanoma_neutrophil":
    raise RuntimeError("ACT_melanoma_neutrophil root not found")
sys.path.append(str(PROJECT_ROOT / "helper_files"))

from paths_parameters import x_axis_tick_labels_flow, bar_colours_hex_blood_flow, lymph_nodes_label, \
    marker_per_organ_dict, neutrophil_label, teff_label, tcm_label, data_repo_path
from functions_stat_test import ratio_paired_test
from functions_vis import single_paired_plot, boxplot_paired_no_brLNr

plt.style.use(str(PROJECT_ROOT / "helper_files" / "journal_style.mplstyle"))

# paths etc
processed_data_path_cd8 = f'{data_repo_path}/processed/flow_cytometry'
processed_data_path_panimmune = f'{data_repo_path}/processed/flow_cytometry'
stats_dir = f'{processed_data_path_panimmune}/statistics'
os.makedirs(stats_dir, exist_ok=True)
fig_path = f'{data_repo_path}/figures/flow_cytometry/no_cpg_conditions/paired'
df_panimmune = pd.read_csv(f'{processed_data_path_panimmune}/ln_pan_immune_normalised_v2_cd45_corrected.csv')
df_cd8 = pd.read_csv(f'{processed_data_path_cd8}/flow_cd8_panel_incl_tumor.csv')
modifier = 'cd45_corrected_no_d7'
# make dir and subdir single_paired if it doesn't exist
os.makedirs(fig_path, exist_ok=True)
#%%
lymph_nodes = ['inLNr', 'inLNl', 'brLNr']  # used for df subset check, so keep. adjusted later.
colours_conditions = [ "#fa815fff","#fa815fff"] # yellow, yellow, orange, orange
x_axis_tick_labels = ['ACT day 14',  'ACT day 14\nno CpG']
#%%
df_cd8 = df_cd8[df_cd8['Organ'].isin(lymph_nodes)]
#%% cd8 and panimmune have different amount of rows. count the amount each mouse_id occurs and print for both
df_cd8_counts = df_cd8['Mouse_ID'].value_counts()
df_panimmune_counts = df_panimmune['Mouse_ID'].value_counts()
# compare the counts per mouse_id
# print the mouse_ids that are not in both dataframes
only_in_cd8 = set(df_cd8_counts.index) - set(df_panimmune_counts.index)
only_in_panimmune = set(df_panimmune_counts.index) - set(df_cd8_counts.index)
print(f'Only in cd8: {only_in_cd8}')
print(f'Only in panimmune: {only_in_panimmune}')
# so there are a few missing, I assume for experimental reasons. Look into later. see obsidian.
#%% paired statistics for cpg-no cpg conditions for neutrophils
# keep only conditions containing d7 and d14
# df_panimmune_slim = df_panimmune[df_panimmune['Condition'].str.contains('d7|d14')]
df_panimmune_slim = df_panimmune[df_panimmune['Condition'].str.contains('d14')]
df_panimmune_slim = df_panimmune_slim.loc[:, ['Experiment', 'File', 'Condition','LN_type',
       'Mouse_ID', 'Neutrophils percent']]
#%% calculate the stats for the neutrophils
value_col = 'Neutrophils percent'
df_condition_names = df_panimmune_slim['Condition'].unique()
group_col = 'LN_type'
ln_order = ['inLNr', 'inLNl']
LN_to_compare = [('inLNr', 'inLNl')]
# df with the results of the statistical tests.
stats_df = pd.DataFrame(columns=['Cell type', 'Condition', 'LN 1', 'LN 2', 'p-value', 'test_used'])
for value_col in [value_col]:
    for condition in df_condition_names:
        for LN1, LN2 in LN_to_compare:
            # for some mice not both lymph nodes were processed. drop these mice
            ln_df_LN1 = df_panimmune_slim[(df_panimmune_slim[group_col] == LN1) & (df_panimmune_slim['Condition'] == condition)]
            ln_df_LN2 = df_panimmune_slim[(df_panimmune_slim[group_col] == LN2) & (df_panimmune_slim['Condition'] == condition)]
            # check if the mouse_ids are the same and in the same order
            if not  ln_df_LN1['Mouse_ID'].tolist() == ln_df_LN2['Mouse_ID'].tolist():
                # drop the mouse ID that is not in both conditions
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
# write to file
stats_df.to_csv(f'{stats_dir}/ln_pan_immune_stats_no_CpG_ratio_paired_{modifier}.csv')
#%% plot the data
n_categories = len(df_condition_names)
ncols = 1
fig_width = ((n_categories*1.5)+1)*ncols  # adapted for 2 categories
markers_LN = [marker_per_organ_dict[ln] for ln in ln_order]
#%% plot the data
for cell_type_i, cell_type in enumerate([value_col]):
    print('Working on', cell_type)
    for LN1, LN2 in LN_to_compare:
        # try:
        fig, ax = plt.subplots(1, ncols, figsize=(fig_width, 6))
        data = df_panimmune_slim
        # get stats df for cell type and LN combination
        stats_df_cell = stats_df[(stats_df['Cell type'] == cell_type) & (stats_df['LN 1'] == LN1) & (stats_df['LN 2']
                                                                                                 == LN2)]
        stats_df_cell = stats_df_cell.drop(['Cell type', 'Condition'],
                                           axis=1)  # drop the LN_type column for the plot
        # get the condition names
        condition_names = data['Condition'].unique()
        # get the values for the y axis
        y_values = []
        for condition in condition_names:
            y_values.append(data[data['Condition'] == condition][cell_type].values)
        # make a boxplot of the data
        boxplot_paired_no_brLNr(data, 'Condition', cell_type, x_group_col='LN_type', stats_df=stats_df_cell, ax=ax,
                       x_group_order=ln_order, markers=markers_LN, condition_colours=colours_conditions
                       )
        ax.set_ylabel(neutrophil_label)
        # ax.set_xlabel('Condition')
        ax.set_xlabel('')
        ax.set_xticklabels(x_axis_tick_labels)
        # save the figure
        plt.savefig(f'{fig_path}/{cell_type}_{LN1}_{LN2}_notitle_{modifier}.pdf', dpi=300)
        plt.close(fig)
#%% same for tcm and teff
# df_cd8 = df_cd8[df_cd8['Condition'].str.contains('d7|d14')]
df_cd8 = df_cd8[df_cd8['Condition'].str.contains('d14')]
# re-order the conditions according to order in df_condition_names
cond_order = ['ACT-d14', 'ACT-d14_no CpG']
cat_dtype = CategoricalDtype(categories=cond_order, ordered=True)
df_cd8['Condition'] = df_cd8['Condition'].astype(cat_dtype)
df_cd8 = df_cd8.sort_values('Condition')

value_cols = ['Teff_freq_pmel', 'Tcm_freq_pmel']
label_names = [teff_label, tcm_label]
group_col = 'Organ'
ln_order = ['inLNr', 'inLNl']
LN_to_compare = [('inLNr', 'inLNl')]
# df with the results of the statistical tests.
stats_df = pd.DataFrame(columns=['Cell type', 'Condition', 'LN 1', 'LN 2', 'p-value', 'test_used'])
for value_col in value_cols:
    for condition in df_condition_names:
        for LN1, LN2 in LN_to_compare:
            # for some mice not both lymph nodes were processed. drop these mice
            ln_df_LN1 = df_cd8[(df_cd8[group_col] == LN1) & (df_cd8['Condition'] == condition)]
            ln_df_LN2 = df_cd8[(df_cd8[group_col] == LN2) & (df_cd8['Condition'] == condition)]
            # check if the mouse_ids are the same and in the same order
            if not  ln_df_LN1['Mouse_ID'].tolist() == ln_df_LN2['Mouse_ID'].tolist():
                # drop the mouse ID that is not in both conditions
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
# write to file
stats_df.to_csv(f'{stats_dir}/ln_cd8_tcm_teff_stats_no_CpG_CpG_ratio_paired_{modifier}.csv')
#%% plot the data
ncols = 1
fig_width = ((n_categories*1.5)+1)*ncols  # adapted for 2 categories
markers_LN = [marker_per_organ_dict[ln] for ln in ln_order]
#%% plot the data
for cell_type_i, cell_type in enumerate(value_cols):
    print('Working on', cell_type)
    for LN1, LN2 in LN_to_compare:
        # try:
        fig, ax = plt.subplots(1, ncols, figsize=(fig_width, 6))
        data = df_cd8
        # get stats df for cell type and LN combination
        stats_df_cell = stats_df[(stats_df['Cell type'] == cell_type) & (stats_df['LN 1'] == LN1) & (stats_df['LN 2']
                                                                                                 == LN2)]
        stats_df_cell = stats_df_cell.drop(['Cell type', 'Condition'], axis=1)  # drop the LN_type column for the plot

        # get the values for the y axis
        y_values = []
        for condition in df_condition_names:
            y_values.append(data[data['Condition'] == condition][cell_type].values)
        # make a boxplot of the data
        boxplot_paired_no_brLNr(data, 'Condition', cell_type, x_group_col='Organ', stats_df=stats_df_cell, ax=ax,
                       x_group_order=ln_order, markers=markers_LN, condition_colours=colours_conditions,

                       )
        ax.set_ylabel(f'{label_names[cell_type_i]}')
        # ax.set_xlabel('Condition')
        ax.set_xlabel('')
        ax.set_xticklabels(x_axis_tick_labels)
        # save the figure
        plt.savefig(f'{fig_path}/{cell_type}_{LN1}_{LN2}_notitle_{modifier}.pdf', dpi=300)
