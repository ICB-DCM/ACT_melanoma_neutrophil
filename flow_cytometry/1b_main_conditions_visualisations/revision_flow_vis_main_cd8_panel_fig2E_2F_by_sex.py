#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
# '''=================================================
# @Project: ACT_melanoma_neutrophil (data: act_neutrophil_data_repository)
# @Author : Gemma van der Voort
# @Desc   : paired plot but without the brLNr for less crowding, plus paired single plots for Tcm and Teff day 14.
# @Desc updated: Main-condition CD8 panel paired plots and single paired highlights.
# update: add counts
# '''=================================================
import os
import sys
from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt

PROJECT_ROOT = Path(__file__).resolve()
while PROJECT_ROOT.name != "ACT_melanoma_neutrophil" and PROJECT_ROOT.parent != PROJECT_ROOT:
    PROJECT_ROOT = PROJECT_ROOT.parent
if PROJECT_ROOT.name != "ACT_melanoma_neutrophil":
    raise RuntimeError("ACT_melanoma_neutrophil root not found")
sys.path.append(str(PROJECT_ROOT / "helper_files"))

from paths_parameters import x_axis_tick_labels_flow, bar_colours_hex_blood_flow, lymph_nodes_label, \
    marker_per_organ_dict, teff_label, tcm_label, data_repo_path
from functions_stat_test import ratio_paired_test
from functions_vis import boxplot_paired_no_brLNr, single_paired_plot

plt.style.use(str(PROJECT_ROOT / "helper_files" / "journal_style.mplstyle"))

# paths etc
processed_data_path = f'{data_repo_path}/processed/flow_cytometry'
stats_dir = f'{processed_data_path}/statistics'
os.makedirs(stats_dir, exist_ok=True)
fig_path = f'{data_repo_path}/figures/flow_cytometry/main_conditions'

# make dir and sub-dirs if they don't exist
os.makedirs(fig_path, exist_ok=True)
os.makedirs(f'{fig_path}/multiple_paired', exist_ok=True)
os.makedirs(f'{fig_path}/single_paired', exist_ok=True)
females_rm = True

df_complete = pd.read_csv(f'{processed_data_path}/flow_cd8_panel_incl_tumor.csv', index_col=0)
# df_condition_names = df_complete['Condition'].unique()  # since CD8 pmels, remove all non-ACT conditions
df_condition_names = ['ACT-d3', 'ACT-d7', 'ACT-d14', 'ACT_Relapse', 'ACT-d14_no CpG', 'ACT-d7_no CpG']
#%% add sex info
d7_female_mice = ['HOE-2390', 'HOE-2391', 'HOE-2388', 'HOE-2389', 'HOE-2387']
d14_female_mice = ['HOE-2384', 'HOE-2385', 'HOE-2386', 'HOE-2393']
female_mice_ids = d7_female_mice + d14_female_mice
#%% version 1: freq pmel only
value_cols_pmel = ['Teff_freq_pmel',
                   'Tcm_freq_pmel']
value_cols = value_cols_pmel
value_cols_y_label = ['pmel-1 $T_{EFF}$ \nFrequency of pmel-1 cells',
                      'pmel-1 $T_{CM}$ \nFrequency of pmel-1 cells']
if females_rm:
    img_modifier='incl_102_d14_pmelfreq_only_fem_removed'  # incl_102_d14
    df_complete = df_complete[~df_complete['Mouse_ID'].isin(female_mice_ids)].copy()
else:
    img_modifier='incl_102_d14_pmelfreq_only'
# for single paired plots, also check if it needs adjustment below.
# %% version 2: counts
# value_cols_pmel = ['Teff_count',
#                    'Tcm_count']
# value_cols = value_cols_pmel
# value_cols_y_label = ['pmel-1 $T_{EFF}$ \nCount normalized',
#                       'pmel-1 $T_{CM}$ \nCount normalized']
# img_modifier='incl_102_d14_counts_only'

group_col = 'Organ'
condition_col = 'Condition'
LN_to_compare = [('inLNr', 'inLNl')]
ln_order = ['inLNr', 'inLNl']
# subset data to pmels only, ln only
df_pmel = df_complete[df_complete['Parent'] == 'pmel']
df_pmel = df_pmel[df_pmel['Organ'].isin(ln_order)]
df = df_pmel.copy()
# print the amount of 0 values
zero_counts = df[value_cols].eq(0).sum()
stats_df = pd.DataFrame(columns=['Cell type', 'Condition', 'LN 1', 'LN 2', 'p-value', 'test_used'])
for value_col in value_cols:
    for condition in df_condition_names:
        for LN1, LN2 in LN_to_compare:
            print(LN1)
            print(LN2)
            # for some mice not both lymph nodes were processed. drop these mice
            ln_df_LN1 = df[(df[group_col] == LN1) & (df['Condition'] == condition)]
            ln_df_LN2 = df[(df[group_col] == LN2) & (df['Condition'] == condition)]
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
stats_df.to_csv(f'{stats_dir}/paired_stats_cd8_panel_pmel_no_brLNr_{img_modifier}.csv')
# todo d3 and relapse conditions have some values at 0, which invalidates the ratio paired test. add small offset.
# I don't show this data in the paper, but leave this for revisions if needed.
#%% single paired plots  (as example/highlight)
condition = 'ACT-d14'
# for counts:
cell_types_single = value_cols_pmel
names_single = value_cols_y_label

y_lims = []
df_d14 = df[df['Condition'] == condition]
# keep lns only (otherwise y_highest_point is incorrect
df_d14 = df_d14[df_d14['Organ'].isin(['inLNr', 'inLNl'])]
for cell_type_i, cell_type_single in enumerate(cell_types_single):
    y_lim = df_d14[cell_type_single].max()
    y_highest_point = df[cell_type_single].max()
    print(f'Highest point in {cell_type_single} is {y_highest_point}')
    y_lims.append(y_lim)
y_lim = max(y_lims)
for cell_type_i, cell_type_single in enumerate(cell_types_single):
    fig, ax = plt.subplots(1, 1, figsize=(3, 6))
    colour = "#fa815fff"
    markers = ['o', 's']
    x_pair = ('inLNr', 'inLNl')
    single_paired_plot(df_d14, 'Condition', cell_type_single, 'Mouse_ID', condition, group_col,
                       x_pair, stats_df, markers = markers, scatter_kwargs={'color': colour, 's': 50}, ax=ax,  y_lim=y_lim,
                       remove_unpaired=False, mouse_id_female=female_mice_ids)
    ax.set_ylabel(f'{names_single[cell_type_i]}', fontsize=20)
    ax.set_xlabel("")
    ax.tick_params(axis='x', labelsize=20)
    plt.savefig(f'{fig_path}/single_paired/{cell_type_single}_{condition}_ratio_paired_counts_{img_modifier}_sex.pdf')
# %% not all female mouse IDs are printed. Print all mouse IDs in df_d14
print(f'Mouse IDs in {condition}: {df_d14["Mouse_ID"].unique()}')
#save to file
df_complete['Mouse_ID'].value_counts().sort_index().to_csv(f'{processed_data_path}/mouse_id_counts_complete_{img_modifier}.csv')
