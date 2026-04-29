#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
# '''=================================================
# @Project: ACT_melanoma_neutrophil (data: act_neutrophil_data_repository)
# @Author : Gemma van der Voort
# @Desc   : replotting of the large region data from the codex data. paired on the anndata object.
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

from paths_parameters import x_axis_tick_labels_codex, bar_colours_hex_codex, marker_per_organ_dict, data_repo_path
from functions_vis import boxplot_paired, boxplot_paired_no_brLNr
from functions_stat_test import ratio_paired_test

plt.style.use(str(PROJECT_ROOT / "helper_files" / "journal_style_old.mplstyle"))

# paths etc
processed_data_path = f'{data_repo_path}/processed/codex'
stats_dir = f'{processed_data_path}/statistics'
os.makedirs(stats_dir, exist_ok=True)
metacluster_key =  'Metacluster v4_8_filtered_0.5_maxit_25_v4_8' # updated in revision
cn_col = 'CN_k50_n20'
neut_df_file = f"Neutrophil_abundance_3regions_{cn_col}_{metacluster_key}_filtered.csv"

fig_dir = f'{data_repo_path}/figures/codex/manuscipt_plots/large_regions'

norm_method = 'percent_z' # norm by cd45+ cells in that metacluster
print(f'no stats, data normalised by {norm_method}')
modifier = f'no_stats_{norm_method}_{metacluster_key}'
# make dir and subdir single_paired if it doesn't exist
os.makedirs(f'{fig_dir}', exist_ok=True)
#%%
neut_df = pd.read_csv(f'{processed_data_path}/{neut_df_file}')
value_cols_part = ['B-cell follicle', 'T-cell zone', 'Medulla-Interfollicular zone-SCS']
value_cols = [f'{col}_{norm_method}' for col in value_cols_part]
value_cols_title = ['B cell follicle', 'T cell zone', 'Medulla/SCS/Interfollicular zone']
y_label = 'Neutrophil\n% of cells in zone'
x_label = ''

condition_column = 'Treatment'
mouse_id_column = 'mouse_id'

# reorder the conditions
condition_order = ['tumour_bearing', 'Cyclo_d1', 'ACT_d3', 'ACT_d7', 'ACT_d14']
neut_df[condition_column] = pd.Categorical(neut_df[condition_column], categories=condition_order)
# also actually reshuffle the df to have the conditions in the right order
neut_df = neut_df.sort_values(by=[condition_column])
df_condition_names = neut_df[condition_column].unique()
group_col = 'Organ'
ln_order = ['inLNr', 'inLNl', 'brLNr']
LN_to_compare = [('inLNr', 'inLNl')]

# df with the results of the statistical tests.
stats_df = pd.DataFrame(columns=['Region', 'Condition', 'LN 1', 'LN 2', 'p-value', 'test_used'])
for value_col in value_cols:
    for condition in condition_order:
        for LN1, LN2 in LN_to_compare:
            # for some mice not both lymph nodes were processed. drop these mice
            ln_df_LN1 = neut_df[(neut_df[group_col] == LN1) & (neut_df[condition_column] == condition)]
            ln_df_LN2 = neut_df[(neut_df[group_col] == LN2) & (neut_df[condition_column] == condition)]
            # check if the mouse_ids are the same and in the same order
            if not  ln_df_LN1[mouse_id_column].tolist() == ln_df_LN2[mouse_id_column].tolist():
                # drop the mouse ID that is not in both conditions
                only_ln1 = set(ln_df_LN1[mouse_id_column].tolist()) - set(ln_df_LN2[mouse_id_column].tolist())
                only_ln2 = set(ln_df_LN2[mouse_id_column].tolist()) - set(ln_df_LN1[mouse_id_column].tolist())
                print(f'Only in {LN1}: {only_ln1}, only in {LN2}: {only_ln2}. Dropped single mice')
                ln_df_LN1 = ln_df_LN1[~ln_df_LN1[mouse_id_column].isin(only_ln1)]
                ln_df_LN2 = ln_df_LN2[~ln_df_LN2[mouse_id_column].isin(only_ln2)]
            else:
                # print(f'Mice in {LN1} and {LN2} are the same')
                pass
            # if ln_df_LN1 is less than 3 rows, skip the test
            if len(ln_df_LN1) < 3 or len(ln_df_LN2) < 3:
                print(f'Skipping {value_col} for condition {condition} in {LN1} and {LN2} because not enough data')
                continue
            else:
                print(f'amount of mice in {LN1} and {LN2} for condition {condition}: {len(ln_df_LN1)} and {len(ln_df_LN2)}')
            t_stat, p_val, test_used = ratio_paired_test(ln_df_LN1, ln_df_LN2, value_col, non_parametric=False)
            stats_df = pd.concat([stats_df, pd.DataFrame([[value_col, condition, LN1, LN2, p_val,
                                                           test_used]],
                                                         columns=['Region', 'Condition', 'LN 1', 'LN 2', 'p-value',
                                                                  'test_used'])])
# write to file
stats_df.to_csv(f'{stats_dir}/codex_regions_stats_{modifier}.csv')

# %% plot the data
n_categories = len(condition_order)
ncols = 1
fig_width = ((n_categories*1.5)-1)*ncols  # a bit broader because we have 3 entries per condition
markers_LN = [marker_per_organ_dict[ln] for ln in ln_order]
for region_type_i, region_type in enumerate(value_cols):
    print('Working on', region_type)
    for LN1, LN2 in LN_to_compare:
        # try:
        fig, ax = plt.subplots(1, ncols)
        stats_df_region = stats_df[(stats_df['Region'] == region_type) & (stats_df['LN 1'] == LN1) & (stats_df['LN 2']
                                                                                                 == LN2)]
        stats_df_region = stats_df_region.drop(['Region', 'Condition'], axis=1)
        boxplot_paired(neut_df, condition_column, region_type, x_group_col='Organ', stats_df=None, ax=ax, y_lim=16,
                       x_group_order=ln_order, markers=markers_LN, condition_colours=bar_colours_hex_codex
                       )
        ax.set_title(value_cols_title[region_type_i], y=1.05)
        ax.set_ylabel(y_label)
        ax.set_xlabel(x_label)
        ax.set_xticklabels(x_axis_tick_labels_codex)
        plt.tight_layout()
        plt.savefig(f'{fig_dir}/{region_type}_{modifier}_stats.pdf', dpi=300)
        plt.close()
# %%repeat the plot but for two lymph nodes only
ln_order = ['inLNr', 'inLNl']
markers_LN = [marker_per_organ_dict[ln] for ln in ln_order]
neut_df_two = neut_df[neut_df['Organ'].isin(ln_order)]
for region_type_i, region_type in enumerate(value_cols):
    print('Working on', region_type)
    for LN1, LN2 in LN_to_compare:
        # try:
        fig, ax = plt.subplots(1, ncols)
        stats_df_region = stats_df[(stats_df['Region'] == region_type) & (stats_df['LN 1'] == LN1) & (stats_df['LN 2']
                                                                                                 == LN2)]
        stats_df_region = stats_df_region.drop(['Region', 'Condition'], axis=1)
        boxplot_paired_no_brLNr(neut_df_two, condition_column, region_type, x_group_col='Organ', stats_df=None, ax=ax, y_lim=16,
                       x_group_order=ln_order, markers=markers_LN, condition_colours=bar_colours_hex_codex
                       )
        ax.set_title(value_cols_title[region_type_i], y=1.05)
        ax.set_ylabel(y_label)
        ax.set_xlabel(x_label)
        ax.set_xticklabels(x_axis_tick_labels_codex)
        plt.tight_layout()
        plt.savefig(f'{fig_dir}/{region_type}_{modifier}_nobrLNr_stats_replotted.pdf', dpi=300)
        plt.close()
