#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
# '''=================================================
# @Project: ACT_melanoma_neutrophil (data: act_neutrophil_data_repository)
# @Author : Gemma van der Voort
# @Desc   : new visualisation for the blood data of the CD8+ melanoma project. plots neut and lymphocytes side by side
# @Desc updated: Hematology blood plots for neutrophils and lymphocytes.
# update 20 feb: take out tumor onset time point, do mtc
# '''=================================================
import sys
from pathlib import Path
import os
import pandas as pd
import matplotlib.pyplot as plt

PROJECT_ROOT = Path(__file__).resolve()
while PROJECT_ROOT.name != "ACT_melanoma_neutrophil" and PROJECT_ROOT.parent != PROJECT_ROOT:
    PROJECT_ROOT = PROJECT_ROOT.parent
if PROJECT_ROOT.name != "ACT_melanoma_neutrophil":
    raise RuntimeError("ACT_melanoma_neutrophil root not found")
sys.path.append(str(PROJECT_ROOT / "helper_files"))

from paths_parameters import x_axis_tick_labels_blood, bar_colours_hex_blood_flow, data_repo_path
from functions_stat_test import unpaired_multiple_conditions # unpaired_test
from functions_vis import boxplot_unpaired

plt.style.use(str(PROJECT_ROOT / "helper_files" / "journal_style.mplstyle"))
# paths etc
data_output_path = f'{data_repo_path}/processed/hematology'
stats_dir = f'{data_output_path}/statistics'
os.makedirs(stats_dir, exist_ok=True)
fig_path = f'{data_repo_path}/figures/hematology'
modifier_file = 'data_corrections_jan2026'

df_blood = pd.read_csv(f'{data_output_path}/neut_lym_stacked_{modifier_file}.csv', index_col=0)
# drop condition tumor onset
df_blood = df_blood[df_blood['Condition'] != 'Onset tumour growth']
mtc = True
force_parametric = True # set after initial normalicy testing
modifier = 'no_onset_mtc_parametric'

# stat test
df_condition_names = df_blood['Condition'].unique()
group_col = 'Condition'
value_cols = ['Neutrophil count (10^9/L)', 'Lymphocyte count (10^9/L)']
value_col_file_name = ['Neutrophil_count', 'Lymphocyte_count']
# modify so that ^9 is superscript
y_label_names = ['Neutrophils in blood [10$^9$/L]', 'Lymphocytes in blood [10$^9$/L]']
conditions_to_compare = [(df_condition_names[i], df_condition_names[i + 1]) for i in
                         range(len(df_condition_names) - 1)]
# df with the results of the statistical tests.
stats_df_columns = ['Cell type', 'Condition 1', 'Condition 2', 'p-value', 'test_used', 'mtc_res']
stats_df = pd.DataFrame(columns=stats_df_columns)
for value_col in value_cols:
    # for each condition, compare it to the next one
    for condition1, condition2 in conditions_to_compare:
        t_stat, p_val, test_used, mtc_res = unpaired_multiple_conditions(df_blood, value_col, group_col,
                                                 (condition1, condition2), force_parametric=True,
                                                                         non_parametric=False,
                                                                multiple_testing_correction=True)
        stats_df = pd.concat([stats_df, pd.DataFrame([[value_col, condition1, condition2, p_val, test_used, mtc_res]],
                                                     columns=stats_df_columns)])
# write to file
stats_df.to_csv(f'{stats_dir}/neut_lym_blood_stats_{modifier_file}_{modifier}.csv')
#%%
stats_df = stats_df.reset_index(drop=True)

#%% plot the data, one plot per value_col instead of side by side
n_categories = len(df_blood['Condition'].unique())
for value_col_i, value_col in enumerate(value_cols):
    print(f'Working on {value_col}')
    fig, ax = plt.subplots(1, 1, figsize=(n_categories, 6))
    stats_df_cell = stats_df[stats_df['Cell type'] == value_col].copy()
    print(f'Stats df for {value_col}:\n{stats_df_cell}')
    stats_df_cell = stats_df_cell.drop('Cell type', axis=1).copy()  # drop the cell type column for the plot
    if mtc:
        p_val_col_name = 'mtc_res'
    else:
        p_val_col_name = 'p-value'
    boxplot_unpaired(df_blood, group_col, value_col, stats_df=stats_df_cell, ax=ax, p_val_col_name=p_val_col_name,
                     strip_kwargs={'palette': bar_colours_hex_blood_flow, 'marker': 'X'})
    ax.set_title('')  # remove title
    ax.set_ylabel(y_label_names[value_col_i])
    ax.set_xlabel(group_col)
    ax.set_xticklabels(x_axis_tick_labels_blood)
    ########### modifications from style sheet for these plots only ###########
    # set x tick label font size smaller
    plt.setp(ax.get_xticklabels(), fontsize=16)
    # omit the highest 2 y tick labels (2, since the highest one is also in the list but not printed)
    # also omit the lowest label (-1, not biologically relevant)
    yticks = ax.get_yticks().tolist()
    if len(yticks) > 4:
        ax.set_yticks(yticks[1:-2])
    # remove '.0' from y tick labels if present (integer when possible)
    yticklabels = [str(int(tick)) if tick.is_integer() else str(tick) for tick in ax.get_yticks()]
    ax.set_yticklabels(yticklabels)
    ######## back to base style ###########
    ax.set_xlabel('')
    plt.tight_layout()
    plt.savefig(f'{fig_path}/{value_col_file_name[value_col_i]}_blood_boxplot_no_title_{modifier_file}_{modifier}.pdf')
