#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
# '''=================================================
# @Project: ACT_melanoma_neutrophil (data: act_neutrophil_data_repository)
# @Author : Gemma van der Voort
# @Desc   : update august: parametric test, remove x axis title
# @Desc updated: Unpaired Pmel-in-blood visualization and stats.
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

from paths_parameters import x_axis_tick_labels_flow, bar_colours_hex_blood_flow, data_repo_path
from functions_stat_test import unpaired_multiple_conditions # unpaired_test
from functions_vis import boxplot_unpaired

plt.style.use(str(PROJECT_ROOT / "helper_files" / "journal_style.mplstyle"))
# paths etc
data_output_path = f'{data_repo_path}/processed/flow_cytometry'
stats_dir = f'{data_output_path}/statistics'
os.makedirs(stats_dir, exist_ok=True)
fig_path = f'{data_repo_path}/figures/flow_cytometry/main_conditions'
os.makedirs(fig_path, exist_ok=True)

df_flow_blood = pd.read_csv(f'{data_output_path}/pmel_t_blood_means.csv', index_col=0)
mtc = True
force_parametric = True  # set after initial normalicy testing
modifier = 'nullcond_rem_mtc_parameteric'
# rename cohort to condition
df_flow_blood = df_flow_blood.rename(columns={'Cohort': 'Condition', 'Mean.2': 'Mean'})
group_col = 'Condition'
value_col = 'Mean'
value_col_file_name = 'pmel_mean_blood'
y_label_name = 'Pmel-1 CD8$^+$ T cells\n% of live CD45.2$^+$ cells'
# drop all columns with Conditions 'Naive', 'Tumour-bearing 3-5 mm', 'Cyclo only'
df_flow_blood = df_flow_blood[~df_flow_blood['Condition'].isin(['Naive', 'Tumour-bearing 3-5 mm', 'Cyclo only'])]

#
# stat test
df_condition_names = ['ACT-d3', 'ACT-d7', 'ACT-d14', 'ACT_Relapse']
conditions_to_compare = [(df_condition_names[i], df_condition_names[i + 1]) for i in
                         range(len(df_condition_names) - 1)]
# df with the results of the statistical tests.
stats_df_columns = ['Cell type', 'Condition 1', 'Condition 2', 'p-value', 'test_used', 'mtc_res']
stats_df = pd.DataFrame(columns=[stats_df_columns])

# %%for each condition, compare it to the next one
for condition1, condition2 in conditions_to_compare:
    print(condition1, condition2)
    t_stat, p_val, test_used, mtc_res = unpaired_multiple_conditions(df_flow_blood, value_col, group_col,
                                             (condition1, condition2), force_parametric=True,
                                                                     non_parametric=False,
                                                            multiple_testing_correction=True)
    print('pval', p_val)
    # print('t stat', t_stat)
    print('mtc res', mtc_res)
    stats_df = pd.concat([stats_df, pd.DataFrame([[value_col, condition1, condition2, p_val, test_used, mtc_res]],
                                                 columns=[stats_df_columns])])
# write to file
stats_df.to_csv(f'{stats_dir}/pmel_flow_blood_stats_{modifier}.csv')

#%% plot the data
bar_colours_hex_blood_flow_pmel = bar_colours_hex_blood_flow[3:]
x_axis_tick_labels_blood_pmel = x_axis_tick_labels_flow[3:]
n_categories = len(df_flow_blood[group_col].unique())
fig, ax = plt.subplots(1, 1, figsize=(n_categories, 6))
df_flow_blood = df_flow_blood.reset_index(drop=True)
# replace nans with 0
df_flow_blood = df_flow_blood.fillna(0)
# drop the CpG conditions
df_flow_blood = df_flow_blood[~df_flow_blood['Condition'].str.contains('CpG')]
boxplot_unpaired(df_flow_blood, group_col, value_col, stats_df=None, ax=ax, p_val_col_name=None,
                 strip_kwargs={'palette': bar_colours_hex_blood_flow_pmel, 'marker': 'X'})
ax.set_title('')  # remove title
ax.set_ylabel(y_label_name)
ax.set_xlabel(group_col)
ax.set_xticklabels(x_axis_tick_labels_blood_pmel)
ax.set_xlabel('')  # remove xlabel
plt.tight_layout()
plt.savefig(f'{fig_path}/{value_col_file_name}_boxplot_no_title_{modifier}.pdf')
