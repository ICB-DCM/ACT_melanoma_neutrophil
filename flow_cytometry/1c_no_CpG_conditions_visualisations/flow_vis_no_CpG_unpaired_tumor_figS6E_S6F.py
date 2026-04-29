#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
# '''=================================================
# @Project: ACT_melanoma_neutrophil (data: act_neutrophil_data_repository)
# @Author : Gemma van der Voort
# @Desc   : perturbation experiment, unpaired (between conditions). Needed some adjustment of plotting code for
# significance brackets. this one is for tumor flow data (only cd8 panel)
# @Desc updated: Unpaired CD8 tumor plots for no-CpG conditions.
# update: parameteric tests, for sup7
# update 20251023: rm d7 from sup7
# '''=================================================
import os
import sys
from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt
from pandas.api.types import CategoricalDtype
import math

PROJECT_ROOT = Path(__file__).resolve()
while PROJECT_ROOT.name != "ACT_melanoma_neutrophil" and PROJECT_ROOT.parent != PROJECT_ROOT:
    PROJECT_ROOT = PROJECT_ROOT.parent
if PROJECT_ROOT.name != "ACT_melanoma_neutrophil":
    raise RuntimeError("ACT_melanoma_neutrophil root not found")
sys.path.append(str(PROJECT_ROOT / "helper_files"))

from paths_parameters import x_axis_tick_labels_flow, bar_colours_hex_blood_flow, organ_names_short, \
    marker_per_organ_dict, teff_label, tcm_label, data_repo_path
from functions_stat_test import unpaired_multiple_conditions  # unpaired_test
from functions_vis import boxplot_unpaired

plt.style.use(str(PROJECT_ROOT / "helper_files" / "journal_style.mplstyle"))

# paths etc
processed_data_path_cd8 = f'{data_repo_path}/processed/flow_cytometry'
stats_dir = f'{processed_data_path_cd8}/statistics'
os.makedirs(stats_dir, exist_ok=True)
fig_path = f'{data_repo_path}/figures/flow_cytometry/no_cpg_conditions/unpaired'
df_cd8 = pd.read_csv(f'{processed_data_path_cd8}/flow_cd8_panel_incl_tumor.csv')
test_function_used = 'unpaired_multiple_conditions'
modifier= 'tumor_no_d7cpg'
mtc = True
non_parametric = False  # decided on parametric after testing
force_parametric = True # set after initial normalicy testing
print(f'Using {test_function_used} with multiple testing correction: {mtc}')
# make dir and subdir single_paired if it doesn't exist
os.makedirs(fig_path, exist_ok=True)
#%%
lymph_nodes = ['Tumour']
colours_conditions = [ "#fa815fff","#fa815fff"]  # orange, orange
x_axis_tick_labels = ['ACT\nday 14',  'ACT\nday 14\nno CpG']  # preserve order!
conditions = ['ACT-d14', 'ACT-d14_no CpG']
conditions_to_compare = [ ('ACT-d14', 'ACT-d14_no CpG')]
group_col = 'Condition'
#%% repeat for tcm and teff
df_cd8_slim = df_cd8[['Mouse_ID', 'Condition', 'Organ','Teff_freq_pmel', 'Tcm_freq_pmel']]
df_cd8_slim = df_cd8_slim[df_cd8_slim['Condition'].isin(conditions)]
df_cd8_slim = df_cd8_slim[df_cd8_slim['Organ'].isin(lymph_nodes)]
#%%
value_cols = ['Teff_freq_pmel', 'Tcm_freq_pmel']
value_cols_y_label = [teff_label, tcm_label]
organ_col = 'Organ'
cat_dtype = CategoricalDtype(categories=conditions, ordered=True)
#%%
df_cd8_slim['Condition'] = df_cd8_slim['Condition'].astype(cat_dtype)
df_cd8_slim = df_cd8_slim.sort_values('Condition')
# %%stats
stats_df_columns = ['Cell type', organ_col, 'Condition 1', 'Condition 2', 'p-value', 'test_used', 'mtc_res']
stats_df_cd8 = pd.DataFrame(columns=stats_df_columns)
for value_col in value_cols:
    for i_LN, LN in enumerate(lymph_nodes):
        # for each condition, compare it to the next one
        for condition1, condition2 in conditions_to_compare:
            df_1ln = df_cd8_slim[df_cd8_slim[organ_col] == LN]
            t_stat, p_val, test_used, mtc_res = unpaired_multiple_conditions(df_1ln, value_col, group_col,
                                                                             (condition1, condition2),
                                                                             non_parametric=non_parametric,
                                                                             multiple_testing_correction=mtc,
                                                                             force_parametric=force_parametric
                                                                             )
            stats_df_cd8 = pd.concat([stats_df_cd8, pd.DataFrame([[value_col, LN, condition1, condition2, p_val,
                                                           test_used, mtc_res]],
                                                         columns=stats_df_columns)])
# write to file
stats_df_cd8.to_csv(f'{stats_dir}/ln_cd8_stats_CpG_perturbations_{test_function_used}_mtc{mtc}_'
                f'np_{non_parametric}_fp_{force_parametric}_{modifier}.csv')
#%% plotting
n_categories = len(conditions)
ncols = len(lymph_nodes)
width_per_plot = (n_categories-1)*(1+n_categories*0.35)
axes_width = ncols * width_per_plot
fixed_margin = 1.5  # space for y-axis and labels
fig_width = axes_width + fixed_margin
markers_LN = 'P'
for cell_type_i, cell_type in enumerate(value_cols):
    fig, axes = plt.subplots(1, ncols, figsize=(4, 6), sharey=True)

    # Set margins manually for equal side independent of subplot amount
    left_margin_fraction = fixed_margin / fig_width
    fig.subplots_adjust(left=left_margin_fraction, right=0.98, wspace=0.3)

    # if more than one axis, flatten the axes
    if ncols > 1:
        axes = axes.flatten()
    else:
        axes = [axes]
    # get the highest value for y-axis limit
    max_y = math.ceil(df_cd8_slim[cell_type].max())
    # for each lymph node in LN_type, plot the cell type count
    for i_LN, LN in enumerate(lymph_nodes):
        ax = axes[i_LN]
        data = df_cd8_slim[(df_cd8_slim[organ_col] == LN)]  # [cell_type]
        # get stats df for cell type and LN
        stats_df_LN = stats_df_cd8[(stats_df_cd8['Cell type'] == cell_type) & (stats_df_cd8[organ_col] == LN)]
        stats_df_LN = stats_df_LN.drop(['Cell type', organ_col], axis=1)  # drop the LN_type column for the plot
        print(stats_df_LN.head())
        if mtc:
            p_val_col_name = 'mtc_res'
        else:
            p_val_col_name = 'p-value'
        boxplot_unpaired(data, group_col, cell_type, stats_df=stats_df_LN, ax=ax, p_val_col_name=p_val_col_name,
                         strip_kwargs = {'palette':colours_conditions, 'marker':markers_LN[i_LN]}, y_lim=max_y,
                         significance_bar_order=['low', 'low', 'high', 'higher'])
        ax.set_title('')  # remove title
        ax.set_ylabel(f'{value_cols_y_label[cell_type_i]}')
        ax.set_xlabel("")
        # set y tick labels (to prevent >100 ticks showing)
        ax.set_yticks([0, 20, 40, 60, 80, 100])
        ax.set_xticklabels(x_axis_tick_labels)  #, rotation=90)
    # plt.tight_layout()
    plt.savefig(f'{fig_path}/{cell_type}_flow_boxplot_CpG_perturbations_{test_used}_'
                f'mtc{mtc}_np_{non_parametric}_fp_{force_parametric}_{modifier}.pdf')
    plt.close()
