#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
# '''=================================================
# @Project: ACT_melanoma_neutrophil (data: act_neutrophil_data_repository)
# @Author : Gemma van der Voort
# @Desc   : paired plot but without the brLNr for less crowding.
# @Desc updated: Main-condition paired pan-immune plots with single-paired highlights.
# update 8.8: us english, drop useless x axis titles, add single paired d7 plot
# update 2026.01.14: narrower plots for act-only cells, rm non-significant lines.
# revision: consider sex as a variable. check what adding them does for the visualisations.
# '''=================================================
import os
import sys
from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
from matplotlib.transforms import ScaledTranslation
import numpy as np

PROJECT_ROOT = Path(__file__).resolve()
while PROJECT_ROOT.name != "ACT_melanoma_neutrophil" and PROJECT_ROOT.parent != PROJECT_ROOT:
    PROJECT_ROOT = PROJECT_ROOT.parent
if PROJECT_ROOT.name != "ACT_melanoma_neutrophil":
    raise RuntimeError("ACT_melanoma_neutrophil root not found")
sys.path.append(str(PROJECT_ROOT / "helper_files"))

from paths_parameters import x_axis_tick_labels_flow, bar_colours_hex_blood_flow, lymph_nodes_label, \
    marker_per_organ_dict, neutrophil_label, data_repo_path
from functions_stat_test import ratio_paired_test
from functions_vis import boxplot_paired_no_brLNr, single_paired_plot

plt.style.use(str(PROJECT_ROOT / "helper_files" / "journal_style.mplstyle"))

# paths etc
processed_data_path = f'{data_repo_path}/processed/flow_cytometry'
stats_dir = f'{processed_data_path}/statistics'
os.makedirs(stats_dir, exist_ok=True)
fig_path = f'{data_repo_path}/figures/flow_cytometry/main_conditions'
ln_df_raw = pd.read_csv(f'{processed_data_path}/ln_pan_immune_normalised_v2_cd45_corrected.csv')
modifier = 'cd45_corrected_narrow'
# make dir and sub-dirs if they don't exist
os.makedirs(fig_path, exist_ok=True)
os.makedirs(f'{fig_path}/multiple_paired', exist_ok=True)
os.makedirs(f'{fig_path}/single_paired', exist_ok=True)

females_rm = False
#%%
# drop condition if it contains the substring 'no CpG'
ln_df = ln_df_raw[~ln_df_raw['Condition'].str.contains('no CpG')]
df_condition_names = ln_df['Condition'].unique()
#%% add sex info
naive_female_mice = ['BAL-3418', 'BAL-3419', 'BAL-3531', 'BAL-3428', 'BAL-3429', 'BAL-3431', 'BAL-3432', 'BAL-3433',
                     'BAL-4023', 'BAL-4024', 'BAL-4025', 'BAL-4026']
tb_female_mice = ['BAL-3727', 'BAL-3729', 'BAL-3728', 'BAL-3656', 'BAL-3725', 'BAL-3726', 'HOE-2235', 'HOE-2237']
cy_female_mice = ['BAL-3611', 'BAL-3614', 'BAL-3612', 'BAL-3615', 'HOE-2251']
untreated_female_mice = ['BAL-3829', 'BAL-3827', 'BAL-3828', 'HOE-2234']
d3_female_mice = ['HOE-2273', 'HOE-2265', 'HOE-2261', 'HOE-2264', 'HOE-2263']
d7_female_mice = ['HOE-2390', 'HOE-2391', 'HOE-2388', 'HOE-2389', 'HOE-2387']
d14_female_mice = ['HOE-2384', 'HOE-2385', 'HOE-2386', 'HOE-2393']
relapse_female_mice = ['HOE-2285', 'HOE-2296', 'HOE-2289', 'HOE-2294', 'HOE-2292', 'HOE-2293', 'HOE-2288', 'HOE-2287']
# female_mice_ids = d7_female_mice + d14_female_mice # IDs are unique, so just generate a long list.
female_mice_ids = (naive_female_mice+tb_female_mice+cy_female_mice+untreated_female_mice+d3_female_mice+
                   d7_female_mice+d14_female_mice+relapse_female_mice)
if females_rm:
    ln_df = ln_df[~ln_df['Mouse_ID'].isin(female_mice_ids)].copy()
    modifier += '_fem_removed'
value_cols = ['Pmel T percent',
              'Neutrophils percent',]
value_cols_y_label = ['Pmel-1 T cells [% of CD45.2$^+$ cells]',
                        'Neutrophils [% of CD45.2$^+$ cells]',
                      ]
modifier_2 = 'percent'
# alternative: run with counts (remember to rename downstream)
# value_cols = ['Pmel T  | Count normalised',
#                 'Neutrophils | Count normalised']
# value_cols_y_label = ['Pmel-1 T cells [normalised count]',
#                         'Neutrophils [normalised count]', ]
# modifier_2 = 'counts'
group_col = 'LN_type'
ln_order = ['inLNr', 'inLNl']
LN_to_compare = [('inLNr', 'inLNl')]  #  ('inLNl', 'brLNr'), ('inLNr', 'brLNr')
# df with the results of the statistical tests.
stats_df = pd.DataFrame(columns=['Cell type', 'Condition', 'LN 1', 'LN 2', 'p-value', 'test_used'])
for value_col in value_cols:
    for condition in df_condition_names:
        for LN1, LN2 in LN_to_compare:
            # for some mice not both lymph nodes were processed. drop these mice
            ln_df_LN1 = ln_df[(ln_df[group_col] == LN1) & (ln_df['Condition'] == condition)]
            ln_df_LN2 = ln_df[(ln_df[group_col] == LN2) & (ln_df['Condition'] == condition)]
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
stats_df.to_csv(f'{stats_dir}/ln_pan_immune_stats_no_CpG_ratio_paired_no_brLNr_{modifier}_{modifier_2}.csv')

#%% modify: no brLnr
n_categories = len(df_condition_names)
ncols = 1
fig_width = ((n_categories * 1.5) - 1) * ncols  # normal, broader than single due to 2 conditions.

markers_LN = [marker_per_organ_dict[ln] for ln in ln_order]

for cell_type_i, cell_type in enumerate(value_cols):  # make single plots for each cell type
    print(f'Working on {cell_type}')
    # adapt for pmel t cells, who are only in the later conditions
    pmel_conditions = ['ACT-d3', 'ACT-d7', 'ACT-d14', 'ACT_Relapse']
    if cell_type == 'Pmel T percent' or cell_type == 'Pmel T  | Count normalised':
        ln_df_pmel = ln_df[ln_df['Condition'].isin(pmel_conditions)]
    for LN1, LN2 in LN_to_compare:  # plot for each LN combination, too crowed otherwise.
        try:
            data = ln_df_pmel if cell_type == 'Pmel T percent' or cell_type == 'Pmel T  | Count normalised' else ln_df  # [cell_type]
            # remove the brLNr from data
            data = data[data['LN_type'] != 'brLNr']
            # remove the brLNr from the categories in LN_type
            ln_order_no_brLNr = [ln for ln in ln_order if ln != 'brLNr']


            # get stats df for cell type and LN combination
            stats_df_cell = stats_df[(stats_df['Cell type'] == cell_type) & (stats_df['LN 1'] == LN1) & (stats_df['LN 2']
                                                                                                         == LN2)]
            if cell_type == 'Pmel T percent' or cell_type == 'Pmel T  | Count normalised':
                fig, ax = plt.subplots(1, ncols, figsize=(5, 6)) # ## updated, let's see.
                stats_df_cell = stats_df_cell[stats_df_cell['Condition'].isin(pmel_conditions)]
                bar_colours_hex_blood_flow_pmel = bar_colours_hex_blood_flow[3:]  # remove the first 3 colours
                x_axis_tick_labels_blood_flow_pmel = x_axis_tick_labels_flow[3:]
                stats_df_cell = stats_df_cell.drop(['Cell type', 'Condition'],
                                                   axis=1)  # drop the LN_type column for the plot
                # troubleshoot female mice colouring.
                boxplot_paired_no_brLNr(data, 'Condition', cell_type, x_group_col='LN_type', stats_df=stats_df_cell, ax=ax,
                               x_group_order=ln_order, markers=markers_LN, condition_colours=bar_colours_hex_blood_flow_pmel,
                                        omit_ns=True, mouse_id_female = female_mice_ids, mouse_id_col='Mouse_ID'
                               )
                ax.set_ylabel(f'{value_cols_y_label[cell_type_i]}', fontsize=20)
                # ax.set_xlabel('Condition')
                ax.set_xlabel('')
                ax.set_xticklabels(x_axis_tick_labels_blood_flow_pmel, fontsize=18)
                # if counts in cell type, format y axis ticks as 10^x
                if 'Count' in cell_type:
                    ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
                    ax.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
            else:
                fig, ax = plt.subplots(1, ncols, figsize=(9, 6))  # ## updated, let's see.
                stats_df_cell = stats_df_cell.drop(['Cell type', 'Condition'],
                                                   axis=1)  # drop the LN_type column for the plot
                boxplot_paired_no_brLNr(data, 'Condition', cell_type, x_group_col='LN_type', stats_df=stats_df_cell, ax=ax,
                               x_group_order=ln_order, markers=markers_LN, condition_colours=bar_colours_hex_blood_flow,
                                        omit_ns=True, mouse_id_female=female_mice_ids, mouse_id_col='Mouse_ID'
                              )
                ax.set_ylabel(f'{value_cols_y_label[cell_type_i]}')
                # ax.set_xlabel('Condition')
                ax.set_xlabel('')
                ax.set_xticklabels(x_axis_tick_labels_flow)
                if 'Count' in cell_type:
                    ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
                    ax.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
                    fig.canvas.draw() # draw the canvas to fill the text object
                    text = ax.yaxis.get_offset_text()
                    dx, dy = 0.0, 10 / 72  # 10 points upward; 72 points per inch
                    offset = ScaledTranslation(dx, dy, fig.dpi_scale_trans)
                    text.set_transform(text.get_transform() + offset)
            # if cell_type contains |, replace it with _ for saving
            cell_type_save = cell_type.replace(' | ', '_')
            plt.savefig(f'{fig_path}/multiple_paired/{cell_type_save}_ratio_paired_boxplot_'
                        f'no_brLNr_{modifier}_{modifier_2}_no_title_sex_black.pdf')
            plt.close()
        except TypeError:
            print(f'No data for {cell_type} in {LN1} and {LN2}')
            pass
#%% single paired plots  (as example/highlight)
neutrophil_label_single = 'Neutrophils % of CD45.2$^+$ cells'
condition = 'ACT-d14'
fig, ax = plt.subplots(1, 1, figsize=(3, 6))
colour = "#fa815fff"
markers = ['o', 's']
x_pair = ('inLNr', 'inLNl')
single_paired_plot(ln_df, 'Condition', 'Neutrophils percent', 'Mouse_ID', condition, 'LN_type',
                   x_pair, stats_df, markers = markers, scatter_kwargs={'color': colour, 's': 50}, ax=ax,
                   remove_unpaired=False, mouse_id_female=female_mice_ids)
ax.set_ylabel(neutrophil_label_single, size=20)
# increase x tick label size
ax.tick_params(axis='x', labelsize=20)
ax.set_xlabel('')
# plt.show()
plt.savefig(f'{fig_path}/single_paired/Neutrophils_{condition}_ratio_paired_{modifier}_{modifier_2}_sex_black.pdf')
# %% same for day 7
condition = 'ACT-d7'
fig, ax = plt.subplots(1, 1, figsize=(3, 6))
colour = "#fde725ff"
markers = ['o', 's']
x_pair = ('inLNr', 'inLNl')
single_paired_plot(ln_df, 'Condition', 'Neutrophils percent', 'Mouse_ID', condition, 'LN_type',
                   x_pair, stats_df, markers = markers, scatter_kwargs={'color': colour, 's': 50}, ax=ax,
                   remove_unpaired=False, mouse_id_female=female_mice_ids)
ax.set_ylabel(neutrophil_label_single, size=20)
ax.set_xlabel('')
ax.tick_params(axis='x', labelsize=20)
# plt.show()
plt.savefig(f'{fig_path}/single_paired/Neutrophils_{condition}_ratio_paired_{modifier}_{modifier_2}_sex_black.pdf')
