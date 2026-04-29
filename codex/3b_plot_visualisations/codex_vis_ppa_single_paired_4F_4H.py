#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
# '''=================================================
# @Project: ACT_melanoma_neutrophil (data: act_neutrophil_data_repository)
# @Author : Gemma van der Voort
# @Desc   : do single paired plots for ppa data, d14, clLN vs tdLN, with mtc for multiple explored threhsolds
# revision: Updated with max_iterations to 25
# '''=================================================
import sys
from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt

import os
from statsmodels.stats.multitest import multipletests


PROJECT_ROOT = Path(__file__).resolve()
while PROJECT_ROOT.name != "ACT_melanoma_neutrophil" and PROJECT_ROOT.parent != PROJECT_ROOT:
    PROJECT_ROOT = PROJECT_ROOT.parent
if PROJECT_ROOT.name != "ACT_melanoma_neutrophil":
    raise RuntimeError("ACT_melanoma_neutrophil root not found")
sys.path.append(str(PROJECT_ROOT / "helper_files"))

from paths_parameters import x_axis_tick_labels_flow, bar_colours_hex_blood_flow, lymph_nodes_label, marker_per_organ_dict, data_repo_path
from functions_stat_test import ratio_paired_test
from functions_vis import boxplot_paired, single_paired_plot

plt.style.use(str(PROJECT_ROOT / "helper_files" / "journal_style_old.mplstyle"))
processed_spatial_data_dir = f'{data_repo_path}/processed/codex'
stats_dir = f'{processed_spatial_data_dir}/statistics'
os.makedirs(stats_dir, exist_ok=True)
modifier_base = f'20241021_all_ln_tree6_rep2'
version_key = 'v4_8'
img_modifier = 'v4_corrected'
force_parametric = False # set after initial normalicy testing. Since this script is used for single paired plots only,
# I don't need to enforce all to a parametric test after initial testing.
metacluster_key = 'Metacluster v4_8_filtered_0.5_maxit_25_v4_8'
modifier = f'{img_modifier}_{metacluster_key}_fp_{force_parametric}'
df_t_zone_proximity = pd.read_csv(f"{processed_spatial_data_dir}/df_t_zone_proximity_neutrophils_{img_modifier}_{metacluster_key}.csv")
fig_dir = f'{data_repo_path}/figures/codex/manuscipt_plots/patch_proximity_analysis'
os.makedirs(fig_dir, exist_ok=True)

#%%
condition_col = 'Treatment'
condition = 'ACT_d14' # only interested in d14
mouse_id_col = 'mouse_id'
value_cols = ['abundance_below_0_px', 'abundance_above_0_px',
              'abundance_below_60_px', 'abundance_above_60_px',
              'abundance_below_180_px', 'abundance_above_180_px',
              'abundance_below_300_px', 'abundance_above_300_px']
value_cols_y_label = ['Neutrophil fraction\ninside T cell zone',
                      'Neutrophil fraction\noutside T cell zone',
                      'Neutrophil fraction\nwithin 20 micron of T cell zone',
                      'Neutrophil fraction\nabove 20 micron outside T cell zone',
                      'Neutrophil fraction\nwithin 60 micron of T cell zone',
                      'Neutrophil fraction\nabove 60 micron outside T cell zone',
                      'Neutrophil fraction\nwithin 100 micron of T cell zone',
                      'Neutrophil fraction\nabove 100 micron outside T cell zone']
#%% df_t_zone_proximity has some 0's and 1's. add a small epsilon offset so ratio tests can be calculated OK.
epsilon = 1e-6
df_t_zone_proximity[value_cols] = df_t_zone_proximity[value_cols].clip(lower=epsilon, upper=1 - epsilon)
df_t_zone_proximity[value_cols] = df_t_zone_proximity[value_cols]*100 # so I can use % instead of fractions
group_col = 'Organ'
ln_order = ['inLNr', 'inLNl']
LN_to_compare = [('inLNr', 'inLNl')]
    # df with the results of the statistical tests.
stats_df = pd.DataFrame(columns=['Cell type', condition_col, 'LN 1', 'LN 2', 'p-value', 'test_used'])
for value_col in value_cols:
    for LN1, LN2 in LN_to_compare:
        # for some mice not both lymph nodes were processed. drop these mice
        ln_df_LN1 = df_t_zone_proximity[(df_t_zone_proximity[group_col] == LN1) & (df_t_zone_proximity[condition_col] == condition)]
        ln_df_LN2 = df_t_zone_proximity[(df_t_zone_proximity[group_col] == LN2) & (df_t_zone_proximity[condition_col] == condition)]
        # check if the mouse_ids are the same and in the same order
        if not  ln_df_LN1[mouse_id_col].tolist() == ln_df_LN2[mouse_id_col].tolist():
            # drop the mouse ID that is not in both conditions
            only_ln1 = set(ln_df_LN1[mouse_id_col].tolist()) - set(ln_df_LN2[mouse_id_col].tolist())
            only_ln2 = set(ln_df_LN2[mouse_id_col].tolist()) - set(ln_df_LN1[mouse_id_col].tolist())
            print(f'Only in {LN1}: {only_ln1}, only in {LN2}: {only_ln2}. Dropped single mice')
            ln_df_LN1 = ln_df_LN1[~ln_df_LN1[mouse_id_col].isin(only_ln1)]
            ln_df_LN2 = ln_df_LN2[~ln_df_LN2[mouse_id_col].isin(only_ln2)]
        else:
            # print(f'Mice in {LN1} and {LN2} are the same')
            pass
        t_stat, p_val, test_used = ratio_paired_test(ln_df_LN1, ln_df_LN2, value_col, non_parametric=False,
                                                     force_parametric=force_parametric)
        stats_df = pd.concat([stats_df, pd.DataFrame([[value_col, condition, LN1, LN2, p_val,
                                                       test_used]],
                                                     columns=['Cell type', condition_col, 'LN 1', 'LN 2', 'p-value',
                                                              'test_used'])])
# write to file
stats_df.to_csv(f'{stats_dir}/ln_ppa_codex_ratio_paired_{modifier}.csv')

#%% do mtc on exploratory ring distances only (60, 180, 300 px);
# 0 px (within region) is the pre-specified primary endpoint and not part of the exploratory family
# subset stats_df for the 3 exploratory distances:
exploratory_distances = ['abundance_above_60_px', 'abundance_above_180_px', 'abundance_above_300_px',
                         'abundance_below_60_px', 'abundance_below_180_px', 'abundance_below_300_px']
stats_df_exploratory = stats_df[stats_df['Cell type'].isin(exploratory_distances)].copy()
# apply holm
stats_df_exploratory['p_adj'] = multipletests(stats_df_exploratory['p-value'], method='holm')[1]
# save to file
stats_df_exploratory.to_csv(f'{stats_dir}/ln_ppa_codex_ratio_paired_exploratory_mtc_{modifier}.csv', index=False)
#%%
condition = 'ACT_d14'
value_col = 'abundance_above_60_px'
value_col_y_label = 'Neutrophil \n[% >20µm outside T cell zone]'
colour = "#fa815fff" # d14 orange
markers = ['o', 's']
x_pair = ('inLNr', 'inLNl')

fig, ax = plt.subplots(1, 1, figsize=(3, 6))
single_paired_plot(df_t_zone_proximity, x_col=condition_col, y_col=value_col, id_col= mouse_id_col,
                   condition=condition, organ_col=group_col,
                   x_pair= x_pair, stats_df=stats_df_exploratory, p_val = 'p_adj',
                   markers = markers, scatter_kwargs={'color': colour, 's': 50}, ax=ax,
                   remove_unpaired=False)
ax.set_ylabel(value_col_y_label, fontsize=22)
ax.set_xlabel('')
# remove 105 tick
ax.set_yticks([70, 80, 90, 100])
plt.savefig(f'{fig_dir}/Neutrophil_ppa_{condition}_{value_col}_ratio_paired_{modifier}.pdf')
#%%
condition = 'ACT_d14'
value_col = 'abundance_below_0_px'
value_col_y_label = 'Neutrophil\n % inside T cell zone'  # of total neutrophils
colour = "#fa815fff" # d14 orange
markers = ['o', 's']
x_pair = ('inLNr', 'inLNl')

fig, ax = plt.subplots(1, 1, figsize=(3, 6))
single_paired_plot(df_t_zone_proximity, x_col=condition_col, y_col=value_col, id_col= mouse_id_col,
                   condition=condition, organ_col=group_col,
                   x_pair= x_pair, stats_df=stats_df, markers = markers, scatter_kwargs={'color': colour, 's': 50}, ax=ax,
                   remove_unpaired=False)
ax.set_ylabel(value_col_y_label)
ax.set_xlabel('')
plt.savefig(f'{fig_dir}/Neutrophil_ppa_{condition}_ratio_paired_{modifier}.pdf')
