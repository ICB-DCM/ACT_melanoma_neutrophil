#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
# '''=================================================
# @Project: ACT_melanoma_neutrophil (data: act_neutrophil_data_repository)
# @Author : Gemma van der Voort
# @Desc   :
# @Desc updated: Sum early-gate counts across panels and export normalized totals.
# calculate the actual all_cell counts based on the beads and add them over 2 or 3 panels.
# sum up the amounts for the three panels per organ.
#  correct for the 86 and 87 missing pmel beads, discard pmel panel info and correct to total LN for other two panels.
# plots a histogram to check processing in plausible, but leaves paper-style figure for vis script.
# '''=================================================
import os
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

PROJECT_ROOT = Path(__file__).resolve()
while PROJECT_ROOT.name != "ACT_melanoma_neutrophil" and PROJECT_ROOT.parent != PROJECT_ROOT:
    PROJECT_ROOT = PROJECT_ROOT.parent
if PROJECT_ROOT.name != "ACT_melanoma_neutrophil":
    raise RuntimeError("ACT_melanoma_neutrophil root not found")
sys.path.append(str(PROJECT_ROOT / "helper_files"))

from paths_parameters import data_repo_path

plt.style.use(str(PROJECT_ROOT / "helper_files" / "journal_style.mplstyle"))
processed_data_path = f'{data_repo_path}/processed/flow_cytometry'
fig_path = f'{data_repo_path}/figures/flow_cytometry/preprocessing'
os.makedirs(fig_path, exist_ok=True)
#%%
conditions_flow = ['Naive', 'Tumour-bearing 3-5 mm', 'Cyclo only', 'ACT-d3', 'ACT-d7',
       'ACT-d14', 'ACT_Relapse', 'ACT-d14_no CpG', 'ACT-d7_no CpG']
panel_amount_flow = [2,2,2,3,3,3,3,3,3]  # control experiments don't yet have pmel panel. fractions for LNs.
dict_cond_panel = {cond: panel for cond, panel in zip(conditions_flow, panel_amount_flow)}
lymph_nodes = ['inLNr', 'inLNl', 'brLNr']
#%% read in the data
df_early_gates = pd.read_csv(f'{processed_data_path}/early_gates_slo_counts.csv', index_col=0)
# drop all not in lymph_nodes
df_early_gates_ln = df_early_gates[df_early_gates['Organ'].isin(lymph_nodes)].copy()
# %%add column bead_factor, 5000 (original amount of beads) / bead_count and organ_factor_LN
df_early_gates_ln['bead_factor'] = 5000 / df_early_gates_ln['bead_count']
df_early_gates_ln['organ_factor'] = 1/ df_early_gates_ln['Organ_fraction_LN']
df_early_gates_ln['count_multiplication_factor'] = df_early_gates_ln['bead_factor'] * df_early_gates_ln['organ_factor']
# %% now loop over count columns (except bead) and create a new column with the normalised counts
count_columns = ['all_count', 'lym_count', 'single_count', 'live_count', 'cd45_count']
for col in count_columns:
    df_early_gates_ln[col + '_norm'] = df_early_gates_ln[col] * df_early_gates_ln['count_multiplication_factor']
#%% now all cells per lymph node. as a check, count the occurrences of each lymph node and mouse_id pair.
# I expect 2 or 3 per combination.
panel_summary = (
    df_early_gates_ln
    .groupby(['Mouse_ID', 'Organ', 'Condition'])['Panel']
    .unique()
    .reset_index()
)
panel_summary = panel_summary.sort_values(['Condition', 'Mouse_ID', 'Organ'])
# %%inspected manually. There are some missing panels/organs. But these are consistent with the processed data
# Maike provided, except for BAL-3605. Remove that mouse from analysis. see flow_codex_panel_successes.ods.
# reason: gating incomplete, probably experimental issue. No note found.
df_early_gates_ln = df_early_gates_ln[~df_early_gates_ln['Mouse_ID'].isin(['BAL-3605'])]
panel_summary = panel_summary[~panel_summary['Mouse_ID'].isin(['BAL-3605'])]
#%% now calculate total live cell counts per organ and mouse_id, sum of each of the panels. For the samples with issues,
# correct.
df_total_counts = panel_summary.copy()
for col in count_columns:
    # sum the counts over the panels
    df_total_counts[col + '_norm'] = 0
norm_columns = [col + '_norm' for col in count_columns]
df_total_counts[norm_columns] = df_total_counts[norm_columns].astype(float)

for index, row in panel_summary.iterrows():
    condition = row['Condition']
    correct_amount_of_panels = dict_cond_panel[condition]
    actual_amount_of_panels = len(row['Panel'])

    print(condition, correct_amount_of_panels, actual_amount_of_panels)

    match_filter = (
        (df_early_gates_ln['Mouse_ID'] == row['Mouse_ID']) &
        (df_early_gates_ln['Organ'] == row['Organ']) &
        (df_early_gates_ln['Condition'] == condition)
    )
    summed_counts = df_early_gates_ln.loc[match_filter, norm_columns].sum()

    if actual_amount_of_panels == correct_amount_of_panels:
        factor = 1
    elif actual_amount_of_panels == 2 and correct_amount_of_panels == 3:
        factor = 1.5
    elif actual_amount_of_panels == 1 and correct_amount_of_panels == 3:
        factor = 3
    elif actual_amount_of_panels == 1 and correct_amount_of_panels == 2:
        factor = 2
    else:
        print('Error: amount of missing panels not taken into account')
        continue

    df_total_counts.loc[index, norm_columns] = summed_counts * factor
# add 'dead_cells' column
df_total_counts['dead_fraction_of_single'] = df_total_counts['all_count_norm'] / df_total_counts['single_count_norm']
norm_columns += ['dead_fraction_of_single']
#%% quick plot of data distributions per condition and organ to check for plausibility and outliers.
# seaborn histplot of df_total_counts, hue by condition, facet by organ
hues = ['Condition', 'Organ']
for hue in hues:
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    axes = axes.flatten()

    for i, col in enumerate(norm_columns):
        sns.histplot(data=df_total_counts, x=col, hue=hue, ax=axes[i], bins=20, legend=True)
        axes[i].set_title(col)
        axes[i].set_xlabel('Cell counts')
        axes[i].set_ylabel('Frequency')
    sns.despine()
    plt.tight_layout()
    plt.savefig(f'{fig_path}/total_counts_distribution_{hue}.pdf', bbox_inches='tight')
    plt.close()
#%% revision: highlight dead cell fraction in a histogram, per condition.
fig, axes = plt.subplots(1, 1, figsize=(10, 6))
sns.histplot(data=df_total_counts, x='dead_fraction_of_single', hue='Condition', ax=axes, bins=100, legend=True,
             kde=True)
# axes.set_title('Dead cell fraction of single cells')
axes.set_xlabel('Dead cells [% of single cells]')
# set x limit to 2
axes.set_xlim(0, 2)
axes.set_ylabel('Frequency')
sns.despine()
plt.tight_layout()
plt.savefig(f'{fig_path}/dead_fraction_of_single_histogram_by_condition_outlier_rm.pdf', bbox_inches='tight')
#%% save the data
df_total_counts.to_csv(f'{processed_data_path}/early_gates_slo_counts_summed_over_panels.csv', index=False)
