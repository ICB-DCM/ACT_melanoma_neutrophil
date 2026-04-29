#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
# '''=================================================
# @Project: ACT_melanoma_neutrophil (data: act_neutrophil_data_repository)
# @Author : Gemma van der Voort
# @Desc   : import all relevant datasets and use them to construct a dataset overview table to accompany the manuscript.
# @Desc updated: Build dataset overview tables from processed flow, tumor growth, and CODEx data.
# note: run this script after preprocessing the tumor growth, codex and flow cytometry modalities, since it uses
# preprocessed data from these.
# '''=================================================
import sys
from pathlib import Path

import pandas as pd
import os
import numpy as np

PROJECT_ROOT = Path(__file__).resolve()
while PROJECT_ROOT.name != "ACT_melanoma_neutrophil" and PROJECT_ROOT.parent != PROJECT_ROOT:
    PROJECT_ROOT = PROJECT_ROOT.parent
if PROJECT_ROOT.name != "ACT_melanoma_neutrophil":
    raise RuntimeError("ACT_melanoma_neutrophil root not found")
sys.path.append(str(PROJECT_ROOT / "helper_files"))

from paths_parameters import data_repo_path

processed_data_path = f'{data_repo_path}/overview_of_experiments'
#%% load, setting the 2nd row as header
# NOTE: v5 has 1 mice treatment cohort updated from unt to TB (already did in data, but label was still old. Fixed.
df_data_overview = pd.read_excel(f'{processed_data_path}/Data_consistency_check_jan2026_v6.ods', engine='odf', header=[1,2])
# some headers have extra info in the first row. where this row is not nan, append the contents to the header
def merge_double_header(columns, sep=" | "):
    new_cols = []

    for lvl0, lvl1 in columns:
        # normalize missing values
        lvl0 = None if pd.isna(lvl0) or str(lvl0).startswith("Unnamed") else str(lvl0).strip()
        lvl1 = None if pd.isna(lvl1) or str(lvl1).startswith("Unnamed") else str(lvl1).strip()

        if lvl0 and lvl1:
            new_cols.append(f"{lvl0}{sep}{lvl1}")
        elif lvl1:
            new_cols.append(lvl1)
        elif lvl0:
            new_cols.append(lvl0)
        else:
            new_cols.append(None)

    return new_cols
df_data_overview.columns = merge_double_header(df_data_overview.columns)

# %%list unique headers in 'Treatment cohort' column
treatment_cohorts = df_data_overview['Treatment cohort'].unique().tolist()
#%% correct cohort spelling
df_data_overview['Treatment cohort'] = df_data_overview['Treatment cohort'].str.replace('Naive', 'Naïve')
df_data_overview['Treatment cohort'] = df_data_overview['Treatment cohort'].str.replace('Tumour-bearing (3-5mm)', 'Tumour-bearing (3-5 mm)')
df_data_overview['Treatment cohort'] = df_data_overview['Treatment cohort'].str.replace('Tumour-bearing (8-10mm)', 'Tumour-bearing (8-10 mm)')
relevant_cohorts = ['Naïve', 'Tumour-bearing (3-5 mm)',  # in data overview file
                    'Tumour-bearing (3-5 mm) + 1 day cyclophosphamide', 'ACT_d3',
                    'ACT_d7', 'ACT_d14', 'ACT_relapse (long-term)',
                    'Tumour-bearing (8-10 mm)',
                    'ACT_d7_no CpG', 'ACT_d14_no CpG']
                    # rm 'Tumour-bearing (3-5 mm) + 4 day cyclophosphamide'
cohort_names =['Naïve', 'Tumor bearing', # these will be used in table, match others to this format
                'Cy', 'ACT day 3',
                'ACT day 7', 'ACT day 14', 'Relapse',
               'Untreated', # ='Tumor bearing 8-10 mm',
               'ACT day 7, no CpG', 'ACT day 14, no CpG']

# create df for table 1. rows: all mice per cohort, mice with tumor growth, mice used for flow, mice used for codex.
# columns: relevant cohorts
table1_df = pd.DataFrame(index=['Total mice', 'Mice with tumor growth', 'Mice used for flow', 'Mice used for codex'],
                         columns=cohort_names+ ['Total'])
# %%tumor growth not in table, but only from tumor growth data. Import from the scripts on this.
tg_df = pd.read_csv(f'{data_repo_path}/processed/tumor_growth/per_treatment/'
                    f'tumor_growth_mice_no_growth_summary_incl_no_CpGincl_tb3_tb8_newstyle_data_corrections_sep2025.csv')
# rename treatment name 'Tumor bearing 8-10mm to Untreated
tg_df['Treatment'] = tg_df['Treatment'].str.replace('Tumor bearing 8-10 mm', 'Untreated')
# reorder treatment column to match cohort names
tg_df['Treatment'] = pd.Categorical(tg_df['Treatment'], categories=cohort_names)
tg_df = tg_df.sort_values(by=['Treatment'])

#%% fill df table 1
df_data_overview['FLOW'].unique()
flow_present = [3, 'inR', 'inR brR', 'inL, brR', 'V']
df_data_overview['CODEX'].unique()
codex_present = [3, '2 (inR, inL)', '3 (CODEX ran twice 100.2 and 100.8)']  # exclude the frozen but not processed one
for cohort_index, cohort in enumerate(cohort_names):
    cohort_name_overview = relevant_cohorts[cohort_index]
    if cohort== 'Naïve':
        total_mice = df_data_overview[df_data_overview['Treatment cohort'] == 'Naïve'].shape[0]
        mice_with_tumor_growth = 0
    else:
        total_mice = tg_df[tg_df['Treatment'] == cohort]['All mice'].values[0]
        mice_without_tumor_growth = tg_df[tg_df['Treatment'] == cohort]['Mice without tumor growth'].values[0]
        mice_with_tumor_growth = total_mice - mice_without_tumor_growth
    flow_mice = df_data_overview[(df_data_overview['Treatment cohort'] == cohort_name_overview) &
                                 (df_data_overview['FLOW'].isin(flow_present))].shape[0]
    codex_mice = df_data_overview[(df_data_overview['Treatment cohort'] == cohort_name_overview) &
                                 (df_data_overview['CODEX'].isin(codex_present))].shape[0]
    table1_df.at['Total mice', cohort] = total_mice
    table1_df.at['Mice with tumor growth', cohort] = mice_with_tumor_growth
    table1_df.at['Mice used for flow', cohort] = flow_mice
    table1_df.at['Mice used for codex', cohort] = codex_mice
# add total column
table1_df['Total'] = table1_df.sum(axis=1)

#%% # check that the total mice per cohort in df_data_overview matches the total mice counted in tgc.
# for cohort in relevant_cohorts:
#     total_mice = df_data_overview[df_data_overview['Treatment cohort'] == cohort].shape[0]
#     print(f'Working on cohort {cohort} with total mice {total_mice}')

#%% check flow amounts by loading the pan-immune flow data and checking mice-IDs per cohort
processed_data_path_flow_pi = f'{data_repo_path}/processed/flow_cytometry'
df_pi = pd.read_csv(f'{processed_data_path_flow_pi}/ln_pan_immune_normalised_v2_cd45_corrected.csv')
# get unique conditions
# df_pi['Condition'].unique()
for condition in df_pi['Condition'].unique():
    mice_in_condition = df_pi[df_pi['Condition'] == condition]['Mouse_ID'].unique().tolist()
    n_mice_in_condition = len(mice_in_condition)
    # print(f'Condition {condition} has {n_mice_in_condition} mice: {mice_in_condition}')
# all correct
#%% repeat for codex
processed_codex_data_dir = f'{data_repo_path}/processed/codex'
modifier_base = f'20241021_all_ln_tree6_rep2'
modifier_finetuning = 'layer_uns_rem_finetuned_DCs_fixed'
modifier_postprocessing = f'transcleaned_postprocessed'
file_name_codex = f'{modifier_base}_{modifier_finetuning}_{modifier_postprocessing}'
df_codex = pd.read_csv(f'{processed_codex_data_dir}/df_grouped_coi_{file_name_codex}_all_cells.csv', index_col=0)
for treatment in df_codex['Treatment'].unique():
    mice_in_treatment = df_codex[df_codex['Treatment'] == treatment]['mouse_id'].unique().tolist()
    n_mice_in_treatment = len(mice_in_treatment)
    # print(f'Treatment {treatment} has {n_mice_in_treatment} mice: {mice_in_treatment}')
# save and export table 1
table1_df.to_csv(f'{processed_data_path}/dataset_overview_table_1_jan2026.csv')
#%% transpose and save as csv (since this one doesn't fit well)
table1_df.T.to_csv(f'{processed_data_path}/dataset_overview_table_1_transposed_jan2026.csv')

############################### table 1, version 2: add average cell counts ############################
# load total counts and count the mean (st dev) per lymph node per condition
processed_data_path_early_gates = f'{data_repo_path}/processed/flow_cytometry'
df_total_counts = pd.read_csv(f'{processed_data_path_early_gates}/early_gates_slo_counts_summed_over_panels.csv')
# df_total_counts['Condition'].unique()
flow_to_table_mapping = {
    'Naive': 'Naïve',
    'Tumour-bearing 3-5 mm': 'Tumor bearing',
    'Cyclo only': 'Cy',
    'ACT-d3': 'ACT day 3',
    'ACT-d7': 'ACT day 7',
    'ACT-d14': 'ACT day 14',
    'ACT_Relapse': 'Relapse',
    'ACT-d7_no CpG': 'ACT day 7, no CpG',
    'ACT-d14_no CpG': 'ACT day 14, no CpG'
}
# %%add row 'Average total cell count per LN (x10^6) ± stdev' to table1_df
avg_counts_row = 'Flow: Mean live cell count per LN [$10^6$] (±sd)'
table1_df.loc[avg_counts_row] = ''
for cohort_index, cohort in enumerate(cohort_names):
    if cohort == 'Untreated':
        continue  # no data for this condition in flow
    # get df for this condition
    condition_name_flow = [key for key, value in flow_to_table_mapping.items() if value == cohort][0]
    print(cohort, condition_name_flow)
    df_condition = df_total_counts[df_total_counts['Condition'] == condition_name_flow]
    # calculate mean and stdev over all lymph nodes
    mean_count = df_condition['live_count_norm'].mean()
    std_count = df_condition['live_count_norm'].std()
    # format as string with 10^6 (currently in units of cells) with 2 decimal places
    mean_std_str = f'{mean_count/1e6:.2f} (±{std_count/1e6:.2f})'
    table1_df.at[avg_counts_row, cohort] = mean_std_str
# calculate total column for this row
mean_total = df_total_counts['all_count_norm'].mean()
std_total = df_total_counts['all_count_norm'].std()
mean_std_total_str = f'{mean_total/1e6:.2f} (±{std_total/1e6:.2f})'
table1_df.at[avg_counts_row, 'Total'] = mean_std_total_str
# %% repeat for codex
# df_codex['Treatment'].unique()
codex_to_table_mapping = {
    'tumour_bearing': 'Tumor bearing',
    'Cyclo_d1': 'Cy',
    'ACT_d3': 'ACT day 3',
    'ACT_d7': 'ACT day 7',
    'ACT_d14': 'ACT day 14'
}
codex_counts_row = 'CODEX: Mean live cell count per LN slice [$10^3$] (±sd)'
table1_df.loc[codex_counts_row] = ''
# codex df needs some preprocessing: rows are per cell type. group df on treatment, mouse_id, ln_type to get total cells per mouse per ln
df_codex_totals = df_codex.groupby(['mouse_id', 'Organ', 'Treatment']).agg(sum).reset_index()
# drop cell type column and population column
df_codex_totals = df_codex_totals.drop(columns=['Cell type', 'Population'])
#%%
for cohort_index, cohort in enumerate(cohort_names):
    if cohort not in codex_to_table_mapping.values():
        continue  # no data for this condition in codex
    # get df for this condition
    condition_name_codex = [key for key, value in codex_to_table_mapping.items() if value == cohort][0]
    print(cohort, condition_name_codex)
    df_condition = df_codex_totals[df_codex_totals['Treatment'] == condition_name_codex]
    print(df_condition.shape)
    # calculate mean and stdev over all lymph nodes
    mean_count = df_condition['count'].mean()
    std_count = df_condition['count'].std()
    # format as string with 10^6 (currently in units of cells) with 2 decimal places
    mean_std_str = f'{mean_count/1e3:.2f} (±{std_count/1e3:.2f})'
    table1_df.at[codex_counts_row, cohort] = mean_std_str
# calculate total column for this row
mean_total = df_codex_totals['count'].mean()
std_total = df_codex_totals['count'].std()
mean_std_total_str = f'{mean_total/1e3:.2f} (±{std_total/1e3:.2f})'
table1_df.at[codex_counts_row, 'Total'] = mean_std_total_str
# %%save updated table 1 transposed
table1_df.T.to_csv(f'{processed_data_path}/dataset_overview_table_1_with_cell_counts_transposed_jan2026.csv')
