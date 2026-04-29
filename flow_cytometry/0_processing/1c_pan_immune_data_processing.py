#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
# '''=================================================
# @Project: ACT_melanoma_neutrophil (data: act_neutrophil_data_repository)
# @Author : Gemma van der Voort
# @Desc   : preprocess the pan-immune data. Don't drop any cell types or conditions, just normalise,
# @Desc updated: Normalize pan-immune flow data and export CD45-corrected counts.
# drop the failed sample and save the data.
# update 28.05.2025: cd45 true numbers and not sum of phenotyped cd45+ cells (misses macrophages)
# '''=================================================
import os
import sys
from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

PROJECT_ROOT = Path(__file__).resolve()
while PROJECT_ROOT.name != "ACT_melanoma_neutrophil" and PROJECT_ROOT.parent != PROJECT_ROOT:
    PROJECT_ROOT = PROJECT_ROOT.parent
if PROJECT_ROOT.name != "ACT_melanoma_neutrophil":
    raise RuntimeError("ACT_melanoma_neutrophil root not found")
sys.path.append(str(PROJECT_ROOT / "helper_files"))

from paths_parameters import data_repo_path

# paths etc
data_path = f'{data_repo_path}/raw/flow_cytometry/pan_immune'
file = 'Cell Counts_Pan-Immune_20240219.xlsx'
sheets = ['Tumour_Pan-Immune', 'Spleen_Pan-Immune', 'LN_Pan-Immune']


data_output_path = f'{data_repo_path}/processed/flow_cytometry'
early_gates_slo = f'{data_repo_path}/processed/flow_cytometry/early_gates_slo_counts.csv'
fig_path = f'{data_repo_path}/figures/flow_cytometry/preprocessing'
os.makedirs(fig_path, exist_ok=True)
# %%read in the data
ln_df = pd.read_excel(f'{data_path}/{file}', sheet_name=sheets[2], skiprows=5, index_col=[0, 1, 2], header=[0])
# %%data preprocessing
# set 2nd index level to a regular column called 'Condition'
ln_df_slim = ln_df.reset_index(level=1)
ln_df_slim = ln_df_slim.rename(columns={'level_1': 'Condition'})
#name the index levels
ln_df_slim.index.names = ['Experiment', 'File']
# drop all col after index 11
ln_df_slim = ln_df_slim.iloc[:, :12]
# we need another column for the type of LN
ln_df_slim['LN_type'] = ln_df_slim.index.get_level_values(1).str.split('_').str[0]
# and a Mouse_ID column for the final part of the '_' split, with the final 4 characters removed
# some entries have two underscores due to the condition being part of the file name.
ln_df_slim['Mouse_ID'] = ln_df_slim.index.get_level_values(1).str.split('_').str[-1].str[:-4]
# we have some spelling variations, so we need to correct them
ln_df_slim['LN_type'] = ln_df_slim['LN_type'].replace({'brLNri': 'brLNr', 'brLN-r': 'brLNr', 'brLN-right': 'brLNr',
                                                       'inLNle': 'inLNl', 'inLN-l': 'inLNl', 'inLN-left': 'inLNl',
                                                       'inLNri': 'inLNr', 'inLN-r': 'inLNr', 'inLN-right': 'inLNr'})
# normalise by all 45+ positive cells
cols_to_normalise_by = ['DCs | Count normalised', 'CD8 T | Count normalised', 'CD4 Tconv | Count normalised',
                        'CD25+ CD4+ T | Count normalised',  # Tregs
                        'Neutrophils | Count normalised', 'Nk cells | Count normalised',
                        'B cells | Count normalised']
cell_count_columns = ['DCs | Count normalised', 'CD11b+ Dcs | Count normalised',
       'Xcr1+ Dcs | Count normalised', 'CD8 T | Count normalised',
       'Endo CD8 T | Count normalised', 'Pmel T  | Count normalised',
       'CD4 Tconv | Count normalised', 'CD25+ CD4+ T | Count normalised',
       'Neutrophils | Count normalised', 'Nk cells | Count normalised',
       'B cells | Count normalised']
# add column 'total CD45+ cells' to the ln_df_slim which is the sum of all the columns in cols_to_normalise_by
ln_df_slim['total CD45+ cells uncorrected'] = ln_df[cols_to_normalise_by].sum(axis=1).tolist()
# correct data errors:
ln_df_slim = ln_df_slim.drop('inLN-l_HOE-2242.fcs', level=1)  # outlier with only 300 cells
# 22.05.2025: found a typo in mouse BAL-3427, (only brLnr, pan-immune): should be BAL-3727 (only brLNr, pan-immune missing)
ln_df_slim.loc[ln_df_slim['Mouse_ID'] == 'BAL-3427', 'Mouse_ID'] = 'BAL-3727'

# %%calculate actual amount of CD45+ cells (count normalised)
df_early_gates = pd.read_csv(early_gates_slo)
lymph_nodes = ['inLNr', 'inLNl', 'brLNr']  # we are only interested in LNs here
# drop all rows that are not lymph nodes
df_early_gates_ln = df_early_gates[df_early_gates['Organ'].isin(lymph_nodes)]
df_early_gates_ln = df_early_gates_ln[df_early_gates_ln['Panel'] == 'panimmune']
# %%add column bead_factor, 5000 (original amount of beads) / bead_count and organ_factor_LN
df_early_gates_ln['bead_factor'] = 5000 / df_early_gates_ln['bead_count']
df_early_gates_ln['organ_factor'] = 1/ df_early_gates_ln['Organ_fraction_LN']
df_early_gates_ln['count_multiplication_factor'] = df_early_gates_ln['bead_factor'] * df_early_gates_ln['organ_factor']
df_early_gates_ln['CD45+ count normalised'] = df_early_gates_ln['cd45_count'] * df_early_gates_ln['count_multiplication_factor']
#%% check for duplicate df_early_gates_ln File names
duplicates = df_early_gates_ln[df_early_gates_ln.duplicated(subset=['File'])]
# drop them
df_early_gates_ln = df_early_gates_ln.drop_duplicates(subset=['File'])
# match 'File' in df_early_gates_ln with 'File' in ln_df_slim, report any mismatches (should be 4)
# make index col in ln_slim into regular column
ln_df_slim = ln_df_slim.reset_index()
mismatches_slim = ln_df_slim[~ln_df_slim['File'].isin(df_early_gates_ln['File'])]
mismatches_early_gates = df_early_gates_ln[~df_early_gates_ln['File'].isin(ln_df_slim['File'])]
# %%checks out, these are mouse 3605 (analysis issue, no late gating) and mouse 2242 inLNl, with very little cells.
# drop these mice from df_early_gates_ln
df_early_gates_ln = df_early_gates_ln[~df_early_gates_ln['File'].isin(mismatches_early_gates['File'])]
# %% now add the CD45+ count normalised to the ln_df_slim by matching the 'File' column (order not the same)
ln_df_slim = ln_df_slim.merge(df_early_gates_ln[['File', 'CD45+ count normalised']],
                                on='File', how='left')
#%% calculate the percentual increase of 'CD45+ count normalised' compared to 'total CD45+ cells uncorrected' and plot
# the fractions in a histogram
ln_df_slim['percentual_increase_cd45'] = ((ln_df_slim['CD45+ count normalised'] / ln_df_slim['total CD45+ cells uncorrected']) * 100) - 100
fig, ax = plt.subplots(figsize=(8, 6))
sns.histplot(ln_df_slim, x='percentual_increase_cd45', bins=30, kde=True, ax=ax, hue='Condition')
ax.set_xlabel('Percentual increase of CD45+ cells compared to uncorrected total CD45+ cells')
ax.set_ylabel('Frequency')
fig.savefig(f'{fig_path}/percentual_increase_cd45_histogram.pdf')
# %%add percent of the total CD45+ cells to the ln_df_slim
for col in cell_count_columns:
    # for the new name, keep only the part before the '|', and add 'percent' to end
    percent_col_name = col.split('|')[0].strip() + ' percent'
    ln_df_slim[percent_col_name] = (ln_df_slim[col] / ln_df_slim['CD45+ count normalised']) * 100

#%% update: rename BAL-4028 in FLOW to BAL-4058 for the lymph nodes
ln_df_slim.loc[(ln_df_slim['LN_type'].isin(['inLNr', 'inLNl', 'brLNr'])) & (ln_df_slim['Mouse_ID'] == 'BAL-4028'),
                'Mouse_ID'] = 'BAL-4058'  # doesn't warrant rerunning, just for correctness-sake in case of later
# paired analysis with blood, spleen, tumour etc.
# no need to rerun later, since only brLNr was changed, and we don't do paired analysis for that organ.
# save data to file
ln_df_slim.to_csv(f'{data_output_path}/ln_pan_immune_normalised_v2_cd45_corrected.csv')
