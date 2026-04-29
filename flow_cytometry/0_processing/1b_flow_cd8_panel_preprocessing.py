#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
# '''=================================================
# @Project: ACT_melanoma_neutrophil (data: act_neutrophil_data_repository)
# @Author : Gemma van der Voort
# @Desc   : preprocessing the cd8-panel of the flow cytometry data, with traditional gating done by Maike.
# @Desc updated: Preprocess CD8 panel flow cytometry data with Maike's traditional gating.
# for recreating figures and double-checking stat. testing results.
# '''=================================================
import os
import sys
from pathlib import Path
import pandas as pd
import warnings
import missingno as msno
import matplotlib.pyplot as plt

PROJECT_ROOT = Path(__file__).resolve()
while PROJECT_ROOT.name != "ACT_melanoma_neutrophil" and PROJECT_ROOT.parent != PROJECT_ROOT:
    PROJECT_ROOT = PROJECT_ROOT.parent
if PROJECT_ROOT.name != "ACT_melanoma_neutrophil":
    raise RuntimeError("ACT_melanoma_neutrophil root not found")
sys.path.append(str(PROJECT_ROOT / "helper_files"))

from paths_parameters import experiments_flow_noCpG, conditions_flow_no_CpG, data_repo_path
from functions_data_preprocessing import preprocess_cd8_panel_flow_endogenous, preprocess_cd8_panel_flow_pmel

plt.style.use(str(PROJECT_ROOT / "helper_files" / "journal_style.mplstyle"))
raw_data_path = f'{data_repo_path}/raw/flow_cytometry/cd8_panel'
processed_data_path = f'{data_repo_path}/processed/flow_cytometry'
fig_path = f'{data_repo_path}/figures/flow_cytometry/preprocessing'
os.makedirs(fig_path, exist_ok=True)
dirs = ['SLO', 'Tumour']
prism_file_path = f'{data_repo_path}/raw/flow_cytometry/cd8_panel_prism_export'
prism_files = ['Pmel Tcm+ [of all live Pmel+].txt', 'Pmel Teff+ [of all live Pmel+].txt']
no_cpg_file = 'flow_cd8_noCpG_incl_tumor.csv'  # from 0c_flow_noCpG_cd8_preprocessing.py, correct format already.
warnings.filterwarnings('ignore')

# traverse both subdirs and list all files
slo_files = os.listdir(f'{raw_data_path}/{dirs[0]}')
slo_files.sort()
tumour_files = os.listdir(f'{raw_data_path}/{dirs[1]}')
tumour_files.sort()
# %%process all files
complete_df = pd.DataFrame(columns=['Experiment', 'File', 'Condition', 'Organ', 'Mouse_ID', 'Parent',
                                     'Teff_count', 'Teff_freq_CD45', 'Teff_freq_endo', 'Teff_freq_pmel',
                                     'Tcm_count', 'Tcm_freq_CD45', 'Tcm_freq_endo', 'Tcm_freq_pmel'])
experiment_to_condition = dict(zip(experiments_flow_noCpG, conditions_flow_no_CpG))
row_count = 0
for file in slo_files + tumour_files:
    if not file.endswith('.xls'):
        continue
    base_name = file.split('.')[:-1]  # remove the .xls
    base_name = '.'.join(base_name)
    print(f'{base_name}')
    # whether the dir is SLO or Tumour
    if file in slo_files:
        dir_type = dirs[0]
    else:
        dir_type = dirs[1]
    file_location = f'{raw_data_path}/{dir_type}/{file}'
    # check how many sheets slo_file has
    slo_df_sheets = pd.ExcelFile(file_location)
    sheets_list = slo_df_sheets.sheet_names
    if len(sheets_list) == 1:
        sheet = sheets_list[0]
        type = 'endogenous'
        print(type)
        df_slim = preprocess_cd8_panel_flow_endogenous(file_location, sheet)
    else:
        for sheet in sheets_list:
            sheet_stripped = sheet.strip()  # strip starting and trailing whitespaces
            if sheet_stripped == 'Endo CD8 T' or sheet_stripped == 'CD8 T' or sheet_stripped == 'Endo CD8':
                type = 'endogenous'
                print(type)
                df_slim = preprocess_cd8_panel_flow_endogenous(file_location, sheet)
            elif sheet_stripped == 'Pmel T' or sheet_stripped == 'Pmel':
                type = 'pmel'
                print(type)
                df_slim = preprocess_cd8_panel_flow_pmel(file_location, sheet)

            else:
                print(f'Unknown sheet: {sheet}')
    row_count += len(df_slim)
    df_slim.to_csv(f'{processed_data_path}/{base_name}_{type}.csv')
    # if file starts with '0', drop 0
    base_new_pruned = base_name[1:] if base_name.startswith('0') else base_name
    experiment_number = base_new_pruned[:2]
    # use mapping of exp number to condition to define condition
    condition = experiment_to_condition[experiment_number]
    df_slim['Experiment'] = experiment_number
    df_slim['File'] = base_name
    df_slim['Condition'] = condition
    df_slim['Parent'] = type
    # add to complete df
    complete_df = pd.concat([complete_df, df_slim])

#%% now add the prism files to the bottom. Experiment: 102. File: one of two prism files, Condition: ACT_d14
# Organ: by header, only 3 LNs. Mouse_ID: 1st column, ".fcs" removed. Parent is 'pmel' in Teff_freq_pmel and Tcm_freq_pmel
# read .txt files as dataframes
df_prism_teff_raw = pd.read_csv(f'{prism_file_path}/{prism_files[1]}', sep='\t')
df_prism_teff = df_prism_teff_raw.rename(columns={'Unnamed: 0': 'Mouse_ID', 'inLNr\r(tdLN)': 'inLNr'})
df_prism_teff['Mouse_ID'] = df_prism_teff['Mouse_ID'].str.replace('.fcs', '', regex=False)
# now melt the df to long format, new column 'Organ' with the organ names, value column 'Teff_freq_pmel'
df_prism_teff = df_prism_teff.melt(id_vars=['Mouse_ID'],
                                      var_name='Organ', value_name='Teff_freq_pmel')
# one value for 2366 inLNr missing, drop this
df_prism_teff = df_prism_teff.dropna(subset=['Teff_freq_pmel'])
# transform from string to float
df_prism_teff['Teff_freq_pmel'] = df_prism_teff['Teff_freq_pmel'].replace(',', '.', regex=True).astype(float)
# %%same for Tcm
df_prism_tcm_raw = pd.read_csv(f'{prism_file_path}/{prism_files[0]}', sep='\t')
df_prism_tcm = df_prism_tcm_raw.rename(columns={'Unnamed: 0': 'Mouse_ID', 'inLNr\r(tdLN)': 'inLNr'})
df_prism_tcm['Mouse_ID'] = df_prism_tcm['Mouse_ID'].str.replace('.fcs', '', regex=False)
df_prism_tcm = df_prism_tcm.melt(id_vars=['Mouse_ID'],
                                        var_name='Organ', value_name='Tcm_freq_pmel')
df_prism_tcm = df_prism_tcm.dropna(subset=['Tcm_freq_pmel'])
# transform from string to float
df_prism_tcm['Tcm_freq_pmel'] = df_prism_tcm['Tcm_freq_pmel'].replace(',', '.', regex=True).astype(float)
# %%unique mouse IDs in complete_df exp 97 and organ brLNr
exp_97_mice = complete_df[complete_df['Experiment'] == '97']
exp_97_mice = exp_97_mice[exp_97_mice['Organ'] == 'brLNr']
exp_97_mice = exp_97_mice['Mouse_ID'].unique()
#%% merge the prism dataframes with each other: add Tcm_freq_pmel to Teff_freq_pmel
df_prism = pd.merge(df_prism_teff, df_prism_tcm, on=['Mouse_ID', 'Organ'])
df_prism = df_prism[~df_prism['Mouse_ID'].isin(exp_97_mice)]
df_prism['Experiment'] = '102'
df_prism['File'] = 'prism export'
df_prism['Condition'] = 'ACT-d14'
df_prism['Parent'] = 'pmel'
#%% merge the prism data with the complete_df
complete_df_2 = pd.concat([complete_df, df_prism], ignore_index=True)
#%% now merge the noCpG data with the complete_df
df_noCpG = pd.read_csv(f'{processed_data_path}/{no_cpg_file}')
df_noCpG = df_noCpG.drop(columns=['Unnamed: 0'])
df_noCpG['Teff_freq_pmel'] = df_noCpG['Teff_freq_pmel'].str.replace(',','.').astype(float)
df_noCpG['Tcm_freq_pmel'] = df_noCpG['Tcm_freq_pmel'].str.replace(',','.').astype(float)
# raname 'ACT-d14_noCpG' to 'ACT-d14_no CpG' and 'ACT-d7_noCpG' to 'ACT-d7_no CpG'
df_noCpG['Condition'] = df_noCpG['Condition'].str.replace('ACT-d14_noCpG', 'ACT-d14_no CpG')
df_noCpG['Condition'] = df_noCpG['Condition'].str.replace('ACT-d7_noCpG', 'ACT-d7_no CpG')
complete_df_3 = pd.concat([complete_df_2, df_noCpG], ignore_index=True)
#%% show df
reduced_matrix = msno.matrix(complete_df_2)
plt.savefig(os.path.join(fig_path, "data_sparsity_matrix_flow_cd8_panel.pdf"), bbox_inches='tight')
plt.close()
#%% count nans and give me the row
complete_df_3 = complete_df_3.reset_index()  # reindex first
#%%
nans = complete_df_3.isna()
# print row where 'Tcm_freq_CD45' is nan
print(complete_df_3[nans['Tcm_freq_CD45']])
# drop sample at index 147
complete_df = complete_df_3.drop(index=147, columns=['index'])

reduced_matrix = msno.matrix(complete_df)
plt.savefig(os.path.join(fig_path, "data_sparsity_matrix_flow_cd8_panel_2256_slo_dropped.pdf"), bbox_inches='tight')
#%% update: rename BAL-4028 in FLOW to BAL-4058 for the lymph nodes
complete_df.loc[(complete_df['Organ'].isin(['inLNr', 'inLNl', 'brLNr'])) & (complete_df['Mouse_ID'] == 'BAL-4028'),
                'Mouse_ID'] = 'BAL-4058'
#%% final check of unique mouse IDs per condition and organ
# save to file
complete_df.to_csv(f'{processed_data_path}/flow_cd8_panel_incl_tumor.csv')

#%%
df = complete_df[complete_df['Parent'] == 'pmel']
df_d14 = df[df['Condition'] == 'ACT-d14']
# inLNl and inLNr only
df_d14 = df_d14[df_d14['Organ'].isin(['inLNr', 'inLNl'])]
# teff only 6 lines, tcm a bit more. Where are the rest of the pores?
# unique mouse IDs
unique_mouse_ids = df_d14['Mouse_ID'].unique()  # 14
#%% how often do they occur?
df_d14['Mouse_ID'].value_counts()
