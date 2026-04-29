#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
# '''=================================================
# @Project: ACT_melanoma_neutrophil (data: act_neutrophil_data_repository)
# @Author : Gemma van der Voort
# @Desc   : flow data for the cpg-no cpg conditions. streamline to process all files from the flowjo exports.
# @Desc updated: Streamlined preprocessing for no-CpG CD8 flow exports.
# '''=================================================
import os
import sys
from pathlib import Path

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve()
while PROJECT_ROOT.name != "ACT_melanoma_neutrophil" and PROJECT_ROOT.parent != PROJECT_ROOT:
    PROJECT_ROOT = PROJECT_ROOT.parent
if PROJECT_ROOT.name != "ACT_melanoma_neutrophil":
    raise RuntimeError("ACT_melanoma_neutrophil root not found")
sys.path.append(str(PROJECT_ROOT / "helper_files"))

from functions_data_preprocessing import split_filename, normalise_organ_names
from paths_parameters import data_repo_path


def extract_cd8_flow_data(filepath, experiment, condition_assigner, teff_col, tcm_col,
                          file_col='Unnamed: 0'):
    df_raw = pd.read_csv(filepath)

    df = pd.DataFrame(columns=['Experiment', 'File', 'Condition', 'Organ', 'Mouse_ID', 'Parent',
                               'Teff_count', 'Teff_freq_CD45', 'Teff_freq_endo', 'Teff_freq_pmel',
                               'Tcm_count', 'Tcm_freq_CD45', 'Tcm_freq_endo', 'Tcm_freq_pmel'])

    df['File'] = df_raw[file_col]
    df['Experiment'] = experiment
    df[['Organ', 'Mouse_ID']] = df['File'].apply(split_filename)
    df['Organ'] = df['Organ'].str.split('_').str[0]
    df['Organ'] = normalise_organ_names(df['Organ'])

    # Print dropped rows with NaN Organ
    dropped = df[df['Organ'].isna()]
    if not dropped.empty:
        print(f"[{experiment}] Dropped rows due to NaN Organ field:")
        for fname in dropped['File']:
            print(f" - {fname}")
    df = df.dropna(subset=['Organ'])

    # Assign condition: function or fixed string
    if callable(condition_assigner):
        df['Condition'] = df['Mouse_ID'].apply(condition_assigner)
    else:
        df['Condition'] = condition_assigner

    df['Parent'] = 'pmel'
    df['Teff_freq_pmel'] = df_raw[teff_col]
    df['Tcm_freq_pmel'] = df_raw[tcm_col]

    return df


def assign_condition_102(mouse_id):
    return 'ACT-d14_no CpG' if mouse_id in no_CpG_d14_mice else 'ACT-d14'


raw_data_path = f'{data_repo_path}/raw/flow_cytometry/early_gates'
processed_data_path = f'{data_repo_path}/processed/flow_cytometry'
files = os.listdir(raw_data_path)

# Define conditions
no_CpG_d14_mice = ['HOE-2479', 'HOE-2489', 'HOE-2474', 'HOE-2482', 'HOE-2488', 'HOE-2473', 'HOE-2477', 'HOE-2470', 'HOE-2481']
# Column mappings SLO files
teff_col = 'All cells/Lymphocytes/Single Cells/Live/CD45.2+/CD8+ Dump-/Pmel 1/Q1: CD62L- , Cd44+ | Freq. of Parent (%)'
tcm_col = 'All cells/Lymphocytes/Single Cells/Live/CD45.2+/CD8+ Dump-/Pmel 1/Q2: CD62L+ , Cd44+ | Freq. of Parent (%)'

# Generate dataframes
df_102 = extract_cd8_flow_data(f'{raw_data_path}/1021_slo_cd8.csv', '102', assign_condition_102,
                               teff_col, tcm_col)
df_102 = df_102[df_102['Condition'] == 'ACT-d14_noCpG']  # drop these because they are already in the bigger dataset.

df_103 = extract_cd8_flow_data(f'{raw_data_path}/1031_slo_cd8.csv', '103', 'ACT-d7_no CpG',
                               teff_col, tcm_col)

# Tumour data
teff_col_tumour = 'All cells/Lymphocytes+Tumour/Single Cells/Live/Immune cells/CD8+ Dump-/Pmel T/Q1: CD62L- , Cd44+ | Freq. of Parent (%)'
tcm_col_tumour = 'All cells/Lymphocytes+Tumour/Single Cells/Live/Immune cells/CD8+ Dump-/Pmel T/Q2: CD62L+ , Cd44+ | Freq. of Parent (%)'

df_102_tumour = extract_cd8_flow_data(f'{raw_data_path}/1021_tumour_cd8.csv', '102', assign_condition_102,
                                      teff_col_tumour, tcm_col_tumour)

df_103_tumour = extract_cd8_flow_data(f'{raw_data_path}/1031_tumour_cd8.csv', '103', 'ACT-d7_no CpG',
                                      teff_col_tumour, tcm_col_tumour)
# also add d14 and d7 tumour. experiments 96 and 97
cd8_96_tumour = pd.read_csv(f'{raw_data_path}/0961_tumour_cd8.csv')
cd8_97_tumour = pd.read_csv(f'{raw_data_path}/0971_tumour_cd8.csv')
#%% for column check (since many spelling variations)
# export cd8_96_tumour columns to text
df_102_raw = pd.read_csv(f'{raw_data_path}/1021_slo_pmel.csv')
with open(f'{processed_data_path}/pmel_102_slo_columns.txt', 'w') as f:
    for col in df_102_raw.columns:
        f.write(f'{col}\n')
#%%
cd8_96_tumour_slim = cd8_96_tumour.loc[:, ~cd8_96_tumour.columns.str.contains(r'\.\d{1,2}$')]
teff_col_tumour_96 = "All cells/Lymphocytes+Tumour cells/Single Cells/Live/CD45.2+ cells/CD8+ T cells/Pmel T cells/Q1: CD62L- , CD44+ | Freq. of Parent (%)"
tcm_col_tumour_96 = "All cells/Lymphocytes+Tumour cells/Single Cells/Live/CD45.2+ cells/CD8+ T cells/Pmel T cells/Q2: CD62L+ , CD44+ | Freq. of Parent (%)"
#%%
df_96_tumour = extract_cd8_flow_data(f'{raw_data_path}/0961_tumour_cd8.csv', '96', 'ACT-d7',
                                      teff_col_tumour_96, tcm_col_tumour_96)
teff_col_tumour_97 = "All cells/Lymphocytes+Tumour cells/Single Cells/Live/CD45.2+ cells/CD8+ T cells/Pmel T cells/Q3: CD44+ , CD62L- | Freq. of Parent (%)"
tcm_col_tumour_97 = "All cells/Lymphocytes+Tumour cells/Single Cells/Live/CD45.2+ cells/CD8+ T cells/Pmel T cells/Q2: CD44+ , CD62L+ | Freq. of Parent (%)"
df_97_tumour = extract_cd8_flow_data(f'{raw_data_path}/0971_tumour_cd8.csv', '97', 'ACT-d14',
                                        teff_col_tumour_97, tcm_col_tumour_97)
#%% merge the dataframes and save
df_noCpG_cd8 = pd.concat([df_102, df_103, df_102_tumour, df_103_tumour, df_96_tumour, df_97_tumour])
df_noCpG_cd8.to_csv(f'{processed_data_path}/flow_cd8_noCpG_incl_tumor.csv')
# this now goes into 1b_flow_cd8_panel_preprocessing.py (line 121)
