#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
# '''=================================================
# @Project: ACT_melanoma_neutrophil (data: act_neutrophil_data_repository)
# @Author : Gemma van der Voort
# @Desc   : preprocess the pmel t cells in the blood, similar format to the pan-immune data.
# @Desc updated: Preprocess averaged Pmel T blood data for downstream visualization.
# '''=================================================
import sys
from pathlib import Path

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve()
while PROJECT_ROOT.name != "ACT_melanoma_neutrophil" and PROJECT_ROOT.parent != PROJECT_ROOT:
    PROJECT_ROOT = PROJECT_ROOT.parent
if PROJECT_ROOT.name != "ACT_melanoma_neutrophil":
    raise RuntimeError("ACT_melanoma_neutrophil root not found")
sys.path.append(str(PROJECT_ROOT / "helper_files"))

from functions_data_preprocessing import normalise_organ_names
from paths_parameters import data_repo_path

# paths etc
data_path = f'{data_repo_path}/raw/flow_cytometry/pmels_in_blood'
file = 'CD8 T and Pmel T_Freq-Count_Mean 3 panels.xlsx'
sheets = ['LN_Mean', 'Spleen_Mean', 'Blood_Mean', 'Tumour_Mean']

data_output_path = f'{data_repo_path}/processed/flow_cytometry'

# %%read in the data
blood_df = pd.read_excel(f'{data_path}/{file}', sheet_name=sheets[2], skiprows=5, index_col=[0, 1, 2], header=[0])
# %%data preprocessing
# set 2nd index level to a regular column called 'Condition'
blood_df_slim = blood_df.reset_index(level=1)
blood_df_slim = blood_df_slim.rename(columns={'level_1': 'Condition'})
#name the index levels
blood_df_slim.index.names = ['Experiment', 'File']
# %%drop all col except the mean and the condition (1st and last column
blood_df_slim = blood_df_slim.iloc[:, [0, -1]]
# mouse ID column:
blood_df_slim['Mouse_ID'] = blood_df_slim.index.get_level_values(1)
# find the occurrence of either 'HOE' or 'BAL' in the mouse ID and keep that substring plus the 4 characters after it
blood_df_slim['Mouse_ID'] = blood_df_slim['Mouse_ID'].str.extract(r'(HOE|BAL)(.{5})', expand=True).agg(''.join, axis=1)
# does blood_df_slim contains nans
print(blood_df_slim.isnull().sum())
blood_df_slim = blood_df_slim.dropna() # drops non-pmel conditions and an empty value in 96.1.
# save to file
blood_df_slim.to_csv(f'{data_output_path}/pmel_t_blood_means.csv')
