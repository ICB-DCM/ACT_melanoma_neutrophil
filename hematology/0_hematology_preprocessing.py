#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
# '''=================================================
# @Project: ACT_melanoma_neutrophil (data: act_neutrophil_data_repository)
# @Author : Gemma van der Voort
# @Desc   : recreate the lymphocyte/neutrophil plot for the blood data from the melanoma project. We might find
# @Desc updated: Preprocess hematology (Hemavet) blood data for downstream plots.
# other interesting results as well.
# update sep 2025: the data contains measurements from cohort 80, clyclo d4. Remove these.
# update 20260109: include again, since they are valid measurements.
# '''=================================================
# imports
import sys
from pathlib import Path

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

PROJECT_ROOT = Path(__file__).resolve()
while PROJECT_ROOT.name != "ACT_melanoma_neutrophil" and PROJECT_ROOT.parent != PROJECT_ROOT:
    PROJECT_ROOT = PROJECT_ROOT.parent
if PROJECT_ROOT.name != "ACT_melanoma_neutrophil":
    raise RuntimeError("ACT_melanoma_neutrophil root not found")
sys.path.append(str(PROJECT_ROOT / "helper_files"))

from paths_parameters import data_repo_path
# paths etc
data_path = f'{data_repo_path}/raw/hematology'
file = 'Hemavet data_FINAL_20240220.xlsx'
sheets = ['WBC_pooled relapse', 'Neu_pooled relapse', 'Lym_pooled relapse']

data_output_path = f'{data_repo_path}/processed/hematology'
fig_path = f'{data_repo_path}/figures/hematology'
modifier = 'data_corrections_jan2026'

# read in the data
wbc = pd.read_excel(f'{data_path}/{file}', sheet_name=sheets[0])
neut = pd.read_excel(f'{data_path}/{file}', sheet_name=sheets[1])
lym = pd.read_excel(f'{data_path}/{file}', sheet_name=sheets[2])
# for each df, drop the final few rows (empty/containing summary info)
wbc = wbc.drop(wbc.tail(3).index)
neut = neut.drop(neut.tail(2).index)
lym = lym.drop(lym.tail(2).index)
# rename the first column from to 'Mouse_ID' in each df
for df in [wbc, neut, lym]:
    df.rename(columns={'Unnamed: 0': 'Mouse_ID'}, inplace=True)
# typo in HOE:-2220 in neut. correct to HOE-2220
neut['Mouse_ID'] = neut['Mouse_ID'].replace('HOE:-2220', 'HOE-2220')
# both dfs contain duplicates for mice 2239,2240, 2241 and 2242. Drop the duplicates
neut = neut.drop_duplicates(subset='Mouse_ID')
lym = lym.drop_duplicates(subset='Mouse_ID')
#%% visualise: plot neutrophils and lymphocytes per condition.
# stack the data to plot it
neut = neut.set_index('Mouse_ID')
lym = lym.set_index('Mouse_ID')
neut_stacked = neut.stack().reset_index()  # same #of rows as print(neut.count().sum())
lym_stacked = lym.stack().reset_index()
# rename the columns
neut_stacked.columns = ['Mouse_ID', 'Condition', 'Neutrophil count (10^9/L)']  # units: Neu # (10^9/L)
lym_stacked.columns = ['Mouse_ID', 'Condition', 'Lymphocyte count (10^9/L)']  # units: Lym # (10^9/L)
# create a shared dataframe, columns Mouse_ID, Condition, Neutrophil count, Lymphocyte count
neut_lym_stacked = pd.merge(neut_stacked, lym_stacked, on=['Mouse_ID', 'Condition'],  validate='1:1')
# export as df
neut_lym_stacked.to_csv(f'{data_output_path}/neut_lym_stacked_{modifier}.csv')
