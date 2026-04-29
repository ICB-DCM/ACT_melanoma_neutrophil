#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
# '''=================================================
# @Project: ACT_melanoma_neutrophil (data: act_neutrophil_data_repository)
# @Author : Gemma van der Voort
# @Desc   : preprocessing of the tumor growth curve data from individual Excel sheets per experiment to one csv.
# @Desc updated: Preprocess tumor growth curves into per-experiment CSVs.
# most have a clear table starting around row 80, but the 102 and 103 need to be parsed separately and have their area
# calculated from L and W columns.
# update 26.09: data correction updates. see data_consistency_check_aug2025.ods/obsidian for details.
# '''=================================================
import sys
from pathlib import Path

import pandas as pd
import numpy as np
import os
from functions_tumor_growth_preprocessing import calculate_area_from_L_W

PROJECT_ROOT = Path(__file__).resolve()
while PROJECT_ROOT.name != "ACT_melanoma_neutrophil" and PROJECT_ROOT.parent != PROJECT_ROOT:
    PROJECT_ROOT = PROJECT_ROOT.parent
if PROJECT_ROOT.name != "ACT_melanoma_neutrophil":
    raise RuntimeError("ACT_melanoma_neutrophil root not found")
sys.path.append(str(PROJECT_ROOT / "helper_files"))

from paths_parameters import data_repo_path

raw_data_path = f'{data_repo_path}/raw/tumor_growth'
processed_data_path = f'{data_repo_path}/processed/tumor_growth'
modifier = 'data_corrections_sep2025'  # to add to the filename when saving
# list all files in the raw data path alphabetically
raw_files = os.listdir(raw_data_path)
raw_files.sort()
# %%open the 1st Excel file to explore
df_70_1 = pd.read_excel(f'{raw_data_path}/070.1_Growth curves.xlsx', sheet_name='Tumour size in mm')
df_70_1 = df_70_1.iloc[77:93, 7:33]
df_70_1.columns = df_70_1.iloc[2] # set the 3nd row as the header
# rename the first column from nan to 'ID' (was on the 2nd row)
df_70_1.rename(columns={np.nan: 'ID'}, inplace=True)
df_70_1 = df_70_1.drop([77, 78, 79]) # drop the first 3 rows
df_70_1.set_index('ID', inplace=True) # set ID as index column
#%% open the 2nd Excel file to explore
df_70_2 = pd.read_excel(f'{raw_data_path}/070.2_Growth curves.xlsx', sheet_name='Tumour size in mm')
df_70_2_part = df_70_2.iloc[70:84, 7:41]
# repeat for 70_2
df_70_2_part.columns = df_70_2_part.iloc[2]
df_70_2_part.rename(columns={np.nan: 'ID'}, inplace=True)
df_70_2_part = df_70_2_part.drop([70, 71, 72])
df_70_2_part.set_index('ID', inplace=True)
#%% open the 3rd Excel file to explore
df_76_1 = pd.read_excel(f'{raw_data_path}/076.1_Groth curves.xlsx', sheet_name='Tumour size in mm')
df_76_1_part = df_76_1.iloc[58:71, 8:35]
# repeat for 76_1
df_76_1_part.columns = df_76_1_part.iloc[2]
df_76_1_part.rename(columns={np.nan: 'ID'}, inplace=True)
df_76_1_part = df_76_1_part.drop([58, 59, 60])
df_76_1_part.set_index('ID', inplace=True)
#%% open the 4th Excel file to explore
df_76_2 = pd.read_excel(f'{raw_data_path}/076.2_Growth curves.xlsx', sheet_name='Tumour size in mm')
df_76_2_part = df_76_2.iloc[21:27, 8:22]
# repeat for 76_2
df_76_2_part.columns = df_76_2_part.iloc[2]
df_76_2_part.rename(columns={np.nan: 'ID'}, inplace=True)
df_76_2_part = df_76_2_part.drop([21, 22, 23])
df_76_2_part.set_index('ID', inplace=True)
#%% open the 5th Excel file to explore
df_76_3 = pd.read_excel(f'{raw_data_path}/076.3_Growth curves.xlsx', sheet_name='Tumour size in mm')
df_76_3_part = df_76_3.iloc[85:98, 7:45]
# repeat for 76_3
df_76_3_part.columns = df_76_3_part.iloc[2]
df_76_3_part.rename(columns={np.nan: 'ID'}, inplace=True)
df_76_3_part = df_76_3_part.drop([85, 86, 87])
df_76_3_part.set_index('ID', inplace=True)
#%% open file 79.1
df_79_1 = pd.read_excel(f'{raw_data_path}/079.1_Growth curves.xlsx', sheet_name='Tumour size in mm')
df_79_1_part = df_79_1.iloc[68:79, 8:39]
# repeat for 79_1
df_79_1_part.columns = df_79_1_part.iloc[2]
df_79_1_part.rename(columns={np.nan: 'ID'}, inplace=True)
df_79_1_part = df_79_1_part.drop([68,69,70])
df_79_1_part.set_index('ID', inplace=True)
# drop mice BAL-3970 and BAL-3971 (mice developed preputial gland abscesses and were sacrificed)
# df_79_1_part = df_79_1_part.drop(['BAL-3970', 'BAL-3971'], axis=0) (removed from the summary table already)
#%% open file 79.2
df_79_2 = pd.read_excel(f'{raw_data_path}/079.2_Growth curves.xlsx', sheet_name='Tumour size in mm')
df_79_2_part = df_79_2.iloc[75:87, 7:45]
# repeat for 79_2
df_79_2_part.columns = df_79_2_part.iloc[2]
df_79_2_part.rename(columns={np.nan: 'ID'}, inplace=True)
df_79_2_part = df_79_2_part.drop([75, 76, 77])
df_79_2_part.set_index('ID', inplace=True)
# drop mouse HOE-2220 (this one was moved to a different experiment due to having been sacrificed too early)
# update: instead of dropping, move to df_70. do this in the 2_ script (since then the 70 group is merged and columns
# match) so keep here for now.
# df_79_2_part = df_79_2_part.drop('HOE-2220', axis=0)
#%% open the 86.1 Excel file to explore
df_86_1 = pd.read_excel(f'{raw_data_path}/086.1_Growth curves.xlsx', sheet_name='Tumour size in mm')
df_86_1_part = df_86_1.iloc[100:117, 7:48]
# repeat for 86_1
df_86_1_part.columns = df_86_1_part.iloc[2]
df_86_1_part.rename(columns={np.nan: 'ID'}, inplace=True)
df_86_1_part = df_86_1_part.drop([100, 101, 102])
df_86_1_part.set_index('ID', inplace=True)
#%% open the 87.1 Excel file to explore
df_87_1 = pd.read_excel(f'{raw_data_path}/087.1_Growth curves.xlsx', sheet_name='Tumour size in mm')
df_87_1_part_overview = df_87_1.iloc[270:297, 7:43]  # this sheet has some empty columns with and without title, double check if
# they have been filled out correctly
# %% not all timepoints have been calculated, so use base L and W for this sheet.
df_87_1_part_1 = df_87_1.iloc[2:30, :22]
df_87_1_part_2 = df_87_1.iloc[32:60, :26]
df_87_1_part_2.at[32, 'Unnamed: 6'] = 'Sex'
df_87_1_part_3 = df_87_1.iloc[61:90, :28]
df_87_1_part_3.at[61, 'Unnamed: 6'] = 'Sex'
df_87_1_part_4 = df_87_1.iloc[90:118, :26]
# add 'L' to cell row 93, col unnamed:24 (missing in sheet)
df_87_1_part_4.at[93, 'Unnamed: 24'] = 'L'
df_87_1_part_4.at[90, 'Unnamed: 6'] = 'Sex'
# the whole Cell line injected column is missing from 4, copy over from 1
df_87_1_part_4.loc[90:118, 'Unnamed: 7'] = df_87_1_part_1.loc[2:30, 'Unnamed: 7'].values  # reset index to match
df_87_1_part_5 = df_87_1.iloc[119:147, :28]
# same for 5 and 6 and 7
df_87_1_part_5.at[119, 'Unnamed: 6'] = 'Sex'
df_87_1_part_5.loc[119:147, 'Unnamed: 7'] = df_87_1_part_1.loc[2:30, 'Unnamed: 7'].values  # reset index to match
df_87_1_part_6 = df_87_1.iloc[148:176, :28]
df_87_1_part_6.at[148, 'Unnamed: 6'] = 'Sex'
df_87_1_part_6.loc[148:176, 'Unnamed: 7'] = df_87_1_part_1.loc[2:30, 'Unnamed: 7'].values  # reset index to match
# add 'L' to cell row 151, col unnamed:8, 10, 12 etc(missing in sheet)
missing_cols = [8, 10, 12, 14, 16, 18, 20, 22, 24, 26]
for col in missing_cols:
    df_87_1_part_6.at[151, f'Unnamed: {col}'] = 'L'
df_87_1_part_7 = df_87_1.iloc[177:205, :26]
df_87_1_part_7.at[177, 'Unnamed: 6'] = 'Sex'
df_87_1_part_7.loc[177:205, 'Unnamed: 7'] = df_87_1_part_1.loc[2:30, 'Unnamed: 7'].values  # reset index to match
# again missing columns, add 'L' to cell row 180, even columns from 8 to 24
missing_cols = [8, 10, 12, 14, 16, 18, 20, 22, 24]
for col in missing_cols:
    df_87_1_part_7.at[180, f'Unnamed: {col}'] = 'L'
# there is some more data, but it is on just one mouse which didn't develop a tumour (HOE-2287), so we can drop it.
#%% now let's clean and aggregate them.
dfs_87_1 = [df_87_1_part_1, df_87_1_part_2, df_87_1_part_3, df_87_1_part_4, df_87_1_part_5, df_87_1_part_6,
            df_87_1_part_7]
area_dfs_87_1 = []
for index_df, df in enumerate(dfs_87_1):
    area_df = calculate_area_from_L_W(df)
    # set ID as index column
    area_df.set_index('ID', inplace=True)
    area_dfs_87_1.append(area_df)
# %%join the dataframes by ID. keep the shared columns from the 1st df and copy unique columns from the others.
shared_cols = ['Strain', 'DOB', 'Age (wks)', 'Ear notch', 'Sex', 'Cell line injected']
df_87_1_cleaned = area_dfs_87_1[0]
for i in range(1, len(area_dfs_87_1)):
    df_87_1_cleaned = df_87_1_cleaned.join(area_dfs_87_1[i].drop(columns=shared_cols), how='outer')

#%% open the 96.1 Excel file to explore
df_96_1 = pd.read_excel(f'{raw_data_path}/096.1_Growth curves.xlsx', sheet_name='Tumour size in mm')
# no area calculated, so use base L and W for this sheet.
df_96_1_part_1 = df_96_1.iloc[2:21, :22]
df_96_1_part_2 = df_96_1.iloc[22:41, :22]
df_96_1_part_3 = df_96_1.iloc[42:61, :22]
df_96_1_part_3 = df_96_1_part_3.replace('?', np.nan)
df_96_1_part_4 = df_96_1.iloc[62:81, :22]
df_96_1_part_5 = df_96_1.iloc[82:101, :22]
df_96_1_part_6 = df_96_1.iloc[102:121, :18]
# %%calculating the area.
dfs_96_1 = [df_96_1_part_1, df_96_1_part_2, df_96_1_part_3, df_96_1_part_4, df_96_1_part_5, df_96_1_part_6]
area_dfs_96_1 = []
for index_df, df in enumerate(dfs_96_1):
    area_df = calculate_area_from_L_W(df)
    area_df.set_index('ID', inplace=True)  # set ID as index column
    area_dfs_96_1.append(area_df)
# %%join the dataframes by ID. keep the shared columns from the 1st df and copy unique columns from the others.
df_96_1_cleaned = area_dfs_96_1[0]
for i in range(1, len(area_dfs_96_1)):
    df_96_1_cleaned = df_96_1_cleaned.join(area_dfs_96_1[i].drop(columns=shared_cols), how='outer')

#%% open the 97.1 Excel file to explore
df_97_1 = pd.read_excel(f'{raw_data_path}/097.1_Growth curves.xlsx', sheet_name='Tumour size in mm')
# df_97_1_part = df_97_1.iloc[100:117, 7:39]
# no area calculated, so use base L and W for this sheet.
df_97_1_part_1 = df_97_1.iloc[2:21, :22]
df_97_1_part_2 = df_97_1.iloc[22:41, :22]
df_97_1_part_3 = df_97_1.iloc[42:61, :24]
# replace ? with NaN
df_97_1_part_3 = df_97_1_part_3.replace('?', np.nan)
df_97_1_part_4 = df_97_1.iloc[62:81, :24]
df_97_1_part_5 = df_97_1.iloc[82:101, :26]
# there are two more columns for mice that never developed a tumour, omitting those.
# %%calculating the area.
dfs_97_1 = [df_97_1_part_1, df_97_1_part_2, df_97_1_part_3, df_97_1_part_4, df_97_1_part_5]
area_dfs_97_1 = []
for index_df, df in enumerate(dfs_97_1):
    area_df = calculate_area_from_L_W(df)
    area_df.set_index('ID', inplace=True)  # set ID as index column
    area_dfs_97_1.append(area_df)
# %%join the dataframes by ID. keep the shared columns from the 1st df and copy unique columns from the others.
df_97_1_cleaned = area_dfs_97_1[0]
for i in range(1, len(area_dfs_97_1)):
    df_97_1_cleaned = df_97_1_cleaned.join(area_dfs_97_1[i].drop(columns=shared_cols), how='outer')
#%% open the 102.1 Excel file to explore
# different format, so I modified the sheet a bit to fit the rest (skipping 9th, original version)
df_102_1 = pd.read_excel(f'{raw_data_path}/102.1_Growth curves.xlsx', sheet_name='Tumour size in mm')
# calculate area from L and W columns, only partially filled out.
df_102_1_part_1 = df_102_1.iloc[2:26, :24]
df_102_1_part_2 = df_102_1.iloc[27:51, :24]
df_102_1_part_3 = df_102_1.iloc[52:76, :24]
df_102_1_part_4 = df_102_1.iloc[77:101, :24]
df_102_1_part_5 = df_102_1.iloc[102:126, :24]
df_102_1_part_6 = df_102_1.iloc[127:151, :14]
# %%calculating the area.
dfs_102_1 = [df_102_1_part_1, df_102_1_part_2, df_102_1_part_3, df_102_1_part_4, df_102_1_part_5, df_102_1_part_6]
area_dfs_102_1 = []
for index_df, df in enumerate(dfs_102_1):
    df = df.replace('x', np.nan)  # this sheet has some x's indicating CpG treatment, not a tumour size measurement
    df = df.replace('?', np.nan)
    area_df = calculate_area_from_L_W(df)
    area_df.set_index('ID', inplace=True)  # set ID as index column
    area_dfs_102_1.append(area_df)
# %%join the dataframes by ID. keep the shared columns from the 1st df and copy unique columns from the others.
df_102_1_cleaned = area_dfs_102_1[0]
for i in range(1, len(area_dfs_102_1)):
    df_102_1_cleaned = df_102_1_cleaned.join(area_dfs_102_1[i].drop(columns=shared_cols), how='outer')
#%% open the 103.1 Excel file to explore
df_103_1 = pd.read_excel(f'{raw_data_path}/103.1_Growth curves.xlsx', sheet_name='Tumour size in mm')
# calculate area from L and W columns, only partially filled out.
df_103_1_part_1 = df_103_1.iloc[2:13, :24]
df_103_1_part_2 = df_103_1.iloc[14:25, :24]
df_103_1_part_3 = df_103_1.iloc[26:37, :24]
df_103_1_part_4 = df_103_1.iloc[38:49, :24]
#%% calculating the area.
dfs_103_1 = [df_103_1_part_1, df_103_1_part_2, df_103_1_part_3, df_103_1_part_4]
area_dfs_103_1 = []
for index_df, df in enumerate(dfs_103_1):
    df = df.replace('*', np.nan)  # * to indicate a note, not a measurement
    df = df.replace('?', np.nan)
    area_df = calculate_area_from_L_W(df)
    area_df.set_index('ID', inplace=True)  # set ID as index column
    area_dfs_103_1.append(area_df)
# %%join the dataframes by ID. keep the shared columns from the 1st df and copy unique columns from the others.
df_103_1_cleaned = area_dfs_103_1[0]
for i in range(1, len(area_dfs_103_1)):
    df_103_1_cleaned = df_103_1_cleaned.join(area_dfs_103_1[i].drop(columns=shared_cols), how='outer')

#%% now export all the cleaned df to file in processed data path
df_70_1.to_csv(f'{processed_data_path}/per_experiment/tumor_growth_curves_70_1_{modifier}.csv')
df_70_2_part.to_csv(f'{processed_data_path}/per_experiment/tumor_growth_curves_70_2_{modifier}.csv')
df_76_1_part.to_csv(f'{processed_data_path}/per_experiment/tumor_growth_curves_76_1_{modifier}.csv')
df_76_2_part.to_csv(f'{processed_data_path}/per_experiment/tumor_growth_curves_76_2_{modifier}.csv')
df_76_3_part.to_csv(f'{processed_data_path}/per_experiment/tumor_growth_curves_76_3_{modifier}.csv')
df_79_1_part.to_csv(f'{processed_data_path}/per_experiment/tumor_growth_curves_79_1_{modifier}.csv')
df_79_2_part.to_csv(f'{processed_data_path}/per_experiment/tumor_growth_curves_79_2_{modifier}.csv')
df_86_1_part.to_csv(f'{processed_data_path}/per_experiment/tumor_growth_curves_86_1_{modifier}.csv')
df_87_1_cleaned.to_csv(f'{processed_data_path}/per_experiment/tumor_growth_curves_87_1_{modifier}.csv')
df_96_1_cleaned.to_csv(f'{processed_data_path}/per_experiment/tumor_growth_curves_96_1_{modifier}.csv')
df_97_1_cleaned.to_csv(f'{processed_data_path}/per_experiment/tumor_growth_curves_97_1_{modifier}.csv')
df_102_1_cleaned.to_csv(f'{processed_data_path}/per_experiment/tumor_growth_curves_102_1_{modifier}.csv')
df_103_1_cleaned.to_csv(f'{processed_data_path}/per_experiment/tumor_growth_curves_103_1_{modifier}.csv')
