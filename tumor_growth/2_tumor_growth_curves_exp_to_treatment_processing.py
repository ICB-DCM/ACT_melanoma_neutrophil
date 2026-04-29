#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
# '''=================================================
# @Project: ACT_melanoma_neutrophil (data: act_neutrophil_data_repository)
# @Author : Gemma van der Voort
# @Desc   : The per-experiment preprocessing needs to be split from a per-experiment into a per-treatment preprocessing,
# @Desc updated: Aggregate tumor growth curves by treatment and export metadata.
# since some experiments have multiple treatments and some treatments are split across multiple experiments.
# update 26.09: data correction updates. see data_consistency_check_aug2025.ods/obsidian for details.
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

from paths_parameters import data_repo_path

processed_data_path = f'{data_repo_path}/processed/tumor_growth'
fig_dir = f'{data_repo_path}/figures/tumor_growth_curves'
modifier = 'data_corrections_sep2025'  # to add to the filename when saving

files_growth_data = os.listdir(processed_data_path)
metadata_cols = ['Strain', 'DOB', 'Age (wks)', 'Ear notch', 'Sex', 'Cell line injected']
exp_with_metadata = ['87_1', '96_1', '97_1', '102_1', '103_1']
# %%load the dfs
df_70_1 = pd.read_csv(f'{processed_data_path}/per_experiment/tumor_growth_curves_70_1_{modifier}.csv', index_col=0)
df_70_2 = pd.read_csv(f'{processed_data_path}/per_experiment/tumor_growth_curves_70_2_{modifier}.csv', index_col=0)
df_76_1 = pd.read_csv(f'{processed_data_path}/per_experiment/tumor_growth_curves_76_1_{modifier}.csv', index_col=0)
df_76_2 = pd.read_csv(f'{processed_data_path}/per_experiment/tumor_growth_curves_76_2_{modifier}.csv', index_col=0)
df_76_3 = pd.read_csv(f'{processed_data_path}/per_experiment/tumor_growth_curves_76_3_{modifier}.csv', index_col=0)
df_79_1 = pd.read_csv(f'{processed_data_path}/per_experiment/tumor_growth_curves_79_1_{modifier}.csv', index_col=0)
df_79_2 = pd.read_csv(f'{processed_data_path}/per_experiment/tumor_growth_curves_79_2_{modifier}.csv', index_col=0)
df_86_1 = pd.read_csv(f'{processed_data_path}/per_experiment/tumor_growth_curves_86_1_{modifier}.csv', index_col=0)
df_87_1 = pd.read_csv(f'{processed_data_path}/per_experiment/tumor_growth_curves_87_1_{modifier}.csv', index_col=0)
df_96_1 = pd.read_csv(f'{processed_data_path}/per_experiment/tumor_growth_curves_96_1_{modifier}.csv', index_col=0)
df_97_1 = pd.read_csv(f'{processed_data_path}/per_experiment/tumor_growth_curves_97_1_{modifier}.csv', index_col=0)
df_102_1 = pd.read_csv(f'{processed_data_path}/per_experiment/tumor_growth_curves_102_1_{modifier}.csv', index_col=0)
df_103_1 = pd.read_csv(f'{processed_data_path}/per_experiment/tumor_growth_curves_103_1_{modifier}.csv', index_col=0)
# drop metadata columns from the dfs that have them
df_87_1 = df_87_1.drop(columns=metadata_cols)
df_96_1 = df_96_1.drop(columns=metadata_cols)
df_97_1 = df_97_1.drop(columns=metadata_cols)
df_102_1 = df_102_1.drop(columns=metadata_cols)
df_103_1 = df_103_1.drop(columns=metadata_cols)
# %%merge the 70 and 76 dfs. add 70.2 below 70.1, combine any columns that have the same name
df_70 = pd.concat([df_70_1, df_70_2], axis=0)
# sort the columns numerically (so treat them as integers)
df_70 = df_70.reindex(sorted(df_70.columns, key=lambda x: int(x)), axis=1)
df_76 = pd.concat([df_76_1, df_76_2, df_76_3], axis=0)
df_76 = df_76.reindex(sorted(df_76.columns, key=lambda x: int(x)), axis=1)
df_79 = pd.concat([df_79_1, df_79_2], axis=0)
df_79 = df_79.reindex(sorted(df_79.columns, key=lambda x: int(x)), axis=1)
#%% df_102_1 needs to be split:
ACT_d14_102 = ['HOE-2478', 'HOE-2486', 'HOE-2480', 'HOE-2485', 'HOE-2472', 'HOE-2476', 'HOE-2483']
ACT_d14_no_CpG_102 = ['HOE-2479', 'HOE-2489', 'HOE-2474', 'HOE-2482', 'HOE-2488', 'HOE-2473', 'HOE-2477', 'HOE-2470',
                      'HOE-2481']
# find indexes in ACT_d14_102 in df_102_1 and split the df
df_102_1_ACT_d14 = df_102_1[df_102_1.index.isin(ACT_d14_102)]
df_102_1_ACT_d14_no_CpG = df_102_1[df_102_1.index.isin(ACT_d14_no_CpG_102)]
#%% add the df_102_1_ACT_d14 to df_97_1 (both ACT d14)
df_act_d14 = pd.concat([df_97_1, df_102_1_ACT_d14], axis=0)
# sort the columns numerically (so treat them as integers)
df_act_d14 = df_act_d14.reindex(sorted(df_act_d14.columns, key=lambda x: int(x)), axis=1)
# df_87_1 has an all nan row at the end, remove it
df_87_1 = df_87_1.dropna(how='all')
#%% update: 26.09.2025 data consistency check corrections
# BAL 3605: move from TB3 (70) to TB8 (79)
df_79 = pd.concat([df_79, df_70.loc[df_70.index == 'BAL-3605']], axis=0)
df_70 = df_70.loc[df_70.index != 'BAL-3605']
# %%BAL 3614 misnamed tb in codex, is cyclo d1. check here: OK, is in cyclo here.
print(df_70.index.unique())
print(df_76.index.unique())
# %%mouse HOE-2220: tumor growth to tb 3-5 (was originally 8-10, but culled early). I dropped it from tb3 originally,
# now added it to tb3 in scripts 0 and 1, in order to swap now.
df_70 = pd.concat([df_70, df_79.loc[df_79.index == 'HOE-2220']], axis=0)
df_79 = df_79.loc[df_79.index != 'HOE-2220']
# mice HOE-2263 and 2276 in act d3 had their final growth cut off. correction in 0_, check here:
# sacrifice day and events correct now for these mice, timeseries as well (used df interactive view)
# HOE-2495: injected with pmels at wrong time. Discard growth curve in act d7 noCpG.
df_103_1 = df_103_1.loc[df_103_1.index != 'HOE-2495']
print(df_103_1.index.unique())
 # %%df list
dfs = [df_70, df_76, df_79, df_86_1, df_96_1, df_act_d14, df_87_1, df_103_1, df_102_1_ACT_d14_no_CpG]
titles_tumor_growth = ["Tumor bearing", "Cyclo", "Tumor bearing 8-10mm", "ACT day 3",
                      "ACT day 7", "ACT day 14", "Relapse", "ACT day 7, no CpG", "ACT day 14, no CpG"]
# file names: replace spaces with underscores and remove commas
df_names = [name.replace(' ', '_').replace(',', '') for name in titles_tumor_growth]
# %%save the dfs
for df_index, df in enumerate(dfs):
    df.to_csv(f'{processed_data_path}/per_treatment/{df_names[df_index]}_{modifier}.csv')

#%% repeat for metadata dfs
files_metadata = os.listdir(os.path.join(processed_data_path, 'metadata/per_experiment'))
files_metadata.sort()
# files_metadata_explicit = ['tumor_growth_curves_076_1_meta.csv', 'tumor_growth_curves_076_2_meta.csv',
#                            'tumor_growth_curves_076_3_meta.csv', 'tumor_growth_curves_086_1_meta.csv',
#                            'tumor_growth_curves_087_1_meta.csv', 'tumor_growth_curves_096_1_meta.csv',
#                            'tumor_growth_curves_097_1_meta.csv', 'tumor_growth_curves_102_1_meta.csv',
#                            'tumor_growth_curves_103_1_meta.csv']
# files_metadata_explicit = ['tumor_growth_curves_070_1_meta.csv', 'tumor_growth_curves_070_2_meta.csv',
#                            'tumor_growth_curves_076_1_meta.csv', 'tumor_growth_curves_076_2_meta.csv',
#                            'tumor_growth_curves_076_3_meta.csv', 'tumor_growth_curves_079_1_meta.csv',
#                            'tumor_growth_curves_079_2_meta.csv', 'tumor_growth_curves_086_1_meta.csv',
#                            'tumor_growth_curves_087_1_meta.csv', 'tumor_growth_curves_096_1_meta.csv',
#                            'tumor_growth_curves_097_1_meta.csv', 'tumor_growth_curves_102_1_meta.csv',
#                            'tumor_growth_curves_103_1_meta.csv']

# %%load the metadata dfs
df_70_1_meta = pd.read_csv(f'{processed_data_path}/metadata/per_experiment/tumor_growth_curves_070_1_meta_{modifier}.csv', index_col=0)
df_70_2_meta = pd.read_csv(f'{processed_data_path}/metadata/per_experiment/tumor_growth_curves_070_2_meta_{modifier}.csv', index_col=0)
df_76_1_meta = pd.read_csv(f'{processed_data_path}/metadata/per_experiment/tumor_growth_curves_076_1_meta_{modifier}.csv', index_col=0)
df_76_2_meta = pd.read_csv(f'{processed_data_path}/metadata/per_experiment/tumor_growth_curves_076_2_meta_{modifier}.csv', index_col=0)
df_76_3_meta = pd.read_csv(f'{processed_data_path}/metadata/per_experiment/tumor_growth_curves_076_3_meta_{modifier}.csv', index_col=0)
df_79_1_meta = pd.read_csv(f'{processed_data_path}/metadata/per_experiment/tumor_growth_curves_079_1_meta_{modifier}.csv', index_col=0)
df_79_2_meta = pd.read_csv(f'{processed_data_path}/metadata/per_experiment/tumor_growth_curves_079_2_meta_{modifier}.csv', index_col=0)
df_86_1_meta = pd.read_csv(f'{processed_data_path}/metadata/per_experiment/tumor_growth_curves_086_1_meta_{modifier}.csv', index_col=0)
df_87_1_meta = pd.read_csv(f'{processed_data_path}/metadata/per_experiment/tumor_growth_curves_087_1_meta_{modifier}.csv', index_col=0)
df_96_1_meta = pd.read_csv(f'{processed_data_path}/metadata/per_experiment/tumor_growth_curves_096_1_meta_{modifier}.csv', index_col=0)
df_97_1_meta = pd.read_csv(f'{processed_data_path}/metadata/per_experiment/tumor_growth_curves_097_1_meta_{modifier}.csv', index_col=0)
df_102_1_meta = pd.read_csv(f'{processed_data_path}/metadata/per_experiment/tumor_growth_curves_102_1_meta_{modifier}.csv', index_col=0)
df_103_1_meta = pd.read_csv(f'{processed_data_path}/metadata/per_experiment/tumor_growth_curves_103_1_meta_{modifier}.csv', index_col=0)
#%% 70,76,79 exp: add concat on rows
df_70_meta = pd.concat([df_70_1_meta, df_70_2_meta], axis=0)
df_76_meta = pd.concat([df_76_1_meta, df_76_2_meta, df_76_3_meta], axis=0)
df_79_meta = pd.concat([df_79_1_meta, df_79_2_meta], axis=0)
#%% 102: split by No CpG and CpG mice (see lists above)
df_102_1_meta_ACT_d14 = df_102_1_meta[df_102_1_meta['mouse_ID'].isin(ACT_d14_102)]
df_102_1_meta_ACT_d14_no_CpG = df_102_1_meta[df_102_1_meta['mouse_ID'].isin(ACT_d14_no_CpG_102)]
#%% add ACT d14 to 97
df_act_d14_meta = pd.concat([df_97_1_meta, df_102_1_meta_ACT_d14], axis=0)
#%% update: 26.09.2025 data consistency check corrections
# BAL 3605: move from TB3 (70) to TB8 (79)
df_79_meta = pd.concat([df_79_meta, df_70_meta[df_70_meta['mouse_ID'] == 'BAL-3605']], axis=0)
df_70_meta = df_70_meta[~(df_70_meta['mouse_ID'] == 'BAL-3605')]
# BAL 3614 misnamed tb in codex, is cyclo d1. check here: OK, is in cyclo here.
print(df_70_meta['mouse_ID'].unique())
print(df_76_meta['mouse_ID'].unique())
# mouse HOE-2220: tumor growth to tb 3-5 (was originally 8-10, but culled early). I dropped it from tb3 originally,
# now added it to tb3 in scripts 0 and 1, in order to swap now.
df_70_meta = pd.concat([df_70_meta, df_79_meta[df_79_meta['mouse_ID'] == 'HOE-2220']], axis=0)
df_79_meta = df_79_meta[~(df_79_meta['mouse_ID'] == 'HOE-2220')]
# mice HOE-2263 and 2276 in act d3 had their final growth cut off. correction in 0_, check here:
# sacrifice day and events correct now for these mice, timeseries as well (used df interactive view)
# HOE-2495: injected with pmels at wrong time. Discard growth curve in act d7 noCpG.
print(df_103_1_meta['mouse_ID'].unique())
df_103_1_meta = df_103_1_meta[~(df_103_1_meta['mouse_ID'] == 'HOE-2495')]
# %%from each of the meta dfs, drop the columns experiment date cell, value, value cell, and sheet. then save
dfs_meta = [df_70_meta, df_76_meta, df_79_meta, df_86_1_meta, df_96_1_meta, df_act_d14_meta, df_87_1_meta, df_103_1_meta,
            df_102_1_meta_ACT_d14_no_CpG]
for df_meta_index, df_meta in enumerate(dfs_meta):
    df_meta = df_meta.drop(columns=['experiment date cell', 'value', 'value cell', 'sheet'])
    df_meta.to_csv(f'{processed_data_path}/metadata/per_treatment/{df_names[df_meta_index]}_{modifier}.csv')

#%% check
df_act_d14_meta['mouse_ID'].value_counts()
