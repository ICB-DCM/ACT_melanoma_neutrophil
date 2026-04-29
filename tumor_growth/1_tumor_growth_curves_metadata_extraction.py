#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
# '''=================================================
# @Project: ACT_melanoma_neutrophil (data: act_neutrophil_data_repository)
# @Author : Gemma van der Voort
# @Desc   : the Excel files contain information on at which day treatment steps were taken by colouring of the cells.
# extract this information to a dataframe.
# @Desc updated: Extract treatment-event metadata from tumor growth workbooks.
# update  20250626 also extract sacrifice dates from the short and long tumor growth experiments (70 and 79)
# update 26.09: data correction updates. see data_consistency_check_aug2025.ods/obsidian for details.
# '''=================================================
# imports
import sys
from pathlib import Path

import openpyxl
import os
import pandas as pd
import numpy as np
from openpyxl.utils.cell import column_index_from_string
from datetime import datetime
from functions_tumor_growth_preprocessing import extract_coloured_events

# paths
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

raw_files = os.listdir(raw_data_path)
raw_files.sort()
# %%update 202050626: add tb short and long
wb_70_1 = openpyxl.load_workbook(f'{raw_data_path}/070.1_Growth curves.xlsx', data_only=True)
ws_70_1 = wb_70_1['Tumour size in mm']
hex_red_70_1 = 'FFFF0000'  # red
# extract the coloured events from the first worksheet
df_70_1_meta, colours_70_1 = extract_coloured_events(
    ws=ws_70_1,
    events_list=['Sacrifice'],
    colour_list=[hex_red_70_1],  # red
    colour_type_list=['font']
)
# we have some double entries, probably due to hidden formatting. remove all entries where the value is nan
df_70_1_meta = df_70_1_meta.dropna(subset=['value'])
# now drop all where the mouse_ID and the event are duplicates. Keep 1st occurrence (due to summary table at the bottom)
df_70_1_meta = df_70_1_meta.drop_duplicates(subset=['mouse_ID', 'event'], keep='first')
#%% same for 70.2
wb_70_2 = openpyxl.load_workbook(f'{raw_data_path}/070.2_Growth curves.xlsx', data_only=True)
ws_70_2 = wb_70_2['Tumour size in mm']
df_70_2_meta, colours_70_2 = extract_coloured_events(
    ws=ws_70_2,
    events_list=['Sacrifice'],
    colour_list=[hex_red_70_1],  # red
    colour_type_list=['fill']
)
df_70_2_meta = df_70_2_meta.dropna(subset=['value'])
df_70_2_meta = df_70_2_meta.drop_duplicates(subset=['mouse_ID', 'event'], keep='first')
# %%open the worksheet 76 (cyclo)
wb_76_1 = openpyxl.load_workbook(f'{raw_data_path}/076.1_Groth curves.xlsx', data_only=True)
ws_76_1 = wb_76_1['Tumour size in mm']
hex_blue_76_1 = 'FF0070C0'
# sacrifice is always the next date, so no need to extract.
df_76_1_meta, colours_76_1 = extract_coloured_events(  # continue here, some extra errors called but seems to work OK.
    ws=ws_76_1,
    events_list=['Cyclo'],
    colour_list=[hex_blue_76_1],         # blue and red
    colour_type_list=['font']
)
# we have some double entries, probably due to hidden formatting. remove all entries where the value is nan
df_76_1_meta = df_76_1_meta.dropna(subset=['value'])
#%% 2nd workbook
hex_blue_76_2 = 'FF00C0F0'
# FF0070C0 blue
# FFFF0000 red
# FF00B0F0 other blue
wb_76_2 = openpyxl.load_workbook(f'{raw_data_path}/076.2_Growth curves.xlsx', data_only=True)
ws_76_2 = wb_76_2['Tumour size in mm']
df_76_2_meta, colours_76_2 = extract_coloured_events(
    ws=ws_76_2,
    events_list=['Cyclo'],
    colour_list=[['FF0070C0', 'FF00B0F0']],
    colour_type_list=['font']
)
df_76_2_meta = df_76_2_meta.dropna(subset=['value'])
#%% 3rd workbook
wb_76_3 = openpyxl.load_workbook(f'{raw_data_path}/076.3_Growth curves.xlsx', data_only=True)
ws_76_3 = wb_76_3['Tumour size in mm']
df_76_3_meta, colours_76_3 = extract_coloured_events(
    ws=ws_76_3,
    events_list=['Cyclo'],
    colour_list=[['FF00B0F0', 'FF0070C0']],
    colour_type_list=['fill']
)
# drop from index 6 (duplicates. no measurement on chemo day is correct)
df_76_3_meta = df_76_3_meta.iloc[:6]
#%% 79.1
wb_79_1 = openpyxl.load_workbook(f'{raw_data_path}/079.1_Growth curves.xlsx', data_only=True)
ws_79_1 = wb_79_1['Tumour size in mm']
df_79_1_meta, colours_79_1 = extract_coloured_events(
    ws=ws_79_1,
    events_list=['Sacrifice'],
    colour_list=['FFFF0000'],  # red
    colour_type_list=['font'],
)
df_79_1_meta = df_79_1_meta.dropna(subset=['value'])
df_79_1_meta = df_79_1_meta.drop_duplicates(subset=['mouse_ID', 'event'], keep='first')
# drop mice BAL-3970 and BAL-3971 (mice developed preputial gland abscesses and were sacrificed)
df_79_1_meta = df_79_1_meta[~df_79_1_meta['mouse_ID'].isin(['BAL-3970', 'BAL-3971'])]
#%% 79.2
wb_79_2 = openpyxl.load_workbook(f'{raw_data_path}/079.2_Growth curves.xlsx', data_only=True)
ws_79_2 = wb_79_2['Tumour size in mm']
df_79_2_meta, colours_79_2 = extract_coloured_events(
    ws=ws_79_2,
    events_list=['Sacrifice'],
    colour_list=[['FFFF0000', 'FFFF0066']],  # red, magenta (keep magenta mouse 2220, swap to tb3 in 2_ script.)
    colour_type_list=['fill'],
)
df_79_2_meta = df_79_2_meta.dropna(subset=['value'])
df_79_2_meta = df_79_2_meta.drop_duplicates(subset=['mouse_ID', 'event'], keep='first')
#%% 4th workbook
wb_86_1 = openpyxl.load_workbook(f'{raw_data_path}/086.1_Growth curves.xlsx', data_only=True)
ws_86_1 = wb_86_1['Tumour size in mm']
df_86_1_meta, colours_86_1 = extract_coloured_events(
    ws=ws_86_1,
    events_list=['Cyclo', 'Pmels+Virus', 'Sacrifice'],
    colour_list=[['FF00B0F0', 'FF0070C0'], ['FFA66BD3', 'FF7030A0'], ['FFFF0000']],
    colour_type_list=['fill', 'fill', 'fill']
)
# drop all rows with a 'value cell' entry that contains a row equal or larger than 100
df_86_1_meta = df_86_1_meta[df_86_1_meta['value cell'].str.extract(r'(\d+)').astype(int)[0] <= 100]
# now drop all where the mouse_ID and the event are duplicates. Keep 1st occurrence
# df_86_1_meta = df_86_1_meta.drop_duplicates(subset=['mouse_ID', 'event'], keep='first')
#%% 5th workbook
wb_87_1 = openpyxl.load_workbook(f'{raw_data_path}/087.1_Growth curves.xlsx', data_only=True)
ws_87_1 = wb_87_1['Tumour size in mm']
df_87_1_meta, colours_87_1 = extract_coloured_events(
    ws=ws_87_1,
    events_list=['Sacrifice', 'Cyclo', 'Pmels+Virus', 'CpG/PolyIC'],
    colour_list=[['FFFF0000', 'FFFF0066'], ['FF00B0F0'], ['FF9751CB'], ['FF00B050']],
    colour_type_list=['fill', 'fill', 'fill', 'fill']
)
df_87_1_meta_cut = df_87_1_meta.drop_duplicates(subset=['mouse_ID', 'event', 'experiment date'], keep='first')
df_87_1_meta_cut = df_87_1_meta_cut[df_87_1_meta_cut['value cell'].str.extract(r'(\d+)').astype(int)[0] <= 271]
# %%print any mouse_ID's that have less than 6 events
mouse_ids = df_87_1_meta_cut['mouse_ID'].value_counts()
df_mouse_2282 = df_87_1_meta_cut[df_87_1_meta_cut['mouse_ID'] == 'HOE-2282']
# has 1 CpG/polyIC treatment missing, is also the case in the raw data. I asked Michelle for verification.
# conclusion: yes, this mouse only had 2 treatments, and developed his tumor late as well. TG and other
# behaviors are similar to mice with 3 treatments, so we will include it for now.
#%% 6th workbook
print(f'6th book: {raw_files[7]}, ACT d7')
wb_96_1 = openpyxl.load_workbook(f'{raw_data_path}/096.1_Growth curves.xlsx', data_only=True)
ws_96_1 = wb_96_1['Tumour size in mm']
df_96_1_meta, colours_96_1 = extract_coloured_events(
    ws=ws_96_1,
    events_list=['Sacrifice', 'Cyclo', 'Pmels+Virus', 'CpG/PolyIC'],
    colour_list=[['FFFF0000'], ['FF4FD1FF', 'FF00B0F0', 'FF0070C0'], ['FF7030A0', 'FF9751CB'], ['FF00B050']],
    colour_type_list=['fill', 'fill', 'fill', 'fill']
)
df_96_1_meta['mouse_ID'].value_counts()  # correct
#%% 7th workbook
print(f'7th book: {raw_files[8]}, ACT d14')
wb_97_1 = openpyxl.load_workbook(f'{raw_data_path}/097.1_Growth curves.xlsx', data_only=True)
ws_97_1 = wb_97_1['Tumour size in mm']
df_97_1_meta, colours_97_1 = extract_coloured_events(
    ws=ws_97_1,
    events_list=['Sacrifice', 'Cyclo', 'Pmels+Virus', 'CpG/PolyIC', 'CpG/PolyIC'],  # extra here because font is green
    colour_list=[['FFFF0000'], ['FF00B0F0', 'FF0070C0'], ['FF7030A0', 'FF9751CB'], ['FF00B050'], ['FF00B050']],
    colour_type_list=['fill', 'fill', 'fill', 'fill', 'font']
)
df_97_1_meta['mouse_ID'].value_counts()
#%% 8th workbook
print(f'8th book: {raw_files[9]}, ACT d14, mixed')
wb_102_1 = openpyxl.load_workbook(f'{raw_data_path}/102.1_Growth curves.xlsx', data_only=True)
ws_102_1 = wb_102_1['Tumour size in mm']
df_102_1_meta, colours_102_1 = extract_coloured_events(
    ws=ws_102_1,
    events_list=['Sacrifice', 'Cyclo', 'Pmels+Virus', 'CpG/PolyIC', "No CpG/polyIC"],
    colour_list=[['FFFF0000', 'FF92D050'], ['FFA162D0', 'FFA86ED4'], ['FF00B0F0', 'FF0070C0'], ['FFFFFF00'], ['FFCC9B00']],
    colour_type_list=['fill', 'fill', 'fill', 'fill', 'font']
)
#%%
df_102_1_meta['mouse_ID'].value_counts()
#%%
df_mouse_2474 = df_102_1_meta[df_102_1_meta['mouse_ID'] == 'HOE-2474']
# drop any rows with a 'date' value later than the date of the 'Sacrifice' event for that mouse
# %%get the date of the sacrifice event for each mouse # 3/20 didn't develop tumors, sacrifice date not recorded.
sacrifice_dates = df_102_1_meta[df_102_1_meta['event'] == 'Sacrifice'].set_index('mouse_ID')['experiment date']
# drop any rows with a 'date' value later than the date of the 'Sacrifice' event for that mouse
df_102_1_meta_cut = df_102_1_meta.join(sacrifice_dates, on='mouse_ID', rsuffix='_sacrifice')
df_102_1_meta_cut = df_102_1_meta_cut[df_102_1_meta_cut['experiment date'] <= df_102_1_meta_cut['experiment date_sacrifice']]
df_102_1_meta_cut['mouse_ID'].value_counts()
#%%
df_mouse_2480 = df_102_1_meta_cut[df_102_1_meta_cut['mouse_ID'] == 'HOE-2480']  # hidden formatting no cgp.
# drop row containing value cell Q93 (hidden formatting)
df_102_1_meta_cut = df_102_1_meta_cut[df_102_1_meta_cut['value cell'] != 'Q93']
df_mouse_2481 = df_102_1_meta_cut[df_102_1_meta_cut['mouse_ID'] == 'HOE-2481'].iloc[:,:7]  # hidden formatting no cgp.
df_102_1_meta_cut = df_102_1_meta_cut[df_102_1_meta_cut['value cell'] != 'O144']
df_mouse_2470 = df_102_1_meta_cut[df_102_1_meta_cut['mouse_ID'] == 'HOE-2470'].iloc[:,:7]  # hidden formatting no cgp.
df_102_1_meta_cut = df_102_1_meta_cut[df_102_1_meta_cut['value cell'] != 'O133']
df_mouse_2472 = df_102_1_meta_cut[df_102_1_meta_cut['mouse_ID'] == 'HOE-2472'].iloc[:,:7]  # OK
# check that mice don't contain the event 'No CpG/polyIC' and the event 'CpG/PolyIC', print event types per mouse
events_per_mouse = df_102_1_meta_cut.groupby('mouse_ID')['event'].unique()
conflict_mice = events_per_mouse[
    events_per_mouse.apply(lambda events: 'CpG/PolyIC' in events and 'No CpG/polyIC' in events)
] # no conflicts
df_102_1_meta_cut = df_102_1_meta_cut.iloc[:, :-1]  # drop the last column, unneeded.
#%% 9th workbook
print(f'9th book: {raw_files[11]}, ACT d7, no CpG')
wb_103_1 = openpyxl.load_workbook(f'{raw_data_path}/103.1_Growth curves.xlsx', data_only=True)
ws_103_1 = wb_103_1['Tumour size in mm']
df_103_1_meta, colours_103_1 = extract_coloured_events(
    ws=ws_103_1,
    events_list=['Sacrifice', 'Cyclo', 'Pmels+Virus'],
    colour_list=[['FFFF0000'], ['FFA162D0'], ['FF00B0F0']],
    colour_type_list=['fill', 'fill', 'fill']
)
df_103_1_meta['mouse_ID'].value_counts() # correct, mouse 2492 never developed a tumour, so 6 mice in cohort

# %%save the dataframes to file
df_70_1_meta.to_csv(f'{processed_data_path}/metadata/per_experiment/tumor_growth_curves_070_1_meta_{modifier}.csv')
df_70_2_meta.to_csv(f'{processed_data_path}/metadata/per_experiment/tumor_growth_curves_070_2_meta_{modifier}.csv')
df_76_1_meta.to_csv(f'{processed_data_path}/metadata/per_experiment/tumor_growth_curves_076_1_meta_{modifier}.csv')
df_76_2_meta.to_csv(f'{processed_data_path}/metadata/per_experiment/tumor_growth_curves_076_2_meta_{modifier}.csv')
df_76_3_meta.to_csv(f'{processed_data_path}/metadata/per_experiment/tumor_growth_curves_076_3_meta_{modifier}.csv')
df_79_1_meta.to_csv(f'{processed_data_path}/metadata/per_experiment/tumor_growth_curves_079_1_meta_{modifier}.csv')
df_79_2_meta.to_csv(f'{processed_data_path}/metadata/per_experiment/tumor_growth_curves_079_2_meta_{modifier}.csv')
df_86_1_meta.to_csv(f'{processed_data_path}/metadata/per_experiment/tumor_growth_curves_086_1_meta_{modifier}.csv')
df_87_1_meta_cut.to_csv(f'{processed_data_path}/metadata/per_experiment/tumor_growth_curves_087_1_meta_{modifier}.csv')
df_96_1_meta.to_csv(f'{processed_data_path}/metadata/per_experiment/tumor_growth_curves_096_1_meta_{modifier}.csv')
df_97_1_meta.to_csv(f'{processed_data_path}/metadata/per_experiment/tumor_growth_curves_097_1_meta_{modifier}.csv')
df_102_1_meta_cut.to_csv(f'{processed_data_path}/metadata/per_experiment/tumor_growth_curves_102_1_meta_{modifier}.csv')
df_103_1_meta.to_csv(f'{processed_data_path}/metadata/per_experiment/tumor_growth_curves_103_1_meta_{modifier}.csv')
