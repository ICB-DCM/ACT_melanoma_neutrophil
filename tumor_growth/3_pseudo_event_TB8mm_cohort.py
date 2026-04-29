#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
# '''=================================================
# @Project: ACT_melanoma_neutrophil (data: act_neutrophil_data_repository)
# @Author : Gemma van der Voort
# @Desc   : calculate the time when TB8 mice reach a tumor area of 22mm2 (mean/median/mode of TB3 mice on sacrifice day)
# @Desc updated: Add pseudo-cyclo events to TB8 cohort metadata based on TB3 distribution.
# and set it as a pseudo-cyclo day in the metadata file for TB8 mice, so we can compare them to the other cohorts.
# update: 22 is too late, since that means we cut out the first half of the sacrifice aerea histogram.
# We will use 16 instead, which is the 4mm2 tumour size, middle between the 3-5mm2 target. We will still miss the
# start of the curve (at 9) but at least don't underestimate too much.
# '''=================================================
import sys
from pathlib import Path

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.lines import Line2D

PROJECT_ROOT = Path(__file__).resolve()
while PROJECT_ROOT.name != "ACT_melanoma_neutrophil" and PROJECT_ROOT.parent != PROJECT_ROOT:
    PROJECT_ROOT = PROJECT_ROOT.parent
if PROJECT_ROOT.name != "ACT_melanoma_neutrophil":
    raise RuntimeError("ACT_melanoma_neutrophil root not found")
sys.path.append(str(PROJECT_ROOT / "helper_files"))

from paths_parameters import titles_tumour_growth, colours_tumour_growth, data_repo_path

plt.style.use(str(PROJECT_ROOT / "helper_files" / "journal_style.mplstyle"))

processed_data_path = f'{data_repo_path}/processed/tumor_growth/per_treatment'
processed_data_path_meta = f'{data_repo_path}/processed/tumor_growth/metadata/per_treatment'
fig_dir = f'{data_repo_path}/figures/tumor_growth_curves/per_condition'
modifier = 'data_corrections_sep2025'
#%% no treatment dfs
df_tb3 = pd.read_csv(f'{processed_data_path}/Tumor_bearing_{modifier}.csv', index_col=0)
df_tb8 = pd.read_csv(f'{processed_data_path}/Tumor_bearing_8-10mm_{modifier}.csv', index_col=0)
df_tb3_meta = pd.read_csv(f'{processed_data_path_meta}/Tumor_bearing_{modifier}.csv', index_col=0)
df_tb8_meta = pd.read_csv(f'{processed_data_path_meta}/Tumor_bearing_8-10mm_{modifier}.csv', index_col=0)
#%% individual plots for the non-treated (TB3 and TB8) cohorts
# we need to make a pseudo-cyclo day for tb8, since we align them by this. For this, let's find out at what size the TB3
# mice were sacrificed.
# get the area on sacrifice day (see meta df for sacrifice date) per mouse
mouse_sacrifice_dict_tb3 = {}

for mouse in df_tb3.index:
    print(mouse)
    sacrifice_event =  df_tb3_meta[(df_tb3_meta['mouse_ID'] == mouse) & (df_tb3_meta['event'] == 'Sacrifice')]
    if sacrifice_event.empty:
        print(f'Skipping {mouse} with no sacrifice event')
        continue
    sacrifice_day = sacrifice_event['experiment date'].values[0]
    mouse_sacrifice_dict_tb3[mouse] = sacrifice_day
# %%now get the area for that mouse on that day (column)
final_areas_tb3 = []
for mouse, sacrifice_day in mouse_sacrifice_dict_tb3.items():
    if str(sacrifice_day) in df_tb3.columns:
        final_areas_tb3.append(df_tb3.loc[mouse, str(sacrifice_day)])
    else:
        print(f'Skipping {mouse} with no area on sacrifice day {sacrifice_day}')
# drop the final areas that are 0
final_areas_tb3 = [area for area in final_areas_tb3 if area > 0]
# plot these in a histogram and print the mean and median (plus their values) in the plot
plt.figure(figsize=(8, 6))
sns.histplot(final_areas_tb3, bins=100, kde=True)
plt.xlabel("Tumour area on sacrifice day (mm$^2$)")
plt.ylabel("Number of mice")
plt.title("Distribution of tumour area on sacrifice day for TB3 mice")
plt.axvline(np.mean(final_areas_tb3), color='red', linestyle='--', label='Mean: {:.2f}'.format(np.mean(final_areas_tb3)))
plt.axvline(np.median(final_areas_tb3), color='blue', linestyle='--', label='Median: {:.2f}'.format(np.median(final_areas_tb3)))
plt.legend()
sns.despine()
plt.savefig(f'{fig_dir}/TB3_tumor_area_on_sacrifice_day.pdf')
plt.close()
# %%mean, median and mode almost coincide; we will use area > 22 as our pseudo-cyclo day for TB8
# pseudo_cyclo_day_tb8 = 22
# correction: the above means all pseudo times are after day 22, disregarding the first half of the distribution and
# thus skewing the pseudo time to be later on average that the actual time. Instead, we will go with the
# experimentally determined middle point: 4 mm, so an area of 16 mm^2.
pseudo_cyclo_day_tb8 = 16
# per mouse, find the day when the area first exceeds pseudo_cyclo_day_tb8
mouse_pseudo_cyclo_day_dict_tb8 = {}
for mouse in df_tb8.index:
    # find the first day when the area exceeds 16
    area_exceeds_22 = df_tb8.loc[mouse][df_tb8.loc[mouse] > pseudo_cyclo_day_tb8]  # lists the days when the area exceeds 22 (index)
    if not area_exceeds_22.empty:
        first_day = area_exceeds_22.index[0]
        mouse_pseudo_cyclo_day_dict_tb8[mouse] = first_day
    else:
        print(f'Skipping {mouse} with no area exceeding {pseudo_cyclo_day_tb8}')
        # done 4x, these mice didn't develop a tumour.
# add this as metadata to the df_tb8_meta. in column 'event' add 'Pseudo-cyclo' and 'experiment date' the first day
for mouse, pseudo_cyclo_day in mouse_pseudo_cyclo_day_dict_tb8.items():
    df_tb8_meta = pd.concat([df_tb8_meta, pd.DataFrame({'mouse_ID': [mouse], 'event': ['Pseudo-cyclo'],
                                                        'experiment date': [int(pseudo_cyclo_day)]})], ignore_index=True)
# %%check in the meta df that the sacrifice and pseudo-cyclo days are not the same
for mouse in df_tb8_meta['mouse_ID'].unique():
    mouse_rows = df_tb8_meta[df_tb8_meta['mouse_ID'] == mouse]
    # if mouse contains pseudo-cyclo, get the pseudo-cyclo day and sacrifice day
    pseudo_cyclo_rows = mouse_rows[mouse_rows['event'] == 'Pseudo-cyclo']
    sacrifice_rows = mouse_rows[mouse_rows['event'] == 'Sacrifice']
    if pseudo_cyclo_rows.empty:
        print(f'Skipping {mouse} with no pseudo-cyclo event.')
        continue  # skip if no pseudo-cyclo
    pseudo_cyclo_day = pseudo_cyclo_rows['experiment date'].values[0]
    sacrifice_day = sacrifice_rows['experiment date'].values[0]
    if pseudo_cyclo_day == sacrifice_day:
        print(f'Warning: Pseudo-cyclo day and sacrifice day are the same for mouse {mouse}.')
# save the updated meta df
df_tb8_meta.to_csv(f'{processed_data_path_meta}/Tumor_bearing_8-10mm_incl_pseudo_cyclo_{pseudo_cyclo_day_tb8}_{modifier}.csv', index=False)
