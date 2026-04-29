#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
# '''=================================================
# @Project: ACT_melanoma_neutrophil (data: act_neutrophil_data_repository)
# @Author : Gemma van der Voort
# @Desc   : visualise the tumor growth data. one plot per experiment with the individual mice trajectories,
# @Desc updated: Individual tumor growth trajectories and summary exports.
# and a generic plot with the mean and stdev per experiment.
# update 21.03.2025, for v2: offset at first tumour growth
# update 25.03.2025, for v3: offset at first treatment: cyclo. Also annotate other treatments.
# update 26.06.2025, for v5: add the individual trajectories for the two TB cohorts (70, 79)
# update 21.01.2026: export mice no growth/all mice numbers for dataset overview table. update linestyles and font size.
# '''=================================================
import os
import sys
from pathlib import Path

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.lines import Line2D
from matplotlib.widgets import EllipseSelector

PROJECT_ROOT = Path(__file__).resolve()
while PROJECT_ROOT.name != "ACT_melanoma_neutrophil" and PROJECT_ROOT.parent != PROJECT_ROOT:
    PROJECT_ROOT = PROJECT_ROOT.parent
if PROJECT_ROOT.name != "ACT_melanoma_neutrophil":
    raise RuntimeError("ACT_melanoma_neutrophil root not found")
sys.path.append(str(PROJECT_ROOT / "helper_files"))

from paths_parameters import titles_tumour_growth, colours_tumour_growth, data_repo_path

# %%style sheet (we are in a subfolder)
plt.style.use(str(PROJECT_ROOT / "helper_files" / "journal_style.mplstyle"))

processed_data_path = f'{data_repo_path}/processed/tumor_growth/per_treatment'
processed_data_path_meta = f'{data_repo_path}/processed/tumor_growth/metadata/per_treatment'
fig_dir = f'{data_repo_path}/figures/tumor_growth_curves/per_condition'
incl_no_CpG = True
if incl_no_CpG:
    modifier = 'incl_no_CpG'
else:
    modifier = ''
modifier += 'incl_tb3_tb8_newstyle'
file_modifier = 'data_corrections_sep2025'
files_growth = os.listdir(processed_data_path)
files_meta = os.listdir(processed_data_path_meta)
df_cyclo = pd.read_csv(f'{processed_data_path}/Cyclo_{file_modifier}.csv', index_col=0)
df_cyclo_meta = pd.read_csv(f'{processed_data_path_meta}/Cyclo_{file_modifier}.csv', index_col=0)
df_ACT_d3 = pd.read_csv(f'{processed_data_path}/ACT_day_3_{file_modifier}.csv', index_col=0)
df_ACT_d3_meta = pd.read_csv(f'{processed_data_path_meta}/ACT_day_3_{file_modifier}.csv', index_col=0)
df_ACT_d7 = pd.read_csv(f'{processed_data_path}/ACT_day_7_{file_modifier}.csv', index_col=0)
df_ACT_d7_meta = pd.read_csv(f'{processed_data_path_meta}/ACT_day_7_{file_modifier}.csv', index_col=0)
df_ACT_d14 = pd.read_csv(f'{processed_data_path}/ACT_day_14_{file_modifier}.csv', index_col=0)
df_ACT_d14_meta = pd.read_csv(f'{processed_data_path_meta}/ACT_day_14_{file_modifier}.csv', index_col=0)
df_relapse = pd.read_csv(f'{processed_data_path}/Relapse_{file_modifier}.csv', index_col=0)
df_relapse_meta = pd.read_csv(f'{processed_data_path_meta}/Relapse_{file_modifier}.csv', index_col=0)
df_ACT_d7_no_CpG = pd.read_csv(f'{processed_data_path}/ACT_day_7_no_CpG_{file_modifier}.csv', index_col=0)
df_ACT_d7_no_CpG_meta = pd.read_csv(f'{processed_data_path_meta}/ACT_day_7_no_CpG_{file_modifier}.csv', index_col=0)
df_ACT_d14_no_CpG = pd.read_csv(f'{processed_data_path}/ACT_day_14_no_CpG_{file_modifier}.csv', index_col=0)
df_ACT_d14_no_CpG_meta = pd.read_csv(f'{processed_data_path_meta}/ACT_day_14_no_CpG_{file_modifier}.csv', index_col=0)
#%% tb files:
df_tb3 = pd.read_csv(f'{processed_data_path}/Tumor_bearing_{file_modifier}.csv', index_col=0)
df_tb8 = pd.read_csv(f'{processed_data_path}/Tumor_bearing_8-10mm_{file_modifier}.csv', index_col=0)
df_tb3_meta = pd.read_csv(f'{processed_data_path_meta}/Tumor_bearing_{file_modifier}.csv', index_col=0)
pseudo_cyclo_day_tb8 = 16
df_tb8_meta = pd.read_csv(f'{processed_data_path_meta}/Tumor_bearing_8-10mm_incl_pseudo_cyclo'
                          f'_{pseudo_cyclo_day_tb8}_{file_modifier}.csv', index_col=0)

#%%  We will plot individual trajectories, adding event information encoded in df**_meta.
# The plot begins per mouse on the day of the 'cyclo' event. Pmel+Virus is the next event, indicated with a round dot
# on the timeline. CpG/PolyIC events (3 per mouse) are indicated with a square on the line.

titles_tumour_growth_incl_tb = titles_tumour_growth[1:] + ["Tumor bearing", "Untreated"] # "Tumor bearing 8-10 mm"
dfs = [df_cyclo, df_ACT_d3, df_ACT_d7, df_ACT_d14, df_relapse, df_ACT_d7_no_CpG, df_ACT_d14_no_CpG, df_tb3, df_tb8]
meta_dfs = [df_cyclo_meta, df_ACT_d3_meta, df_ACT_d7_meta, df_ACT_d14_meta, df_relapse_meta, df_ACT_d7_no_CpG_meta,
            df_ACT_d14_no_CpG_meta, df_tb3_meta, df_tb8_meta]
dfs_realigned = []
meta_dfs_realigned = []
sacrifice_dates = [0,3,7,14,'Relapse',7,14, -1, 14]
df_no_growth_summary = pd.DataFrame(columns=['Treatment', 'Mice without tumor growth', 'All mice'])
for df_index, (data_df, meta_df) in enumerate(zip(dfs, meta_dfs)):
    print(f'processing: {titles_tumour_growth_incl_tb[df_index]}')
    # Filter out mice that remained below an area of 9 (3x3 mm) for the entire experiment
    n_mice = data_df.shape[0]
    mice_no_growth = data_df[data_df.max(axis=1) < 9].index
    print(f'{titles_tumour_growth_incl_tb[df_index]}: {mice_no_growth} mice with no growth')
    n_mice_zero = len(mice_no_growth)
    data_df = data_df.drop(index=mice_no_growth)
    max_col = data_df.columns[data_df.max() > 0][-1]  # last time-point with live mice
    # create empty columns with the timepoints where no measurements were taken.
    # remove columns after max_col for removing redundant timepoints without live mice:
    data_df = data_df.loc[:, :max_col]
    # add columns for the missing dates. Some days have no data, so we need to add them to the df so we can shift the
    # values up by the first non-zero value per mouse
    cols_as_int = data_df.columns.astype(int)
    missing_days = set(range(cols_as_int.min(), cols_as_int.max() + 1)) - set(cols_as_int)
    for day in missing_days:
        data_df[str(day)] = np.nan
    data_df.columns = data_df.columns.astype(int)
    # reorder the columns by day order
    data_df = data_df.reindex(sorted(data_df.columns), axis=1)
    # Step 1: Build realigned data
    realigned_data = []
    realigned_meta = []
    cyclo_days = {}

    for mouse in data_df.index:
        print(mouse)
        # Get cyclo day for this mouse
        if df_index == 8:  # TB8 cohort
            event = 'Pseudo-cyclo'
        elif df_index == 7:  # TB3 cohort
            event = 'Sacrifice'  # these mice were sacrificed at the cyclo date
        else:
            event = 'Cyclo'
        cyclo_event = meta_df[(meta_df['mouse_ID'] == mouse) & (meta_df['event'] == event)]
        if cyclo_event.empty:
            print(f'Skipping {mouse} with no cyclo event')
            continue

        cyclo_day = (cyclo_event['experiment date'].values[0])
        cyclo_days[mouse] = cyclo_day

        # Get mouse series, reindex to days since cyclo
        mouse_series = data_df.loc[mouse]  #.dropna()
        aligned_series = mouse_series.copy()
        aligned_series.index = aligned_series.index.astype(int)
        aligned_series.index = aligned_series.index - (cyclo_day +1)  # shift so pmel is day 0. doesn't shift data but changes
        # index, so lowest is 0-pmel day (e.g. -16)
        # since now aligned to treatment start, there is not data after index 14 (sacrifice).
        # assert all indices >14 contain only nans
        if df_index not in [4,8]:  # skip for relapse and tb 8
            assert aligned_series[aligned_series.index > sacrifice_dates[df_index]].isna().all()
            aligned_series = aligned_series[aligned_series.index <= sacrifice_dates[df_index]]  # remove days after 14

        realigned_data.append(pd.DataFrame({
            'mouse_ID': mouse,
            'day': aligned_series.index,
            'value': aligned_series.values
        }))

        # now same for metadata: subtract the pmel day from each experiment date for the mouse
        meta_mouse = meta_df[meta_df['mouse_ID'] == mouse].copy()
        meta_mouse['aligned_day'] = meta_mouse['experiment date'] - (cyclo_day +1) # pmel at 0
        realigned_meta.append(meta_mouse)

    df_aligned = pd.concat(realigned_data, ignore_index=True)
    meta_aligned = pd.concat(realigned_meta, ignore_index=True)
    dfs_realigned += [df_aligned]
    meta_dfs_realigned += [meta_aligned]

    #  Plot
    fig, ax = plt.subplots(figsize=(10, 6))

    # Plot mouse lines
    sns.lineplot(data=df_aligned, x='day', y='value', hue='mouse_ID', ax=ax, palette='tab20', dashes=False)

    # Add vertical lines with different dotted styles
    if df_index not in [7, 8]:  # skip for TB3 and TB8 cohorts
        ax.axvline(x=-1, color='gray', linestyle=(0, (1, 5)), linewidth=1.5, label='Cy')  #
    if df_index not in [0,7,8]:  # skip for cyclo and tb
        ax.axvline(x=0, color='gray', linestyle=(0, (5, 5)), linewidth=1.5, label='Pmels+Virus')
    if df_index == 2:  # act d7 has 2 CpG events
        for CpGday in [3, 6]:
            ax.axvline(x=CpGday, color='gray', linestyle=(0, (3, 5, 1, 5)), linewidth=1.5, label='CpG/PolyIC' if CpGday == 3 else None)
    elif df_index not in [0, 1, 5, 6, 7, 8]:  # skip for cyclo, act d3, no CpGs
        for CpGday in [3, 6, 9]:
            ax.axvline(x=CpGday, color='gray', linestyle=(0, (3, 5, 1, 5)), linewidth=1.5, label='CpG/PolyIC' if CpGday == 3 else None)
    else:
        print(f'{titles_tumour_growth_incl_tb[df_index]} has no CpG events')
    if df_index not in [4, 8]:  # skip for relapse, TB8
        ax.axvline(x=sacrifice_dates[df_index], color='gray', linestyle='-', linewidth=1.5, label='Sacrifice')  # solid

    if df_index == 8:  # TB8 cohort,
        ax.set_xlabel('Day (0=est. treatment start)')
    else:
        ax.set_xlabel("Days (0=treatment start)")
    ax.set_ylabel("Tumor area [mm$^2$]")
    sns.despine()
    # legend outside of plot
    # ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), title='Mouse ID')
    # remove legend
    ax.legend().remove()
    if df_index not in [4]:  # skip for relapse
        ax.set_xlim(-5, (sacrifice_dates[df_index]+1))
        # set specific x ticks, aligned between plots for comparison
        x_ticks = [-4, -1, 0, 3, 6, 9, 14]
        ax.set_xticks(x_ticks)
        # Add text annotation for mice that remained at zero to the top left
        # ax.text(x=-4.8, y=int(max(data_df.max()) * 0.9),
        #         s=f'Mice without\ntumour growth:\n{n_mice_zero}/{n_mice}', color='black', verticalalignment='bottom')
    else:
        ax.set_xlim(-5, 60)
        x_ticks = [-4, -1, 0, 3, 6, 9, 14, 20,30,40,50]
        ax.set_xticks(x_ticks)
        # Add text annotation for mice that remained at zero to the top left
        # ax.text(x=9.2, y=int(max(data_df.max()) * 0.9),
        #         s=f'Mice without\ntumour growth:\n{n_mice_zero}/{n_mice}', color='black', verticalalignment='bottom')

    plt.tight_layout()
    plt.savefig(f'{fig_dir}/{titles_tumour_growth_incl_tb[df_index]}_tumor_growth_cyclo_offset_{file_modifier}_'
                f'{modifier}.pdf')
    plt.close()
    # add to summary df for table
    df_no_growth_summary = pd.concat([df_no_growth_summary,
                                     pd.DataFrame({'Treatment': [titles_tumour_growth_incl_tb[df_index]],
                                                   'Mice without tumor growth': [n_mice_zero],
                                                   'All mice': [n_mice]})],
                                    ignore_index=True)
# save summary df
df_no_growth_summary.to_csv(f'{processed_data_path}/tumor_growth_mice_no_growth_summary_{modifier}_{file_modifier}.csv', index=False)
