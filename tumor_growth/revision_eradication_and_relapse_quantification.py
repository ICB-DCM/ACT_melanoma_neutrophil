#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
# '''=================================================
# @Project -> File :ACT_melanoma_neutrophil -> revision_eradication_and_relapse_quantification.py
# @Author : Gemma van der Voort
# @Time   : 3/30/26 6:22 PM
# @Desc   : quantification of ACT efficacy in terms of tumor eradication, reduction and relapse,
# based on aligned tumor growth curves.
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

# %%style sheet
plt.style.use(str(PROJECT_ROOT / "helper_files" / "journal_style.mplstyle"))

processed_data_path = f'{data_repo_path}/processed/tumor_growth/per_treatment'
processed_data_path_meta = f'{data_repo_path}/processed/tumor_growth/metadata/per_treatment'
fig_dir = f'{data_repo_path}/figures/tumor_growth_curves/revision'
#%% load data
file_modifier = 'data_corrections_sep2025'
# files_growth = os.listdir(processed_data_path)
# files_meta = os.listdir(processed_data_path_meta)
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

#%% same as 4a and 4b: align data first. add extra checks here for tumor eradication, reduction and relapse.
titles_tumour_growth_incl_tb = titles_tumour_growth[1:] + ["Tumor bearing", "Untreated"]
dfs = [df_cyclo, df_ACT_d3, df_ACT_d7, df_ACT_d14, df_relapse, df_ACT_d7_no_CpG, df_ACT_d14_no_CpG, df_tb3, df_tb8]
meta_dfs = [df_cyclo_meta, df_ACT_d3_meta, df_ACT_d7_meta, df_ACT_d14_meta, df_relapse_meta, df_ACT_d7_no_CpG_meta,
            df_ACT_d14_no_CpG_meta, df_tb3_meta, df_tb8_meta]
dfs_realigned = []
meta_dfs_realigned = []
sacrifice_dates = [0,3,7,14,'Relapse',7,14, -1, 14]
df_updated_growth_summary = pd.DataFrame(
    columns=[
        'Treatment',
        'Total mice',
        'Mice without tumor growth',
        'Mice with tumor eradication',
        'Mice with tumor reduction',
        'Mice with relapse',
        'Relapse time frame'
    ]
)

dfs_realigned = []
meta_dfs_realigned = []

for df_index, (data_df, meta_df) in enumerate(zip(dfs, meta_dfs)):
    print(f'processing: {titles_tumour_growth_incl_tb[df_index]}')
    # testing: skip if not d14
    # if df_index != 4:
    #     continue

    # Filter out mice that remained below an area of 9 (3x3 mm) for the entire experiment
    n_mice = data_df.shape[0]
    mice_no_growth = data_df[data_df.max(axis=1) < 9].index
    print(f'{titles_tumour_growth_incl_tb[df_index]}: {list(mice_no_growth)} mice with no growth')
    n_mice_zero = len(mice_no_growth)

    data_df = data_df.drop(index=mice_no_growth)

    # keep only until last time-point with any live mice
    max_col = data_df.columns[data_df.max() > 0][-1]
    data_df = data_df.loc[:, :max_col]

    # add missing days so shifting/alignment is well-defined
    cols_as_int = data_df.columns.astype(int)
    missing_days = set(range(cols_as_int.min(), cols_as_int.max() + 1)) - set(cols_as_int)
    for day in missing_days:
        data_df[str(day)] = np.nan

    data_df.columns = data_df.columns.astype(int)
    data_df = data_df.reindex(sorted(data_df.columns), axis=1)

    # alignment outputs
    realigned_data = []
    realigned_meta = []
    cyclo_days = {}

    mice_eradication = []
    mice_reduction = []

    for mouse in data_df.index:
        print(mouse)
        # if mouse != 'HOE-2285':
        #     continue # check 1 mouse trajectory

        # Get alignment event for this mouse
        if df_index == 8:  # TB8 cohort
            event = 'Pseudo-cyclo'
        elif df_index == 7:  # TB3 cohort
            event = 'Sacrifice'
        else:
            event = 'Cyclo'

        cyclo_event = meta_df[
            (meta_df['mouse_ID'] == mouse) &
            (meta_df['event'] == event)
        ]

        if cyclo_event.empty:
            print(f'Skipping {mouse} with no cyclo event')
            continue

        cyclo_day = cyclo_event['experiment date'].values[0]
        cyclo_days[mouse] = cyclo_day

        # Align tumor series so treatment start is day 0
        mouse_series = data_df.loc[mouse].copy()
        mouse_series.index = mouse_series.index.astype(int)
        mouse_series.index = mouse_series.index - (cyclo_day + 1)

        # Since aligned to treatment start, trim after sacrifice where appropriate
        if df_index not in [4, 8]:  # skip for relapse and TB8
            assert mouse_series[mouse_series.index > sacrifice_dates[df_index]].isna().all()
            mouse_series = mouse_series[mouse_series.index <= sacrifice_dates[df_index]]

        # Save aligned tumor data
        realigned_data.append(
            pd.DataFrame({
                'mouse_ID': mouse,
                'day': mouse_series.index,
                'value': mouse_series.values
            })
        )

        # Align metadata too
        meta_mouse = meta_df[meta_df['mouse_ID'] == mouse].copy()
        meta_mouse['aligned_day'] = meta_mouse['experiment date'] - (cyclo_day + 1)
        realigned_meta.append(meta_mouse)

        # --------------------------------------------------
        # response classification on aligned series
        # --------------------------------------------------
        mouse_nonan = mouse_series.dropna()
        # print(mouse_nonan)

        if len(mouse_nonan) == 0:
            continue

        # Only consider mice that reached treatment-sized tumor at any point
        treatment_size_days = mouse_nonan[mouse_nonan > 9].index
        if len(treatment_size_days) == 0:
            continue

        # --------------------
        # tumor eradication
        # definition:
        # tumor has reached treatment size (9mm) and reduced to at any time 0mm for at
        # least 2 consecutive measurements after that
        # --------------------
        last_treatment_size_day = treatment_size_days[-1] # last day the tumor was above treatment size.
        post_treatment = mouse_series[mouse_series.index > last_treatment_size_day]
        post_treatment_nonan = post_treatment.dropna()

        if len(post_treatment_nonan) > 0 and (post_treatment_nonan == 0).any():
            # check for at least 2 consecutive zeros
            is_zero = post_treatment_nonan == 0
            has_consecutive_zeros = (
                (is_zero.astype(int) + is_zero.astype(int).shift(1, fill_value=0)) == 2
            ).any()

            if has_consecutive_zeros:
                print(f'{mouse} classified as eradication')
                mice_eradication.append(mouse)

        # --------------------
        # tumor reduction
        # definition:
        # take the maximum tumor size at aligned day <= 9 (day of final treatment)
        # then require >= 2 consecutive later measurements
        # that are each at least 10% lower than that maximum
        # later regrowth does not matter
        # --------------------
        minimal_reduction = 0.20 # frac of tumor size decrease required to count as reduction
        max_treatment_day = 14 #9
        up_to_d9 = mouse_series[mouse_series.index <= max_treatment_day].dropna()

        if len(up_to_d9) > 0:
            max_tumor_size = up_to_d9.max()

            # use the last occurrence of the max at or before day 9
            # (in case of two of the same peaks or a plateau, but doesn't happen)
            max_day = up_to_d9[up_to_d9 == max_tumor_size].index[-1]

            post_max = mouse_series[mouse_series.index > max_day].dropna()

            if len(post_max) >= 2: # we cannot assess tumor reduction if the mouse was sacrificed after the peak.
                reduction_threshold = max_tumor_size * (1 - minimal_reduction)
                is_reduced = post_max < reduction_threshold

                has_consecutive_reduction = (
                    (is_reduced.astype(int) + is_reduced.astype(int).shift(1, fill_value=0)) == 2
                ).any()

                if has_consecutive_reduction:
                    mice_reduction.append(mouse)

    df_aligned = pd.concat(realigned_data, ignore_index=True)
    meta_aligned = pd.concat(realigned_meta, ignore_index=True)

    dfs_realigned.append(df_aligned)
    meta_dfs_realigned.append(meta_aligned)

    print(f'{titles_tumour_growth_incl_tb[df_index]}: {mice_eradication} mice with tumor eradication')
    print(f'{titles_tumour_growth_incl_tb[df_index]}: {mice_reduction} mice with tumor reduction')

    n_mice_eradication = len(mice_eradication)
    n_mice_reduction = len(mice_reduction)

    # --------------------
    # tumor relapse
    # definition:
    # take the maximum tumor size at aligned day <= 9 (day of final treatment)
    # Relapse: a mouse that has had either reduction or eradication has a 20% increase of their minimum tumor size after
    # the minimum tumor size measurement, and the tumor is bigger than 9 mm2 at the relapse measurement.
    # --------------------
    minimal_increase_for_relapse = 0.20  # fraction of tumor size increase from nadir to count as relapse.
    relapse_min_absolute = 9 # mm2 of tumor to count as relapse.

    mice_relapse = []
    relapse_time_frames = {}

    responders = set(mice_eradication).union(set(mice_reduction))  # so smaller set now with 20%

    for mouse in responders:
        # Align mouse series for responders so treatment start is day 0 (repeat for current subset)
        mouse_series = data_df.loc[mouse].copy()
        mouse_series.index = mouse_series.index.astype(int)
        mouse_series.index = mouse_series.index - (cyclo_days[mouse] + 1)

        mouse_nonan = mouse_series.dropna()
        if len(mouse_nonan) < 3:
            continue

        # define response anchor as before:
        # max tumor size at aligned day <= 9
        up_to_d9 = mouse_series[mouse_series.index <= max_treatment_day].dropna()
        if len(up_to_d9) == 0:
            continue

        max_tumor_size = up_to_d9.max()
        max_day = up_to_d9[up_to_d9 == max_tumor_size].index[-1]

        # look for first response after max_day:
        # 2 consecutive measurements at least 20% below the max
        post_max = mouse_series[mouse_series.index > max_day].dropna()
        if len(post_max) < 2:
            continue

        reduction_threshold = max_tumor_size * (1 - minimal_reduction)
        is_reduced = post_max < reduction_threshold

        response_day = None
        post_max_idx = post_max.index.to_list()

        for i in range(1, len(post_max)):
            if is_reduced.iloc[i - 1] and is_reduced.iloc[i]:
                response_day = post_max_idx[i - 1]  # first day of the first qualifying pair
                break

        # relapse is evaluated after response_day
        post_response = mouse_series[mouse_series.index >= response_day].dropna() # update: now includes the first
        # day of the response. if that is the lowest, it can be the nadir.
        if len(post_response) < 2:
            continue

        nadir = post_response.min() # can be any day post max, including the response day itself
        nadir_day = post_response[post_response == nadir].index[-1] # take last in case of plateau nadir.

        relapse_threshold = nadir * (1 + minimal_increase_for_relapse)

        # only evaluate relapse strictly after nadir
        post_nadir = post_response[post_response.index > nadir_day].dropna()
        if len(post_nadir) < 2:
            continue

        is_relapse = (post_nadir >= relapse_threshold) & (post_nadir >= relapse_min_absolute)

        relapse_day = None
        post_nadir_idx = post_nadir.index.to_list()

        # for i in range(1, len(post_nadir)):
        #     if is_relapse.iloc[i - 1] and is_relapse.iloc[i]:
        #         relapse_day = post_nadir_idx[i - 1]  # first day of first qualifying relapse pair
        #         break
        for i in range(1, len(post_nadir)):
            if is_relapse.iloc[i - 1] and is_relapse.iloc[i]:
                candidate_relapse_day = post_nadir_idx[i - 1]  # first day of first qualifying relapse pair

                # require monotonic increase (plateau allowed) from candidate relapse day onward
                relapse_tail = post_nadir.loc[candidate_relapse_day:]
                is_monotonic_non_decreasing = (relapse_tail.diff().dropna() >= 0).all()

                if is_monotonic_non_decreasing:
                    relapse_day = candidate_relapse_day
                    break

        if relapse_day is not None:
            mice_relapse.append(mouse)
            relapse_time_frames[mouse] = {
                'response_day': response_day,
                'relapse_day': relapse_day,
                'days_to_relapse': relapse_day - response_day,
                'nadir': nadir,
                'relapse_threshold': relapse_threshold
            }

    print(f'Relapsed mice: {mice_relapse}')
    print('Relapse details:')
    for mouse, info in relapse_time_frames.items():
        print(mouse, info)

    n_mice_relapse = len(mice_relapse)

    if n_mice_relapse > 0:
        relapse_durations = [info['days_to_relapse'] for info in relapse_time_frames.values()]
        relapse_mean = np.mean(relapse_durations)
        relapse_sd = np.std(relapse_durations, ddof=1) if len(relapse_durations) > 1 else 0.0
        relapse_time_frame_summary = f'{relapse_mean:.1f} ± {relapse_sd:.1f}'
    else:
        relapse_time_frame_summary = np.nan

    df_updated_growth_summary = pd.concat(
        [
            df_updated_growth_summary,
            pd.DataFrame({
                'Treatment': [titles_tumour_growth_incl_tb[df_index]],
                'Total mice': [n_mice],
                'Mice without tumor growth': [n_mice_zero],
                'Mice with tumor eradication': [n_mice_eradication],
                'Mice with tumor reduction': [n_mice_reduction],
                'Mice with relapse': [n_mice_relapse],
                'Relapse time frame': [relapse_time_frame_summary]
            })
        ],
        ignore_index=True
    )
treatment_order_final = ['Tumor bearing', 'Cy', 'ACT day 3', 'ACT day 7', 'ACT day 14', 'Relapse',  'Untreated', 'ACT day 7, no CpG', 'ACT day 14, no CpG']
df_updated_growth_summary['Treatment'] = pd.Categorical(df_updated_growth_summary['Treatment'], categories=treatment_order_final, ordered=True)
df_updated_growth_summary = df_updated_growth_summary.sort_values('Treatment')
#%% add row 'total' with the sums of the counts
total_row = {
    'Treatment': 'Total',
    'Total mice': df_updated_growth_summary['Total mice'].sum(),
    'Mice without tumor growth': df_updated_growth_summary['Mice without tumor growth'].sum(),
    'Mice with tumor eradication': df_updated_growth_summary['Mice with tumor eradication'].sum(),
    'Mice with tumor reduction': df_updated_growth_summary['Mice with tumor reduction'].sum(),
    'Mice with relapse': df_updated_growth_summary['Mice with relapse'].sum(),
    'Relapse time frame': df_updated_growth_summary['Relapse time frame'].dropna().tolist()[0]  # only relapse.
}
# add to df as final row
df_updated_growth_summary = pd.concat([df_updated_growth_summary, pd.DataFrame(total_row, index=[0])], ignore_index=True)

# save summary table to file
df_updated_growth_summary.to_csv(f'{data_repo_path}/processed/tumor_growth/tumor_growth_response_summary_eradication_relapse_response_until_{max_treatment_day}_{file_modifier}.csv', index=False)