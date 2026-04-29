#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
# '''=================================================
# @Project: ACT_melanoma_neutrophil (data: act_neutrophil_data_repository)
# @Author : Gemma van der Voort
# @Desc   : visualise the tumor growth data. one plot per experiment with the individual mice trajectories,
# @Desc updated: Merged tumor growth plots and summary lines.
# and a generic plot with the mean and stdev per experiment.
# update 21.03.2025, for v2: offset at first tumour growth
# update 25.03.2025, for v3: offset at first treatment: cyclo. Also annotate other treatments.
# update 02.04.2025, for v4: merge the individual curves in the combined plot to prevent overcrowding. also make a
# version with all lines together
# update 30.6.2025, for v5: include the TB cohorts.
# '''=================================================
import os
import sys
from pathlib import Path

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

# %%style sheet (we are in a subfolder)
plt.style.use(str(PROJECT_ROOT / "helper_files" / "journal_style.mplstyle"))

processed_data_path = f'{data_repo_path}/processed/tumor_growth/per_treatment'
processed_data_path_meta = f'{data_repo_path}/processed/tumor_growth/metadata/per_treatment'
fig_dir = f'{data_repo_path}/figures/tumor_growth_curves/merged'
incl_no_CpG = False
if incl_no_CpG:
    modifier = 'incl_no_CpG'
else:
    modifier = ''
modifier = f'{modifier}_v5_linesstyle_updated'
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
#%% align the data
titles_tumour_growth_incl_tb = titles_tumour_growth[1:] + ["Tumor bearing", "Untreated"]  #"Tumour bearing 8-10mm"]
dfs = [df_cyclo, df_ACT_d3, df_ACT_d7, df_ACT_d14, df_relapse, df_ACT_d7_no_CpG, df_ACT_d14_no_CpG, df_tb3, df_tb8]
meta_dfs = [df_cyclo_meta, df_ACT_d3_meta, df_ACT_d7_meta, df_ACT_d14_meta, df_relapse_meta, df_ACT_d7_no_CpG_meta,
            df_ACT_d14_no_CpG_meta, df_tb3_meta, df_tb8_meta]
dfs_realigned = []
meta_dfs_realigned = []
sacrifice_dates = [0,3,7,14,'Relapse',7,14, -1, 14] # 14 not actual sacrifice date for tb8, but used as axis limit here
for df_index, (data_df, meta_df) in enumerate(zip(dfs, meta_dfs)):
    # Filter out mice that remained below an area of 9 (3x3 mm) for the entire experiment
    n_mice = data_df.shape[0]
    # mice_zero = df[df.max(axis=1) == 0].index
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
        if df_index not in [4]:  # skip for relapse and tb 8
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
#%% plotting: merge all tumor growth data in 1 plot, different variants
# variant 1: individual trajectories as thinner lines and one mean line. (relapse shortened): figS1A
if incl_no_CpG:
    dfs_realigned = dfs_realigned  # keep as is.
    colours_tumour_growth_incl_tb = colours_tumour_growth[1:] + ["#fde725ff", "#fa815fff"]+["#5a5eb8ff", "#5a5eb8ff"]
    titles_tumour_growth_incl_tb = titles_tumour_growth_incl_tb # keep as is.
else:
    dfs_realigned = dfs_realigned[:-4] +  dfs_realigned[-2:]  # remove no CpG experiments
    colours_tumour_growth_incl_tb = colours_tumour_growth[1:] + ["#5a5eb8ff", "#6d8a9fff"]  # tb8 light slate blue
    titles_tumour_growth_incl_tb = titles_tumour_growth_incl_tb[:-4] + titles_tumour_growth_incl_tb[-2:]
fig, ax = plt.subplots(figsize=(10, 6))
treatment_handles = []
treatment_labels = []
for df_index, df_aligned in reversed(list(enumerate(dfs_realigned))):  # so shorter experiments are on top
    print(df_index)
    for mouse_id, mouse_df in df_aligned.groupby('mouse_ID'):
        sns.lineplot(
            data=mouse_df,
            x='day',
            y='value',
            color=colours_tumour_growth_incl_tb[df_index],  # one color per experiment
            linewidth=0.8,
            alpha=0.5,
            ax=ax
        )
    # Plot mean line (ticker)
    mean_df = df_aligned.groupby('day', as_index=False)['value'].mean()
    sns.lineplot(
        data=mean_df,
        x='day',
        y='value',
        color=colours_tumour_growth_incl_tb[df_index],
        label=titles_tumour_growth_incl_tb[df_index],
        linewidth=2.5,
        ax=ax
    )
    treatment_handles.append(ax.lines[-1])
    treatment_labels.append(titles_tumour_growth_incl_tb[df_index])
# Add vertical lines with different dotted styles
event_handles = []
event_labels = []
ax.axvline(x=-1, color='gray', linestyle=(0, (1, 5)), linewidth=1.5, label='Cy')  # dotted
event_handles.append(Line2D([0], [0], color='gray', linestyle=(0, (1, 5)), linewidth=1.5))
event_labels.append('Cy')
ax.axvline(x=0, color='gray', linestyle=(0, (5, 5)), linewidth=1.5, label='Pmels+Virus')
event_handles.append(Line2D([0], [0], color='gray', linestyle=(0, (5, 5)), linewidth=1.5))
event_labels.append('Pmels+Virus')
for CpGday in [3, 6, 9]:
    ax.axvline(x=CpGday, color='gray', linestyle=(0, (3, 5, 1, 5)), linewidth=1.5, label='CpG/PolyIC' if CpGday == 3 else None)
event_handles.append(Line2D([0], [0], color='gray', linestyle=(0, (3, 5, 1, 5)), linewidth=1.5))
event_labels.append('CpG/PolyIC')
# x axis limit
ax.set_xlim(-5, 14)
x_ticks = [-4, -1, 0, 3, 6, 9, 14]
ax.set_xticks(x_ticks)
plt.xlabel("Days (0=treatment start)")
plt.ylabel("Tumor area [mm$^2$]")
# plot legend in original order and add legend headers
phantom = Line2D([0], [0], linestyle='none', color='none')
handles = [phantom] + treatment_handles[::-1] + [phantom] + event_handles
labels = ['Treatment'] + treatment_labels[::-1] + ['Treatment events'] + event_labels
legend = ax.legend(handles=handles, labels=labels, loc='center left', bbox_to_anchor=(1, 0.5))
# Bold the section headers (labels at index 0 and 1 + len(treatment_handles))
legend_texts = legend.get_texts()
header_indices = [0, len(treatment_handles) + 1]
for idx in header_indices:
    legend_texts[idx].set_weight('bold')

plt.tight_layout()
sns.despine()
plt.savefig(f'{fig_dir}/tumor_growth_mean_sd_cyclo_offset_short_time_interval__{file_modifier}_{modifier}_individual_lines.pdf')
plt.close()

#%% variant 2: average the lines of all df_aligned except tb8. plot tb8 separately. add sd as error bar.
fig, ax = plt.subplots(figsize=(10, 6))
df_total_aligned = pd.concat(dfs_realigned[:-2], ignore_index=True)  #remove tb cohorts
line_treated = sns.lineplot(
    data=df_total_aligned,
    x='day',
    y='value',
    color='grey',
    errorbar='sd',
    ax=ax,
    label='Treated'
).lines[-1]  # get the last line object for the legend
line_untreated = sns.lineplot(
dfs_realigned[-1],  # tb8
    x='day',
    y='value',
    # color='blue',  # tb8 light slate blue
    errorbar='sd',
    ax=ax,
    label='Untreated'
).lines[-1]
# Add vertical lines with different dotted styles
event_handles = []
event_labels = []
ax.axvline(x=-1, color='gray', linestyle=(0, (1, 5)), linewidth=1.5, label='Cy')  # dotted
event_handles.append(Line2D([0], [0], color='gray', linestyle=(0, (1, 5)), linewidth=1.5))
event_labels.append('Cy')
ax.axvline(x=0, color='gray', linestyle=(0, (3, 5, 1, 5)), linewidth=1.5, label='Pmels+Virus')
event_handles.append(Line2D([0], [0], color='gray', linestyle=(0, (3, 5, 1, 5)), linewidth=1.5))
event_labels.append('Pmels+Virus')
for CpGday in [3, 6, 9]:
    ax.axvline(x=CpGday, color='gray', linestyle=(5, (10, 3)), linewidth=1.5, label='CpG/PolyIC' if CpGday == 3 else None)
event_handles.append(Line2D([0], [0], color='gray', linestyle=(5, (10, 3)), linewidth=1.5))
event_labels.append('CpG/PolyIC')

# x axis limit
ax.set_xlim(-5, 14)
x_ticks = [-4, -1, 0, 3, 6, 9, 14]
ax.set_xticks(x_ticks)
plt.xlabel("Days (0=treatment start)")
plt.ylabel("Tumor area (mm$^2$)")

# plot legend in original order and add legend headers
treatment_handles_merged = [line_treated, line_untreated]
treatment_labels_merged = ['Treated', 'Untreated']
phantom = Line2D([0], [0], linestyle='none', color='none')
handles = [phantom] + treatment_handles_merged+ [phantom] + event_handles
labels = ['Treatment group'] + treatment_labels_merged + ['Treatment events'] + event_labels
legend = ax.legend(handles=handles, labels=labels, loc='center left', bbox_to_anchor=(1, 0.5))
# Bold the section headers (labels at index 0 and 1 + len(treatment_handles))
legend_texts = legend.get_texts()
header_indices = [0, len(treatment_handles_merged) + 1]
for idx in header_indices:
    legend_texts[idx].set_weight('bold')
plt.tight_layout()
sns.despine()
plt.savefig(f'{fig_dir}/tumor_growth_mean_sd_cyclo_offset_short_time_interval__{file_modifier}_{modifier}_grouped.pdf')
plt.close()
# %% variant 3: add n numbers to variant 2 (fig 1C)

df_aligned_tb8 = dfs_realigned[-1]
dfs_treated = dfs_realigned[:-2]  # all except tb cohorts
# from dfs_realigned[-1], per mouse, drop nans after the final day with data
last_valid_days = (
    df_aligned_tb8.dropna(subset=['value'])            # keep only non-NaN rows
      .groupby('mouse_ID')['day']          # group by mouse
      .max()                               # last day with data per mouse
)

# Merge this info back to the dataframe
df_merged = df_aligned_tb8.merge(
    last_valid_days.rename('last_valid_day'),
    on='mouse_ID',
    how='left'
)

# Keep only rows with day <= last valid day for each mouse
df_cleaned_tb8 = df_merged[df_merged['day'] <= df_merged['last_valid_day']].drop(columns='last_valid_day')

# optional: reset index
df_cleaned_tb8 = df_cleaned_tb8.reset_index(drop=True)
#%%
fig, ax = plt.subplots(figsize=(10, 6))
df_total_aligned = pd.concat(dfs_treated, ignore_index=True)  #remove tb cohorts
line_treated = sns.lineplot(
    data=df_total_aligned,
    x='day',
    y='value',
    color='grey',
    errorbar='sd',
    ax=ax,
    label='Treated'
).lines[-1]  # get the last line object for the legend
line_untreated = sns.lineplot(
    df_cleaned_tb8,  # tb8
    x='day',
    y='value',
    # color='blue',  # tb8 light slate blue
    errorbar='sd',
    ax=ax,
    label='Untreated'
).lines[-1]
# Add vertical lines with different dotted styles
event_handles = []
event_labels = []
ax.axvline(x=-1, color='gray', linestyle=(0, (1, 5)), linewidth=1.5, label='Cy')  # dotted
event_handles.append(Line2D([0], [0], color='gray', linestyle=(0, (1, 5)), linewidth=1.5))
event_labels.append('Cy')
ax.axvline(x=0, color='gray', linestyle=(0, (3, 5, 1, 5)), linewidth=1.5, label='Pmels+Virus')
event_handles.append(Line2D([0], [0], color='gray', linestyle=(0, (3, 5, 1, 5)), linewidth=1.5))
event_labels.append('Pmels+Virus')
for CpGday in [3, 6, 9]:
    ax.axvline(x=CpGday, color='gray', linestyle=(5, (10, 3)), linewidth=1.5, label='CpG/PolyIC' if CpGday == 3 else None)
event_handles.append(Line2D([0], [0], color='gray', linestyle=(5, (10, 3)), linewidth=1.5))
event_labels.append('CpG/PolyIC')
# x axis limit
ax.set_xlim(-5, 14)
x_ticks = [-4, -1, 0, 3, 6, 9, 14]
ax.set_xticks(x_ticks)
plt.xlabel("Days (0=treatment start)")
plt.ylabel("Tumor area [mm$^2$]")

# Compute n numbers at each tick
x_ticks_for_numbers = x_ticks[2:]  # skip -4, as it is not a treatment day
x_ticks_for_numbers = x_ticks # Jan prefers all ticks, so we keep them all
treated_ns = []
untreated_ns = []
# set y lim to -5
ymin, ymax = ax.get_ylim()
ax.set_ylim(0, ymax)
for day in x_ticks_for_numbers:
    n_treated = df_total_aligned[df_total_aligned['day'] == day]['mouse_ID'].nunique()
    treated_ns.append(n_treated)
    n_untreated = df_cleaned_tb8[df_cleaned_tb8['day'] == day]['mouse_ID'].nunique()
    untreated_ns.append(n_untreated)

# Add rows below the x-axis
ymin, ymax = ax.get_ylim()
y_ntext = ymin - 0.37 * (ymax - ymin)  # n row placement, 7% below plot
dy = 0.08 * (ymax - ymin)              # vertical distance between rows

# Add "n Treated:" row
ax.text(x_ticks_for_numbers[0] - 0.5, y_ntext, "n Treated:", ha='right', va='center', fontsize=14, fontweight='bold')
for xi, ni in zip(x_ticks_for_numbers, treated_ns):
    ax.text(xi, y_ntext, f"{ni}", ha='center', va='center', fontsize=14, )

# Add "n Untreated:" row
ax.text(x_ticks_for_numbers[0] - 0.5, y_ntext - dy, "n Untreated:", ha='right', va='center', fontweight='bold', fontsize=14, )
for xi, ni in zip(x_ticks_for_numbers, untreated_ns):
    ax.text(xi, y_ntext - dy, f"{ni}", ha='center', va='center', fontsize=14, )

# extend the y-axis to make space for n numbers
ax.set_ylim(ymin - 0.15 * (ymax - ymin), ymax)


# plot legend in original order and add legend headers
treatment_handles_merged = [line_treated, line_untreated]
treatment_labels_merged = ['Treated', 'Untreated']
phantom = Line2D([0], [0], linestyle='none', color='none')
handles = [phantom] + treatment_handles_merged+ [phantom] + event_handles
labels = ['Treatment group'] + treatment_labels_merged + ['Treatment events'] + event_labels
legend = ax.legend(handles=handles, labels=labels, loc='center left', bbox_to_anchor=(1, 0.5))
# Bold the section headers (labels at index 0 and 1 + len(treatment_handles))
legend_texts = legend.get_texts()
header_indices = [0, len(treatment_handles_merged) + 1]
for idx in header_indices:
    legend_texts[idx].set_weight('bold')
plt.tight_layout()
sns.despine()
plt.savefig(f'{fig_dir}/tumor_growth_mean_sd_cyclo_offset_short_time_interval__{file_modifier}_{modifier}_grouped_incl_n.pdf')
# plt.show()
plt.close()

#%% now average tumor size per experiment. one version with, one without CpG.
# so 6/8 lines in a plot. confidence intervals as a shaded area
fig, ax = plt.subplots(figsize=(10, 6))
treatment_handles = []
treatment_labels = []
if incl_no_CpG:
    dfs_realigned = dfs_realigned  # keep as is.
    colours_tumour_growth_incl_tb = colours_tumour_growth[1:] + ["#fde725ff", "#fa815fff"]+["#5a5eb8ff", "#6d8a9fff"]
    titles_tumour_growth_incl_tb = titles_tumour_growth_incl_tb # keep as is.
else:
    dfs_realigned = dfs_realigned[:-4] +  dfs_realigned[-2:]  # remove no CpG experiments
    colours_tumour_growth_incl_tb = colours_tumour_growth[1:] + ["#5a5eb8ff", "#6d8a9fff"]  # tb8 light slate blue
    titles_tumour_growth_incl_tb = titles_tumour_growth_incl_tb[:-4] + titles_tumour_growth_incl_tb[-2:]

for df_index, df_aligned in reversed(list(enumerate(dfs_realigned))):  # so shorter experiments are on top
    print(df_index)
    # print the name of the experiment
    print(titles_tumour_growth_incl_tb[df_index])
    # Group by day and compute mean and standard deviation
    summary = (
        df_aligned
        .groupby('day')
        .agg(mean_value=('value', 'mean'), std_value=('value', 'std'))
        .reindex(range(-5, 15))  # fill missing days
        .reset_index()  # now 'day' is a column
    )
    if df_index == 3:  # act d14 has some interruptions in the band, interpolate
        summary['mean_value'] = summary['mean_value'].interpolate()
        summary['std_value'] = summary['std_value'].interpolate()
        # print(summary)
    if incl_no_CpG and df_index == 5: # dashed line, yellow for no CpG d7
        sns.lineplot(
            data=df_aligned,
            x='day',
            y='value',
            errorbar='sd',  # measure of spread
            label=titles_tumour_growth_incl_tb[df_index],
            color=colours_tumour_growth_incl_tb[2],
            linestyle= ':',
            ax=ax
    )
    elif incl_no_CpG and df_index == 6: # dashed line, orange for no CpG d14
        sns.lineplot(
            data=df_aligned,
            x='day',
            y='value',
            errorbar='sd',  # measure of spread
            label=titles_tumour_growth_incl_tb[df_index],
            color=colours_tumour_growth_incl_tb[3],
            linestyle=':',
            ax=ax
        )
    else:
        sns.lineplot(
            data=df_aligned,
            x='day',
            y='value',
            errorbar='sd',  # measure of spread
            label=titles_tumour_growth_incl_tb[df_index],
            color=colours_tumour_growth_incl_tb[df_index],
            ax=ax
        )
    ax.fill_between(
        summary['day'],
        summary['mean_value'] - summary['std_value'],
        summary['mean_value'] + summary['std_value'],
        alpha=0.2
    )
    treatment_handles.append(ax.lines[-1])
    treatment_labels.append(titles_tumour_growth_incl_tb[df_index])
# Add vertical lines with different dotted styles
event_handles = []
event_labels = []
ax.axvline(x=-1, color='gray', linestyle=(0, (1, 5)), linewidth=1.5, label='Cy')  # dotted
event_handles.append(Line2D([0], [0], color='gray', linestyle=(0, (1, 5)), linewidth=1.5))
event_labels.append('Cy')
ax.axvline(x=0, color='gray', linestyle=(0, (5, 5)), linewidth=1.5, label='Pmels+Virus')
event_handles.append(Line2D([0], [0], color='gray', linestyle=(0, (5, 5)), linewidth=1.5))
event_labels.append('Pmels+Virus')
for CpGday in [3, 6, 9]:
    ax.axvline(x=CpGday, color='gray', linestyle=(0, (3, 5, 1, 5)), linewidth=1.5, label='CpG/PolyIC' if CpGday == 3 else None)
event_handles.append(Line2D([0], [0], color='gray', linestyle=(0, (3, 5, 1, 5)), linewidth=1.5))
event_labels.append('CpG/PolyIC')
# x axis limit
ax.set_xlim(-20, 60)
ax.set_xticks(x_ticks)
# remove y axis labels below 0 (since growth can't be negative)
ax.set_ylim(bottom=0)
plt.xlabel("Days (0=treatment start)")
plt.ylabel("Tumor area [mm$^2$]")
# plot legend in original order
# handles, labels = ax.get_legend_handles_labels()
# ax.legend(handles[::-1], labels[::-1], title='Treatment')
# reorder labels and handles so TB is at the start
# treatment_labels_reordered = treatment_labels[::-1][] # not finished, nobody was borthered by this so far
phantom = Line2D([0], [0], linestyle='none' , color='none')
handles = [phantom] + treatment_handles[::-1] + [phantom] + event_handles
labels = ['Treatment'] + treatment_labels[::-1] + ['Treatment events'] + event_labels
legend = ax.legend(handles=handles, labels=labels, loc='center left', bbox_to_anchor=(1, 0.5))
# Bold the section headers (labels at index 0 and 1 + len(treatment_handles))
legend_texts = legend.get_texts()
header_indices = [0, len(treatment_handles) + 1]

for idx in header_indices:
    legend_texts[idx].set_weight('bold')
plt.tight_layout()
sns.despine()
plt.savefig(f'{fig_dir}/tumor_growth_mean_sd_cyclo_offset_{modifier}_{file_modifier}.pdf')
plt.close()
#%% with shorter relapse: figS1B/figS1L
fig, ax = plt.subplots(figsize=(10, 6))
treatment_handles = []
treatment_labels = []
for df_index, df_aligned in reversed(list(enumerate(dfs_realigned))):  # so shorter experiments are on top
    print(df_index)
    # Group by day and compute mean and standard deviation
    summary = (
        df_aligned
        .groupby('day')
        .agg(mean_value=('value', 'mean'), std_value=('value', 'std'))
        .reindex(range(-5, 15))  # fill missing days
        .reset_index()  # now 'day' is a column
    )
    if df_index == 3:  # act d14 has some interruptions in the band, interpolate
        summary['mean_value'] = summary['mean_value'].interpolate()
        summary['std_value'] = summary['std_value'].interpolate()



    if incl_no_CpG and df_index == 5: # dashed line, yellow for no CpG d7
        sns.lineplot(
            data=df_aligned,
            x='day',
            y='value',
            errorbar='sd',  # measure of spread
            label=titles_tumour_growth_incl_tb[df_index],
            color=colours_tumour_growth_incl_tb[2],
            linestyle= ':',
            ax=ax
    )
    elif incl_no_CpG and df_index == 6: # dashed line, orange for no CpG d14
        sns.lineplot(
            data=df_aligned,
            x='day',
            y='value',
            errorbar='sd',  # measure of spread
            label=titles_tumour_growth_incl_tb[df_index],
            color=colours_tumour_growth_incl_tb[3],
            linestyle=':',
            ax=ax
        )
    else:
        sns.lineplot(
            data=df_aligned,
            x='day',
            y='value',
            errorbar='sd',  # measure of spread
            label=titles_tumour_growth_incl_tb[df_index],
            color=colours_tumour_growth_incl_tb[df_index],
            ax=ax
        )
    ax.fill_between(
        summary['day'],
        summary['mean_value'] - summary['std_value'],
        summary['mean_value'] + summary['std_value'],
        color = colours_tumour_growth_incl_tb[df_index],
        alpha=0.2
    )
    treatment_handles.append(ax.lines[-1])
    treatment_labels.append(titles_tumour_growth_incl_tb[df_index])
# Add vertical lines with different dotted styles
event_handles = []
event_labels = []
ax.axvline(x=-1, color='gray', linestyle=(0, (1, 5)), linewidth=1.5, label='Cy')  # dotted
event_handles.append(Line2D([0], [0], color='gray', linestyle=(0, (1, 5)), linewidth=1.5))
event_labels.append('Cy')
ax.axvline(x=0, color='gray', linestyle=(0, (5, 5)), linewidth=1.5, label='Pmels+Virus')
event_handles.append(Line2D([0], [0], color='gray', linestyle=(0, (5, 5)), linewidth=1.5))
event_labels.append('Pmels+Virus')
for CpGday in [3, 6, 9]:
    ax.axvline(x=CpGday, color='gray', linestyle=(0, (3, 5, 1, 5)), linewidth=1.5, label='CpG/PolyIC' if CpGday == 3 else None)
event_handles.append(Line2D([0], [0], color='gray', linestyle=(0, (3, 5, 1, 5)), linewidth=1.5))
event_labels.append('CpG/PolyIC')
# x axis limit
ax.set_xlim(-5, 14)
x_ticks = [-4, -1, 0, 3, 6, 9, 14]
ax.set_xticks(x_ticks)
# remove y axis labels below 0 (since growth can't be negative) -5 helps show tumors close to 0
ax.set_ylim(bottom=-5)
plt.xlabel("Days (0=treatment start)")
plt.ylabel("Tumor area [mm$^2$]")
# plot legend in original order
# handles, labels = ax.get_legend_handles_labels()
# ax.legend(handles[::-1], labels[::-1], title='Treatment')
phantom = Line2D([0], [0], linestyle='none', color='none')
handles = [phantom] + treatment_handles[::-1] + [phantom] + event_handles
labels = ['Treatment'] + treatment_labels[::-1] + ['Treatment events'] + event_labels
legend = ax.legend(handles=handles, labels=labels, loc='center left', bbox_to_anchor=(1, 0.5))
# Bold the section headers (labels at index 0 and 1 + len(treatment_handles))
legend_texts = legend.get_texts()
header_indices = [0, len(treatment_handles) + 1]

for idx in header_indices:
    legend_texts[idx].set_weight('bold')
plt.tight_layout()
sns.despine()
plt.savefig(f'{fig_dir}/tumor_growth_mean_sd_cyclo_offset_short_time_interval_{modifier}_{file_modifier}.pdf')
plt.close()
