#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
# '''=================================================
# @Project: ACT_melanoma_neutrophil (data: act_neutrophil_data_repository)
# @Author : Gemma van der Voort
# @Desc   : preprocessing using the flowjo exported files for the early gates: all cells, all live cells, and all cd45+ cells.
# @Desc updated: Preprocess early-gate flow cytometry exports into consolidated count tables.
# this will be shown as a part of figure 5.
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


def extract_early_gates(filepath, experiment, panel, condition_assigner, organ_fraction_LN,
                        file_col='Unnamed: 0',
                        all_cell_col="All cells | Count",
                        all_lym = "All cells/Lymphocytes | Count",
                        all_single = "All cells/Lymphocytes/Single Cells | Count",
                        all_live = "All cells/Lymphocytes/Single Cells/Live | Count",
                        all_cd45 = "All cells/Lymphocytes/Single Cells/Live/CD45.2+ | Count",
                        all_beads = 'Beads | Count'
                        ):
    df_raw = pd.read_csv(filepath)

    df = pd.DataFrame(columns=['Experiment', 'File', 'Condition', 'Organ', 'Mouse_ID', 'Panel',
                               'all_count', 'lym_count', 'single_count', 'live_count', 'cd45_count',
                               'bead_count', 'Organ_fraction_LN'])

    df['File'] = df_raw[file_col]
    df['Experiment'] = experiment
    df['Panel'] = panel
    df[['Organ', 'Mouse_ID']] = df['File'].apply(split_filename)
    df['Organ'] = df['Organ'].str.split('_').str[0]
    df['Organ'] = normalise_organ_names(df['Organ'])
    df['Organ_fraction_LN'] = organ_fraction_LN

    # Print dropped rows with NaN Organ
    # dropped = df[df['Organ'].isna()]
    # if not dropped.empty:
    #     print(f"[{experiment}] Dropped rows due to NaN Organ field:")
    #     for fname in dropped['File']:
    #         print(f" - {fname}")
    df = df.dropna(subset=['Organ'])

    # Assign condition: function or fixed string
    if callable(condition_assigner):
        df['Condition'] = df['Mouse_ID'].apply(condition_assigner)
    else:
        df['Condition'] = condition_assigner

    df['all_count'] = convert_column_to_int(df_raw[all_cell_col])
    df['lym_count'] = convert_column_to_int(df_raw[all_lym])
    df['single_count'] = convert_column_to_int(df_raw[all_single])
    df['live_count'] = convert_column_to_int(df_raw[all_live])
    df['cd45_count'] = convert_column_to_int(df_raw[all_cd45])
    df['bead_count'] = convert_column_to_int(df_raw[all_beads])
    return df

def assign_condition_102(mouse_id):
    no_CpG_d14_mice = ['HOE-2479', 'HOE-2489', 'HOE-2474', 'HOE-2482', 'HOE-2488', 'HOE-2473', 'HOE-2477', 'HOE-2470',
                       'HOE-2481']
    return 'ACT-d14_no CpG' if mouse_id in no_CpG_d14_mice else 'ACT-d14'

def convert_column_to_int(col):
    col = col.astype(str).str.replace(',', '.', regex=False)
    return pd.to_numeric(col, errors='coerce') # .astype('Int64') nullable integer dtype


raw_data_path = f'{data_repo_path}/raw/flow_cytometry/early_gates'
processed_data_path = f'{data_repo_path}/processed/flow_cytometry'
files = os.listdir(raw_data_path)
#%%
experiments_flow = ['069', '070', '076', '086', '096', '097', '087', '102', '103']
conditions_flow = ['Naive', 'Tumour-bearing 3-5 mm', 'Cyclo only', 'ACT-d3', 'ACT-d7',
       'ACT-d14', 'ACT_Relapse', assign_condition_102, 'ACT-d7_no CpG']
panel_amount_flow = [2,2,2,3,3,3,3,3,3]  # control experiments don't yet have pmel panel. fractions for LNs.
dict_conditions = {exp: cond for exp, cond in zip(experiments_flow, conditions_flow)}
dict_exp_panel = {exp: panel for exp, panel in zip(experiments_flow, panel_amount_flow)}
#%%
df_to_concat = pd.DataFrame(columns=['Experiment', 'File', 'Condition', 'Organ', 'Mouse_ID', 'Panel',
                               'all_count', 'lym_count', 'single_count', 'live_count', 'cd45_count',
                                     'bead_count', 'Organ_fraction_LN'])
for file in files:
    if file.endswith('.csv') and 'tumour' not in file:  # only interested in the LN, not tumour
        experiment = file[:3]
        panel = file[:-4].split('_')[-1]
        organ_fraction = 1/int(dict_exp_panel[experiment])
        print(file)
        file_path = f'{raw_data_path}/{file}'
        condition_assigner = dict_conditions[experiment]
        if file == '0871_slo_cd8.csv':  # typo in col names
            all_cell_col = "All cellls | Count"
            all_lym = "All cellls/Lymphocytes | Count"
            all_single = "All cellls/Lymphocytes/Single Cells | Count"
            all_live = "All cellls/Lymphocytes/Single Cells/Live | Count"
            all_cd45 = "All cellls/Lymphocytes/Single Cells/Live/Cd45.2+ | Count"
            df_1file = extract_early_gates(file_path, experiment, panel, condition_assigner,
                                     all_cell_col=all_cell_col,
                                     all_lym=all_lym,
                                     all_single=all_single,
                                     all_live=all_live,
                                     all_cd45=all_cd45,
                                     organ_fraction_LN=organ_fraction)
        elif file == '0861_slo_panimmune.csv':  # typo in col name CD45
            all_cd45 = "All cells/Lymphocytes/Single Cells/Live/Cd45.2+ | Count"
            df_1file = extract_early_gates(file_path, experiment, panel, condition_assigner,
                                     all_cd45=all_cd45, organ_fraction_LN=organ_fraction)
        elif file == '0702_slo_cd8.csv':
            all_beads = 'beads | Count'
            df_1file = extract_early_gates(file_path, experiment, panel, condition_assigner,
                                     all_beads=all_beads, organ_fraction_LN=organ_fraction)
        elif file == '0861_slo_pmel.csv' or file == '0871_slo_pmel.csv':
            print('missing beads, not added in experiment. Filling column with nans.')
            all_beads = 'Unnamed: 0'  # filling with file names to prevent error, replaced with nans later.
            df_1file = extract_early_gates(file_path, experiment, panel, condition_assigner,
                                     all_beads=all_beads, organ_fraction_LN=organ_fraction)
            df_1file['bead_count'] = None
        else:
            df_1file = extract_early_gates(file_path, experiment, panel, condition_assigner,
                                           organ_fraction_LN=organ_fraction)

    df_to_concat = pd.concat([df_to_concat, df_1file])
#%% correct organ fraction for some very large LNs in the relapse cohort.
# either half or 2/3rd was kept as material for other experiments. fraction is the same for each panel.
corrected_organ_fraction_dict = {'brLN-r_HOE-2283.fcs': 1 / 9,
                               'brLN-r_HOE-2288.fcs': 1/6,
                               'inLN-l_HOE-2283.fcs': 1/6,
                               'inLN-r_HOE-2270.fcs': 1/6,
                               'inLN-r_HOE-2272.fcs': 1/6,
                               'inLN-r_HOE-2277.fcs': 1/6,
                               'inLN-r_HOE-2283.fcs': 1/9,
                               'inLN-r_HOE-2284.fcs': 1/6,
                               'inLN-r_HOE-2288.fcs': 1/6,
                               'inLN-r_HOE-2292.fcs': 1/6,
                               'inLN-r_HOE-2294.fcs':1/6,
                                 }
for file, fraction in corrected_organ_fraction_dict.items():
    df_to_concat.loc[df_to_concat['File'] == file, 'Organ_fraction_LN'] = fraction
    print(f'corrected {file} to {fraction}')
lymph_nodes = ['inLNr', 'inLNl', 'brLNr']
# change organ fractions to nan when organ is not in lymph nodes (since only meaningful for them)
df_to_concat['Organ_fraction_LN'] = df_to_concat['Organ_fraction_LN'].where(df_to_concat['Organ'].isin(lymph_nodes))

#%% found a typo in mouse BAL-3427, (only brLnr, pan-immune): should be BAL-3727 (only brLNr, pan-immune missing)
df_to_concat.loc[df_to_concat['Mouse_ID'] == 'BAL-3427', 'Mouse_ID'] = 'BAL-3727'
#%% update: rename BAL-4028 in FLOW to BAL-4058 for the lymph nodes
df_to_concat.loc[(df_to_concat['Organ'].isin(['inLNr', 'inLNl', 'brLNr'])) & (df_to_concat['Mouse_ID'] == 'BAL-4028'),
                'Mouse_ID'] = 'BAL-4058'
# %% save
df_to_concat = df_to_concat.reset_index(drop=True)
df_to_concat.to_csv(f'{processed_data_path}/early_gates_slo_counts.csv', index=False)

#%% troubleshooting column names
df_issue = pd.read_csv(f'{raw_data_path}/0861_slo_pmel.csv')
with open(f'{processed_data_path}/0861_slo_pmel_columns.txt', 'w') as f:
    for col in df_issue.columns:
        f.write(f'{col}\n')
df_86_slo_pmel = pd.read_csv(f'{raw_data_path}/0861_slo_pmel.csv')
df_86_slo_pmel_slim = df_86_slo_pmel.loc[:, ~df_86_slo_pmel.columns.str.contains(r'\.\d{1,2}$')]
