#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
# '''=================================================
# @Project: ACT_melanoma_neutrophil (data: act_neutrophil_data_repository)
# @Author : Gemma van der Voort
# @Desc   : used for preprocessing the flow cytometry data the cd8/pmel panels with Maike's traditional gating.
# @Desc updated: Flow cytometry preprocessing utilities.
# '''=================================================
import numpy as np
import pandas as pd
import re


def normalise_organ_names(organ_col: pd.Series):
    """
    replaces unusual spelling of organ names in the organ_col with the standard spelling.
    :param organ_col: pd.Series with the organ names.
    :return: organ_col with normalised names.
    """
    spelling_changes = {'brLNr': ['brLN-r', 'brLN-right', 'brLNri', 'br-LN-r'], 'inLNl': ['inLN-l', 'inLN-left', 'inLNle', 'in-LN-l'],
                        'inLNr': ['inLN-r', 'inLN-right', 'inLNri', 'in-LN-r'], 'Spleen': ['S', 'spleen'], 'Blood': ['B', 'blood'],
                        'Tumour': ['T', 'tumour']}
    reverse_map = {variant: canon for canon, variants in spelling_changes.items() for variant in variants}
    return organ_col.replace(reverse_map, inplace=False)


def split_filename(s, nan_mismatches=True):
    """
    split the filename into organ and mouse_id. example: 'brLN-r_HOE-2228.fcs' -> 'brLN-r', 'HOE-2228'.
    but also something like "B-HOE-2226_d26_Sac.fcs" -> 'B', 'HOE-2226'. Code is robust to extra notes and dashes
    instead of underscores.
    :param s: the filename to split.
    :param nan_mismatches: if True, return NaN if the split was not successful. If False, raise ValueError.
    :return:
    """
    match = re.split(r'(HOE|BAL)', s)
    if len(match) >= 3:
        organ = match[0].rstrip('-_')
        mouse_id = match[1] + match[2].split('.')[0]
        # keep only first 8 characters of mouse_id in case there are extra notes in the file name (blood day etc.)
        mouse_id = mouse_id[:8]
    else:
        # return NaN if the split was not successful
        if nan_mismatches:
            print(f'Warning: could not split {s} into organ and mouse_id, returning NaN')
            return pd.Series([np.nan, np.nan])
        else:
            raise ValueError(f'Could not split {s} into organ and mouse_id, missing substring HOE or BAL')
    return pd.Series([organ, mouse_id])


def preprocess_cd8_panel_flow_endogenous(file_location, sheet_name='Endo CD8 T'):
    """
    Preprocess the flow cytometry data from Maike's traditional gating, for the cd8 panel, endogenous cells.
    :param file_location:
    :return:
    """
    cd8_df = pd.read_excel(file_location, sheet_name=sheet_name)
    cd8_df = cd8_df.rename(columns={cd8_df.columns[0]: 'File_name'})
    # split file name into a part before substring HOE or BAL and after
    cd8_df[['Organ', 'Mouse_ID']] = cd8_df['File_name'].apply(split_filename)
    cd8_df['Organ'] = normalise_organ_names(cd8_df['Organ'])

    # shorten column names by removing common gating prefix. (easier to read but not needed.
    # cd8_df = cd8_df.rename(
    #     columns={col: col.replace('All cells/Lymphocytes/Single Cells/Live/CD45.2+/CD8+ T cells/', '')
    #              for col in cd8_df.columns})
    # for columns that contain the word 'Freq.' in the name, remove all '%' from the values and convert to float
    freq_cols = [col for col in cd8_df.columns if 'Freq.' in col]
    cd8_df[freq_cols] = cd8_df[freq_cols].replace(' %', '', regex=True).replace(',', '.', regex=True).astype(float)
    # Define column patterns of interest in normalized form (marker order can be swapped)
    marker_patterns = [
        ('cd62l-', 'cd44+', '| count', 'Teff_count'),
        ('cd62l-', 'cd44+', '| freq. of cd45.2+', 'Teff_freq_CD45'),
        ('cd62l-', 'cd44+', '| freq. of endogenous cd8 t', 'Teff_freq_endo'),
        ('cd62l-', 'cd44+', '| freq. of parent', 'Teff_freq_endo'), # alternative naming
        ('cd44+', 'cd62l+', '| count', 'Tcm_count'),
        ('cd44+', 'cd62l+', '| freq. of cd45.2+', 'Tcm_freq_CD45'),
        ('cd44+', 'cd62l+', '| freq. of endogenous cd8 t', 'Tcm_freq_endo'),
        ('cd44+', 'cd62l+', '| freq. of parent', 'Tcm_freq_endo') # alternative naming
    ]
    # Base columns to keep
    cols_of_interest = ['Organ', 'Mouse_ID']
    cd8_df_slim = cd8_df[cols_of_interest].copy()

    column_map = {col.lower(): col for col in cd8_df.columns}  # map lowercase -> actual
    lowercase_cols = list(column_map.keys())
    # Try to find and match column names even if marker order is swapped
    matched_columns = set()
    for m1, m2, suffix, new_col in marker_patterns:
        for lc_col in lowercase_cols:
            if all(part in lc_col for part in [m1, m2, suffix]):
                actual_col = column_map[lc_col]
                cd8_df_slim.loc[:, new_col] = cd8_df[actual_col]
                matched_columns.add(new_col)
                break  # only assign once per target

    # Warn about missing matches
    expected_columns = set([m[-1] for m in marker_patterns])
    for col in expected_columns - matched_columns:
        print(f'Warning: {col} not found in columns of {file_location}')
    return cd8_df_slim


def preprocess_cd8_panel_flow_pmel(file_location, sheet_name='Pmel T'):
    """
    Preprocess the flow cytometry data from Maike's traditional gating, for the cd8 panel, pmel T cells.
    :param file_location:
    :return:
    """
    pmel_df = pd.read_excel(file_location, sheet_name=sheet_name)
    pmel_df = pmel_df.rename(columns={pmel_df.columns[0]: 'File_name'})
    # make column organ by taking the first part of the file_name (before the first underscore)
    pmel_df[['Organ', 'Mouse_ID']] = pmel_df['File_name'].apply(split_filename)
    pmel_df['Organ'] = normalise_organ_names(pmel_df['Organ'])
    # shorten column names by removing common gating prefix.
    # pmel_df = pmel_df.rename(
    #     columns={col: col.replace('All cells/Lymphocytes/Single Cells/Live/CD45.2+/CD8+ T cells/', '')
    #              for col in pmel_df.columns})
    # for columns that contain the word 'Freq.' in the name, remove all '%' from the values and convert to float
    freq_cols = [col for col in pmel_df.columns if 'Freq.' in col]
    pmel_df[freq_cols] = pmel_df[freq_cols].replace(' %', '', regex=True).replace(',', '.', regex=True).astype(float)
    # now make a slimmer df with only
    # Define column patterns of interest in normalized form (marker order can be swapped)
    marker_patterns = [
        ('cd62l-', 'cd44+', '| count', 'Teff_count'),
        ('cd62l-', 'cd44+', '| freq. of cd45.2+', 'Teff_freq_CD45'),
        ('cd62l-', 'cd44+', '| freq. of pmel-1 t', 'Teff_freq_pmel'),
        ('cd62l-', 'cd44+', '| freq. of parent', 'Teff_freq_pmel'),  # alternative naming
        ('cd62l-', 'cd44+', '| freq. of pmel t', 'Teff_freq_pmel'),  # alternative naming
        ('cd44+', 'cd62l+', '| count', 'Tcm_count'),
        ('cd44+', 'cd62l+', '| freq. of cd45.2+', 'Tcm_freq_CD45'),
        ('cd44+', 'cd62l+', '| freq. of pmel-1 t', 'Tcm_freq_pmel'),
        ('cd44+', 'cd62l+', '| freq. of parent', 'Tcm_freq_pmel'),  # alternative naming
        ('cd44+', 'cd62l+', '| freq. of pmel t', 'Tcm_freq_pmel')  # alternative naming
    ]
    # Base columns to keep
    cols_of_interest = ['Organ', 'Mouse_ID']
    pmel_df_slim = pmel_df[cols_of_interest].copy()

    column_map = {col.lower(): col for col in pmel_df.columns}  # map lowercase -> actual
    lowercase_cols = list(column_map.keys())
    # Try to find and match column names even if marker order is swapped
    matched_columns = set()
    for m1, m2, suffix, new_col in marker_patterns:
        for lc_col in lowercase_cols:
            if all(part in lc_col for part in [m1, m2, suffix]):
                actual_col = column_map[lc_col]
                pmel_df_slim.loc[:, new_col] = pmel_df[actual_col]
                matched_columns.add(new_col)
                break  # only assign once per target

    # Warn about missing matches
    expected_columns = set([m[-1] for m in marker_patterns])
    for col in expected_columns - matched_columns:
        print(f'Warning: {col} not found in columns of {file_location}')
    return pmel_df_slim
