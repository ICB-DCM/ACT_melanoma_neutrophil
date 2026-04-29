#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
# '''=================================================
# @Project -> File :cd8_mel_codex_analysis -> nbh_based_regions_final_incl_vis.py
# @Author : Gemma van der Voort
# @Time   : 06.04.2026
# @Desc   : for revision: metacluster filtering sensitivy analysis. visualisation of checks in script
# revision_1-nbh_sensitivity_analysis.py.
# '''=================================================
# %% imports
import anndata as ad
import pandas as pd
import os
import sys
from pathlib import Path
import warnings

PROJECT_ROOT = Path(__file__).resolve()
while PROJECT_ROOT.name != "ACT_melanoma_neutrophil" and PROJECT_ROOT.parent != PROJECT_ROOT:
    PROJECT_ROOT = PROJECT_ROOT.parent
if PROJECT_ROOT.name != "ACT_melanoma_neutrophil":
    raise RuntimeError("ACT_melanoma_neutrophil root not found")
sys.path.append(str(PROJECT_ROOT / "helper_files"))

from paths_parameters import (processed_data_dir, fig_dir, modifier, modifier_base, parent_dir,
                              processed_spatial_data_dir, data_repo_path)
from utils_codex import spatial_metacluster_filtering_incl_diagnostics

# %% load data
warnings.filterwarnings("ignore")
img_modifier = 'v4_corrected'  # for the data corrections done at the phenotyping level.

version_key = 'v4_8' # post metacluster assigment. This script tests the region smoothing algorithm.
metacluster_key_orig = f'Metacluster {version_key}'
nbh_zone_dict = {'B-cell follicle': [5, 11, 12, 15], 'T-cell zone': [1, 4, 6,9,13,19, 10],  # moved 10 to T cell zone
                      'Medulla-Interfollicular zone-SCS': [0,2,3,7,8,14,16,17, 18]}
tissue_zones = ['B-cell follicle', 'T-cell zone', 'Medulla-Interfollicular zone-SCS']
adata_mod = ad.read_h5ad(f"{processed_spatial_data_dir}/sa_{modifier_base}_inclregions_{version_key}.h5ad")

#%% smoothing parameter sensitivity analysis (really only need iteration 50 since it has all information, but run
# it once to see runtime differences
majority_threshold = 0.5
max_iterations = [11, 12, 13, 14, 15, 20, 25, 30, 50]
diagnostic_dfs = {}
for max_iteration in max_iterations:
    metacluster_key = f"{metacluster_key_orig}_filtered_{majority_threshold}_max_it_{max_iteration}_{version_key}"
    if metacluster_key not in adata_mod.obs.columns: # skip is possible, since it takes about 15 min.
        print(metacluster_key)

        adata_mod, diagnostics_df = spatial_metacluster_filtering_incl_diagnostics(adata_mod, metacluster_key=metacluster_key_orig,
                                                                                   spatial_key='position', n_neighbors=10,
                                                                                   majority_threshold=majority_threshold, max_iterations=max_iteration,
                                                                                   key_added=metacluster_key)
        # save the diagnostics dataframe to csv
        diagnostics_df.to_csv(f"{processed_spatial_data_dir}/revision/metacluster_filtering_diagnostics_{metacluster_key}.csv", index=False)

        diagnostic_dfs[metacluster_key] = diagnostics_df

#%% diagnositics of filtering steps
# iterations needed per image
iterations_per_image = (
    diagnostics_df.groupby('image', as_index=False)['iterations_to_stability']
    .max()
)

# total number of label changes per image
total_changes_per_image = (
    diagnostics_df.groupby('image', as_index=False)['n_changes']
    .sum()
    .rename(columns={'n_changes': 'total_changes'})
)

# only iterations where actual changes happened
diagnostics_changes_only = diagnostics_df[diagnostics_df['n_changes'] > 0].copy()