#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
# '''=================================================
# @Project: ACT_melanoma_neutrophil (data: act_neutrophil_data_repository)
# @Author : Gemma van der Voort
# @Desc   : for each image, construct a spatial graph and filter it, then save the adata object
# @Desc updated: Spatial graph preprocessing and filtering for CODEx SA.
# '''=================================================
# imports
import os
import sys
from pathlib import Path
import squidpy as sq
import anndata as ad
import matplotlib.pyplot as plt

PROJECT_ROOT = Path(__file__).resolve()
while PROJECT_ROOT.name != "ACT_melanoma_neutrophil" and PROJECT_ROOT.parent != PROJECT_ROOT:
    PROJECT_ROOT = PROJECT_ROOT.parent
if PROJECT_ROOT.name != "ACT_melanoma_neutrophil":
    raise RuntimeError("ACT_melanoma_neutrophil root not found")
sys.path.append(str(PROJECT_ROOT / "helper_files"))

from paths_parameters import fig_dir, modifier, processed_data_dir, modifier_base, processed_spatial_data_dir
from utils_codex import identify_disconnected_points

#%% paths & parameters
fig_dir = f'{fig_dir}/spatial_analysis/nbh_inspection_not_in_manuscript'
os.makedirs(fig_dir, exist_ok=True)
stats_dir = f'{processed_spatial_data_dir}/statistics'
os.makedirs(stats_dir, exist_ok=True)
graph_percentile = 99.5
node_fraction_threshold = 0.01
#%%
adata = ad.read_h5ad(f"{processed_data_dir}/phenotyping_{modifier}.h5ad")
adata_orig = adata.copy()
# adata = adata_orig.copy()  # for troubleshooting
#%% drop unnecessary columns used in phenotyping finetuning
unnecessary_columns = ['artifact_score', 'region', 'leiden', 'Compartment', 'Level_1', 'Level_2', 'Level_3',
                       'leiden_round1', 'cell_type_round1',
                       'leiden_round2', 'cell_type_round2', 'leiden_round3', 'cell_type_round3', 'leiden_round4',
                       'cell_type_round4', 'leiden_round5', 'cell_type_round5', 'leiden_round6', 'cell_type_round6',
                       'leiden_round7', 'cell_type_round7', 'leiden_round8',  'cell_type_round8', 'leiden_round10']
adata.obs = adata.obs.drop(columns=unnecessary_columns)
# %%graph, 6 neighbours default
sq.gr.spatial_neighbors(adata, spatial_key="position", percentile=graph_percentile, library_key='dataset_name')

# %% filter graph
adata = identify_disconnected_points(adata, library_key='dataset_name', node_fraction_threshold=node_fraction_threshold)

# %% plot the first 6 images in 3 cols, 2 rows, the spatial graph with the disconnected points highlighted
fig, axes = plt.subplots(2, 3, figsize=(15, 10))
axes = axes.flatten()
for idx in range(len(adata.obs['dataset_name'].unique()[:6])):
    name = adata.obs['dataset_name'].unique()[idx]
    ax = axes[idx]
    adata_img = adata[adata.obs['dataset_name'] == name]
    sq.pl.spatial_scatter(adata_img, spatial_key='position', shape=None, color='disconnected_point', ax=ax,
                          connectivity_key='spatial_connectivities', library_id=name,
                                  size=0.1,  # dot size
                                  edges_width=0.1,
                          )
    ax.set_title(name)
    handles, labels = ax.get_legend_handles_labels()
    # remove legend and axis labels
    ax.get_legend().remove()
    ax.set_xlabel('')
    ax.set_ylabel('')
# add a legend
fig.legend(handles, labels, loc='upper right')
# set title
fig.suptitle('Spatial graphs with filtering decisions highlighted. Only connected points are considered for '
             'further analysis')
plt.savefig(f"{fig_dir}/{modifier}_disconnected_points.pdf", bbox_inches='tight', dpi=300)
# %%save value counts for disconnected points to csv
adata.obs['disconnected_point'].value_counts().to_csv(f"{stats_dir}/disconnected_points.csv")
# %% drop rows where disconnected_point is not 'connected'
adata = adata[adata.obs['disconnected_point'] == 'connected']
# %%count all cells that are not unclassified
# adata.obs['Cell type'].value_counts()
# %% save
adata.write(f"{processed_spatial_data_dir}/sa_{modifier}_filtered.h5ad")
# %% check if corrections worked.
for id in adata.obs['mouse_id'].unique():  # no HOE-2257, check
    print(id)
# if mouse ID is BAL-3614, print Treatment group
print(adata.obs.loc[adata.obs["mouse_id"] == "BAL-3614", "Treatment"].unique())  # should be cyclo d1
# if mouse id is BAL-3721, print Organs
print(adata.obs.loc[adata.obs["mouse_id"] == "BAL-3721", "Organ"].unique()) # should be inLNr and inLNl
# check if all data was added correctly:
print(adata)
# figs are the same, corrections worked. OK, next step.
