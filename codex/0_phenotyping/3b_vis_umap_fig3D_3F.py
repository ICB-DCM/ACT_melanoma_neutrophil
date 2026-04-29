#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
# '''=================================================
# @Project: ACT_melanoma_neutrophil (data: act_neutrophil_data_repository)
# @Author : Gemma van der Voort
# @Desc   : umap with new colours and grouped marker channels for fig 2
# '''=================================================
# %% imports
import sys
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import os
import pandas as pd
import scanpy as sc
import matplotlib as mpl
import matplotlib.gridspec as gridspec

PROJECT_ROOT = Path(__file__).resolve()
while PROJECT_ROOT.name != "ACT_melanoma_neutrophil" and PROJECT_ROOT.parent != PROJECT_ROOT:
    PROJECT_ROOT = PROJECT_ROOT.parent
if PROJECT_ROOT.name != "ACT_melanoma_neutrophil":
    raise RuntimeError("ACT_melanoma_neutrophil root not found")
sys.path.append(str(PROJECT_ROOT / "helper_files"))

from paths_parameters import cell_colours_hex, processed_spatial_data_dir, fig_dir, modifier_base, modifier
#%%
fig_dir = f'{fig_dir}/phenotyping/'
os.makedirs(fig_dir, exist_ok=True)
adata = ad.read_h5ad(f"{processed_spatial_data_dir}/phenotyping_{modifier}.h5ad")
#%% plot umap with new colours
# rename 'trans T cycling' to 'Cycling pmel-1 CD8 T' in cell type column etc
replace_dict = {'trans T cycling': 'Cycling pmel-1 CD8 T', 'trans CD8': 'Pmel-1 CD8 T',
                'endo T cycling': 'Cycling endo CD8 T', 'endo CD8': 'Endo CD8 T',
                'endo T exhausted': 'Exhausted endo CD8 T',
                'CD4': 'CD4 T', }
adata.obs['Cell type'] = adata.obs['Cell type'].replace(replace_dict)
# also do for cell_colours_hex dict
cell_colours_hex = {replace_dict.get(k, k): v for k, v in cell_colours_hex.items()}
# update cell_types order after renaming
cell_types = list(cell_colours_hex.keys())
# reorder adata to match the cell type order
adata.obs['Cell type'] = pd.Categorical(
    adata.obs['Cell type'],
    categories=cell_types,
    ordered=True
)

# plot the umap
fig, ax = plt.subplots(1, 1)
ax.axis('off')
sc.pl.umap(adata, color=['Cell type'], show=False, wspace=0.3, palette=cell_colours_hex, frameon=False, ax=ax, title="",
           # add_outline=True
           na_in_legend = False,
           )
plt.title('')
plt.axis('off')
plt.savefig(f'{fig_dir}/umap_new_colours_new_naming_{modifier}.pdf', dpi=300, bbox_inches='tight',)
plt.close()

#%% now a grid with marker umap. no spines or axes, just the umap in a 2x3 grid with marker expression and shared legend
# aspect ratio remains messed-up and individual plots are better anyway, so use this only for the shared legend
markers_to_plot = []# 'Ly6G', 'MHCII',] # 'CD8', 'CD90.1', 'Ki67', 'Tox']
nrows = 2
ncols = 3

# Create figure with GridSpec for better layout control
fig = plt.figure(figsize=(18,10))
gs = gridspec.GridSpec(nrows, ncols + 1, width_ratios=[1] * ncols + [0.08])

vmin, vmax = 0, 1
# flip the plasma colormap to have high expression in dark blue
cmap = mpl.cm.plasma_r
norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
axes = [fig.add_subplot(gs[i // ncols, i % ncols]) for i in range(len(markers_to_plot))]

# # Plot UMAPs
# for i, marker in enumerate(markers_to_plot):
#     ax = axes[i]
#     sc.pl.umap(adata, color=[marker], show=False, frameon=False, ax=ax, title=marker, vmax=vmax, vmin=vmin,
#                colorbar_loc=None)
#     ax.set_box_aspect(1)
#     ax.axis('off')
#     ax.set_title(marker)

# Create shared colorbar
cbar_ax = fig.add_subplot(gs[:, -1])  # Last column, spanning all rows
cb = mpl.colorbar.ColorbarBase(cbar_ax, cmap=cmap, norm=norm, orientation='vertical')
cb.set_label("Expression Level")

# Save figure
plt.savefig(f'{fig_dir}/umap_marker_colormap_only_{modifier}.pdf', dpi=300, bbox_inches='tight')
plt.close()
#%% create a plot per marker
markers_to_plot = ['Ly6G', 'MHCII', 'CD8', 'CD90.1', 'Ki67', 'Tox']
for i_marker, marker in enumerate(markers_to_plot):
    fig, ax = plt.subplots(1, 1)
    sc.pl.umap(adata, color=[marker], show=False, frameon=False, ax=ax, title=marker, vmax=vmax, vmin=vmin,
               colorbar_loc=None, cmap=cmap)
    plt.savefig(f'{fig_dir}/umap_markers_shared_legend_square_{marker}_{modifier}.pdf', dpi=300, bbox_inches='tight')
    plt.close()
