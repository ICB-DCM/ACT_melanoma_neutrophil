#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
# '''=================================================
# @Project: ACT_melanoma_neutrophil (data: act_neutrophil_data_repository)
# @Author : Gemma van der Voort
# @Desc   : a subset of the dotplot so that the paper figure is not overcrowded.
# '''=================================================
import sys
from pathlib import Path

import scanpy as sc
import anndata as ad
import matplotlib.pyplot as plt
import seaborn as sns
import os

PROJECT_ROOT = Path(__file__).resolve()
while PROJECT_ROOT.name != "ACT_melanoma_neutrophil" and PROJECT_ROOT.parent != PROJECT_ROOT:
    PROJECT_ROOT = PROJECT_ROOT.parent
if PROJECT_ROOT.name != "ACT_melanoma_neutrophil":
    raise RuntimeError("ACT_melanoma_neutrophil root not found")
sys.path.append(str(PROJECT_ROOT / "helper_files"))

from paths_parameters import cell_colours_hex, processed_data_dir, fig_dir, modifier_base, modifier

fig_dir = f'{fig_dir}/phenotyping/'
os.makedirs(fig_dir, exist_ok=True)
#%%
adata_raw = ad.read_h5ad(f"{processed_data_dir}/phenotyping_{modifier}.h5ad")
#%%
adata = adata_raw.copy()
adata.obs['Cell type'] = adata.obs['Cell type'].replace({'trans T cycling': 'Cycling\npmel-1 CD8 T',
                                                         'trans CD8': 'Pmel-1 CD8 T',
                                                         'endo T cycling': 'Cycling endo CD8 T',
                                                         'endo CD8': 'Endo CD8 T'})
cells = ['B-cell', 'Neutrophil', 'Endo CD8 T', 'Pmel-1 CD8 T', 'Cycling\npmel-1 CD8 T']
markers = ['CD45', 'B220', 'Ly6G', 'CD8', 'CD90.1', 'Ki67']
groupby_key = 'Cell type'
# subset adata to only include the cells of interest
adata_subset = adata[adata.obs['Cell type'].isin(cells)].copy()
#%%
fig, ax = plt.subplots()
sc.pl.dotplot(adata_subset, var_names=markers, groupby=groupby_key, show=False)
sns.despine()
plt.savefig(f"{fig_dir}/dotplot_subsection_new_names_no_space_{modifier}.pdf", dpi=300, bbox_inches='tight')
plt.close()
