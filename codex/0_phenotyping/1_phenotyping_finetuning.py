#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
# '''=================================================
# @Project: ACT_melanoma_neutrophil (data: act_neutrophil_data_repository)
# @Author : Gemma van der Voort
# @Desc   : finetuning of initial clustering results based on cluster coherence and visual inspection
# note: based on initial phenotyping using in-house code that will soon be available as a preprint.
# in the meantime, please load an anndata object of the data with initial phenotyping layers included.
# '''=================================================

import os
import sys
from pathlib import Path
import matplotlib.pyplot as plt
import scanpy as sc
import random
import anndata as ad

PROJECT_ROOT = Path(__file__).resolve()
while PROJECT_ROOT.name != "ACT_melanoma_neutrophil" and PROJECT_ROOT.parent != PROJECT_ROOT:
    PROJECT_ROOT = PROJECT_ROOT.parent
if PROJECT_ROOT.name != "ACT_melanoma_neutrophil":
    raise RuntimeError("ACT_melanoma_neutrophil root not found")
sys.path.append(str(PROJECT_ROOT / "helper_files"))

from paths_parameters import processed_data_dir, fig_dir, parent_dir, data_fraction, modifier_base, modifier_finetuning

adata_raw = ad.read_h5ad(f'{processed_data_dir}/codex_all_ln_initial_phenotyping.h5ad')

seed = 42
random.seed(seed)

adata = adata_raw.copy()
phenotype_column = 'Compartment'
# %%
print(f'length of adata before umap-based filtering: {len(adata)}')
adata = adata[adata.obsm['X_umap'][:, 0] < 10]
adata = adata[adata.obsm['X_umap'][:, 0] > -10]
adata = adata[adata.obsm['X_umap'][:, 1] < 10]
adata = adata[adata.obsm['X_umap'][:, 1] > -10]
print(f'length of adata after umap-based filtering: {len(adata)}')
#%% add reference to connectivities if missing
adata.uns["neighbors"] = {
    "connectivities_key": "connectivities",
    "distances_key": "distances",
    "params": {}  # optional
}
# %%make fig dir if not exists
fig_dir = f'{fig_dir}/phenotyping/finetuning_not_in_manuscript'
os.makedirs(fig_dir, exist_ok=True)

drop_markers = ['NPR', 'NCR', 'CPR', 'nuclear_concavity', 'area']
all_markers = list(adata.var_names)
remaining_markers = [m for m in all_markers if m not in drop_markers]
# %% create dotplot and umap for the original clustering

sc.pl.dotplot(adata, var_names=remaining_markers, groupby=phenotype_column, show=False)
plt.savefig(f'{fig_dir}/dotplot_original_clustering.pdf')
plt.close()

sc.pl.umap(adata, color=phenotype_column, return_fig=True)
plt.savefig(f'{fig_dir}/umap_original_clustering.pdf', bbox_inches='tight')
plt.close()
# %%print the mixed cell types, the cell types that contain a / in adata.obs['Cell type']
print('Mixed cell types:')
mixed_cell_types = list(adata.obs[phenotype_column].unique())
mixed_cell_types = [c for c in mixed_cell_types if '/' in c]
print(mixed_cell_types)  # --> none, so we can skip finetuning these.

# %% finetune the clustering round 1: split the higher level cell types. here only Lymphoid and T-cell
key = 'leiden_round1'
# remove key if it already exists
if key in adata.obs.keys():
    adata.obs = adata.obs.drop(key, axis=1)
# for high_level_cluster in tqdm(['Lymphoid', 'T-cell']): (loop only adds the last entry)
sc.tl.leiden(adata, seed=seed, restrict_to=(phenotype_column, ['Lymphoid']), resolution=0.3, key_added=key)
sc.tl.leiden(adata, seed=seed, restrict_to=(key, ['T-cell']), resolution=0.3, key_added=key)
# in addition, we saw some strange expression in some of the lower-level clusters, so we will split those as well.
# tackle high ki67 expression in trans T exhausted
sc.tl.leiden(adata, seed=seed, restrict_to=(key, ['trans T exhausted']), resolution=0.3, key_added=key)
# look at cd8 signal in RPMs
sc.tl.leiden(adata, seed=seed, restrict_to=(key, ['RPM']), resolution=0.3, key_added=key)
# split NK compartment to see if there are separate eomes and tbet expressing NK cells, since they shouldn't be both
# high.
sc.tl.leiden(adata, seed=seed, restrict_to=(key, ['NK']), resolution=0.3, key_added=key)
# split neutrophil compartment to find apc neutrophils
sc.tl.leiden(adata, seed=seed, restrict_to=(key, ['Neutrophil']), resolution=0.3, key_added=key)
# split endo T exhausted, the eomes expression seems to be too high.
sc.tl.leiden(adata, seed=seed, restrict_to=(key, ['endo T exhausted']), resolution=0.3, key_added=key)
# finally, split the unclassified compartment
sc.tl.leiden(adata, seed=seed, restrict_to=(key, ['Unclassified']), resolution=0.3, key_added=key)
# show dotplot and umap for the first round of clustering
sc.pl.dotplot(adata, var_names=remaining_markers, groupby=key, show=False)
plt.savefig(f'{fig_dir}/dotplot_{key}.pdf')
plt.close()

sc.pl.umap(adata, color=key, return_fig=True)
plt.savefig(f'{fig_dir}/umap_{key}.pdf', bbox_inches='tight')
plt.close()
# %% dict for new phenotypes, as defined in round 1
# print(list(sorted(adata.obs[key].unique())))
# make a dict with as keys and a values list(sorted(adata.obs[key].unique()))
round_1_dict = {'B-cell': 'B-cell',
                'CD169+ Macrophage': 'CD169+ Macrophage',
                'CD4': 'CD4',
                'FDC': 'FDC',
                'FRC': 'FRC',
                'Lymphoid,0': 'recluster',
                'Lymphoid,1': 'recluster',
                'Lymphoid,2': 'recluster',
                'Lymphoid,3': 'endo CD8',
                'MSM': 'MSM',
                'Mature B-cell': 'Mature B-cell',
                'NK,0': 'NK',
                'NK,1': 'NK',
                'NK,2': 'recluster',
                'NK,3': 'mature NK',
                'NK,4': 'NK',
                'Neutrophil,0': 'Neutrophil',
                'Neutrophil,1': 'Neutrophil',
                'Neutrophil,2': 'Neutrophil',
                'Neutrophil,3': 'Neutrophil',
                'Neutrophil,4': 'Neutrophil',
                'Neutrophil,5': 'Neutrophil',
                'Neutrophil,6': 'Neutrophil',
                'Neutrophil,7': 'APC Neutrophil',
                'RPM,0': 'RPM',
                'RPM,1': 'RPM',
                'RPM,2': 'RPM',
                'RPM,3': 'recluster',
                'RPM,4': 'RPM',
                'RPM,5': 'RPM',
                'T-cell,0': 'recluster',
                'T-cell,1': 'recluster',
                'T-cell,2': 'recluster',
                'T-cell,3': 'recluster',
                'T-cell,4': 'endo T exhausted',
                'T-cell,5': 'recluster',
                'Treg': 'Treg',
                'Unclassified,0': 'recluster',
                'Unclassified,1': 'recluster',
                'Unclassified,10': 'recluster',
                'Unclassified,11': 'recluster',
                'Unclassified,12': 'recluster',
                'Unclassified,2': 'recluster',
                'Unclassified,3': 'recluster',
                'Unclassified,4': 'recluster',
                'Unclassified,5': 'recluster',
                'Unclassified,6': 'recluster',
                'Unclassified,7': 'recluster',
                'Unclassified,8': 'recluster',
                'Unclassified,9': 'recluster',
                'cDC1': 'cDC1',
                'cDC2': 'cDC2',
                'endo CD8': "endo CD8",
                'endo T cycling': "endo T cycling",
                'endo T exhausted,0': 'recluster',
                'endo T exhausted,1': 'endo T exhausted',
                'endo T exhausted,2': 'endo T exhausted',
                'endo T exhausted,3': 'recluster',
                'endo T exhausted,4': 'recluster',
                'endo T exhausted,5': 'endo T exhausted',
                'trans CD8': 'trans CD8',
                'trans T cycling': 'trans T cycling',
                'trans T exhausted,0': 'recluster',
                'trans T exhausted,1': 'recluster',
                'trans T exhausted,2': 'recluster',
                'trans T exhausted,3': 'Noise',
                'trans T exhausted,4': 'Noise',
                }
adata.obs['cell_type_round1'] = (adata.obs[key].map(round_1_dict).astype('category'))
# %%drop 'Noise' from the cell types
# count the amount of cells in noise
print(f'amount of cells in Noise: {len(adata[adata.obs["cell_type_round1"] == "Noise"])}')
adata = adata[adata.obs['cell_type_round1'] != 'Noise']
# %%dotplot and umap for the first round of clustering. add cell amount to dotplot
fig = sc.pl.dotplot(adata, var_names=remaining_markers, groupby='cell_type_round1', show=False)
plt.savefig(f'{fig_dir}/dotplot_round1_finished.pdf')
plt.close()
sc.pl.umap(adata, color='cell_type_round1', return_fig=True)
plt.savefig(f'{fig_dir}/umap_round1_finished.pdf', bbox_inches='tight')
plt.close()

# %%round 2: the dotplot looks pretty good, I would maybe split the NKs again to see if we can find a split between
# eomes and tbet expressing NK cells. apart from that, spit the recluster compartments and see if I can recover the
# trans T exhausted cells.
key = 'leiden_round2'
# remove key if it already exists
if key in adata.obs.keys():
    adata.obs = adata.obs.drop(key, axis=1)
sc.tl.leiden(adata, seed=seed, restrict_to=('cell_type_round1', ['recluster']), resolution=0.3, key_added=key)
sc.tl.leiden(adata, seed=seed, restrict_to=(key, ['NK']), resolution=0.3, key_added=key)
# dotplot and umap for the second round of clustering
sc.pl.dotplot(adata, var_names=remaining_markers, groupby=key, show=False)
plt.savefig(f'{fig_dir}/dotplot_round2.pdf')
plt.close()
sc.pl.umap(adata, color=key, return_fig=True)
plt.savefig(f'{fig_dir}/umap_round2.pdf', bbox_inches='tight')
plt.close()

# %% round 2 dict
round_2_dict = {'B-cell': 'B-cell',
                'APC Neutrophil': 'APC Neutrophil',
                'CD4': 'CD4',
                'CD169+ Macrophage': 'CD169+ Macrophage',
                'FDC': 'FDC',
                'FRC': 'FRC',
                'MSM': 'MSM',
                'Mature B-cell': 'Mature B-cell',
                'NK,0': 'NK',
                'NK,1': 'NK',
                'NK,2': 'NK',
                'NK,3': 'NK',
                'Neutrophil': 'Neutrophil',
                'RPM': 'RPM',
                'Treg': 'Treg',
                'cDC1': 'cDC1',
                'cDC2': 'cDC2',
                'endo CD8': "endo CD8",
                'endo T cycling': "endo T cycling",
                'endo T exhausted': 'endo T exhausted',
                'mature NK': 'mature NK',
                'recluster,0': 'recluster',
                'recluster,1': 'endo CD8',
                'recluster,2': 'Monocyte',
                'recluster,3': 'FRC',
                'recluster,4': 'endo CD8',
                'recluster,5': 'CD169+ Macrophage',
                'recluster,6': 'trans T cycling',
                'recluster,7': 'Mature B-cell',
                'recluster,8': 'trans T exhausted',
                'recluster,9': 'B-cell',
                'recluster,10': 'recluster',
                'recluster,11': 'B-cell',
                'recluster,12': 'recluster',
                'recluster,13': 'recluster',
                'recluster,14': 'CD4',
                'trans CD8': 'trans CD8',
                'trans T cycling': 'trans T cycling',
                }
adata.obs['cell_type_round2'] = (adata.obs[key].map(round_2_dict).astype('category'))
# %%dotplot and umap for the result of second round of clustering. add cell amount to dotplot
fig = sc.pl.dotplot(adata, var_names=remaining_markers, groupby='cell_type_round2', show=False)
plt.savefig(f'{fig_dir}/dotplot_round2_finished.pdf')
plt.close()
sc.pl.umap(adata, color='cell_type_round2', return_fig=True)
plt.savefig(f'{fig_dir}/umap_round2_finished.pdf', bbox_inches='tight')
plt.close()

# %% final round: only split the recluster compartments
key = 'leiden_round3'
# remove key if it already exists
if key in adata.obs.keys():
    adata.obs = adata.obs.drop(key, axis=1)
sc.tl.leiden(adata, seed=seed, restrict_to=('cell_type_round2', ['recluster']), resolution=0.3, key_added=key)

# dotplot and umap for the third round of clustering
sc.pl.dotplot(adata, var_names=remaining_markers, groupby=key, show=False)
plt.savefig(f'{fig_dir}/dotplot_round3.pdf')
plt.close()
sc.pl.umap(adata, color=key, return_fig=True)
plt.savefig(f'{fig_dir}/umap_round3.pdf', bbox_inches='tight')
#%% dict
round_3_dict = {'B-cell': 'B-cell',
                'APC Neutrophil': 'APC Neutrophil',
                'CD4': 'CD4',
                'CD169+ Macrophage': 'CD169+ Macrophage',
                'FDC': 'FDC',
                'FRC': 'FRC',
                'MSM': 'MSM',
                'Mature B-cell': 'Mature B-cell',
                'NK': 'NK',
                'Neutrophil': 'Neutrophil',
                'RPM': 'RPM',
                'Treg': 'Treg',
                'cDC1': 'cDC1',
                'cDC2': 'cDC2',
                'endo CD8': "endo CD8",
                'endo T cycling': "endo T cycling",
                'endo T exhausted': 'endo T exhausted',
                'mature NK': 'mature NK',
                'Monocyte': 'Monocyte',
                'recluster,0': 'trans CD8',
                'recluster,1': 'cDC2',
                'recluster,2': 'recluster',
                'recluster,3': 'recluster',
                'recluster,4': 'recluster',
                'recluster,5': 'recluster',
                'recluster,6': 'recluster',
                'recluster,7': 'CD169+ Macrophage',
                'recluster,8': 'FRC',
                'recluster,9': 'recluster',
                'recluster,10': 'recluster',
                'trans CD8': 'trans CD8',
                'trans T cycling': 'trans T cycling',
                'trans T exhausted': 'trans T exhausted'
                }
adata.obs['cell_type_round3'] = (adata.obs[key].map(round_3_dict).astype('category'))
# %%dotplot and umap for the result of third round of clustering. add cell amount to dotplot
fig = sc.pl.dotplot(adata, var_names=remaining_markers, groupby='cell_type_round3', show=False)
plt.savefig(f'{fig_dir}/dotplot_round3_finished.pdf')
plt.close()
sc.pl.umap(adata, color='cell_type_round3', return_fig=True)
plt.savefig(f'{fig_dir}/umap_round3_finished.pdf', bbox_inches='tight')
plt.close()

# %%print value counts for cell_type_round3
print(adata.obs['cell_type_round3'].value_counts()) # still 46.768 to recluster, give it one last round. What is left
# will be in unclassified compartment
key = 'leiden_round4'
# remove key if it already exists
if key in adata.obs.keys():
    adata.obs = adata.obs.drop(key, axis=1)
sc.tl.leiden(adata, seed=seed, restrict_to=('cell_type_round3', ['recluster']), resolution=0.2, key_added=key)

#%%dotplot and umap for 4th round of clustering
sc.pl.dotplot(adata, var_names=remaining_markers, groupby=key, show=False)
plt.savefig(f'{fig_dir}/dotplot_round4.pdf')
plt.close()
sc.pl.umap(adata, color=key, return_fig=True)
plt.savefig(f'{fig_dir}/umap_round4.pdf', bbox_inches='tight')
plt.close()
# %% dict
round_4_dict = {'B-cell': 'B-cell',
                'APC Neutrophil': 'APC Neutrophil',
                'CD4': 'CD4',
                'CD169+ Macrophage': 'CD169+ Macrophage',
                'FDC': 'FDC',
                'FRC': 'FRC',
                'MSM': 'MSM',
                'Mature B-cell': 'Mature B-cell',
                'NK': 'NK',
                'Neutrophil': 'Neutrophil',
                'RPM': 'RPM',
                'Treg': 'Treg',
                'cDC1': 'cDC1',
                'cDC2': 'cDC2',
                'endo CD8': "endo CD8",
                'endo T cycling': "endo T cycling",
                'endo T exhausted': 'endo T exhausted',
                'mature NK': 'mature NK',
                'Monocyte': 'Monocyte',
                'recluster,0': 'Unclassified',
                'recluster,1': 'Unclassified',
                'recluster,2': 'Unclassified',
                'recluster,3': 'Unclassified',
                'recluster,4': 'Unclassified',
                'recluster,5': 'trans CD8',
                'recluster,6': 'Unclassified',
                'recluster,7': 'endo CD8',
                'trans CD8': 'trans CD8',
                'trans T cycling': 'trans T cycling',
                'trans T exhausted': 'trans T exhausted'
                }
adata.obs['cell_type_round4'] = (adata.obs[key].map(round_4_dict).astype('category'))
print(f'cell_type_round4: {adata.obs["cell_type_round4"].value_counts()}')
print(f'amount of NA in cell_type_round4: {len(adata[adata.obs["cell_type_round4"].isna()])}')
# 43598 unclassified finally (2%), no NA's. I am fine with this.

# %%dotplot and umap
fig = sc.pl.dotplot(adata, var_names=remaining_markers, groupby='cell_type_round4', show=False)
plt.savefig(f'{fig_dir}/dotplot_round4_finished.pdf')
plt.close()
sc.pl.umap(adata, color='cell_type_round4', return_fig=True)
plt.savefig(f'{fig_dir}/umap_round4_finished.pdf', bbox_inches='tight')
plt.close()
#%% save the adata to h5ad
# adata_save = adata.copy()
# adata_save.layers = {}
# adata_save.uns = {}
# adata_save.write(f"{processed_data_dir}/phenotyping_{modifier_base}_layer_uns_rem_finetuned.h5ad")

#%% upon inspection of the postprocessed data, we see that a too high amount of t cycling cells are in the control
# treatments (10 and 20%), and that there are even more trans T exhausted cells in the control treatments than in
# some of the treated samples. I will recluster these two compartments once more to filter out some of these outliers
# if I end up removing all trans T exhausted cells, so be it.
key = 'leiden_round5'
# remove key if it already exists
if key in adata.obs.keys():
    adata.obs = adata.obs.drop(key, axis=1)
    # sc.tl.leiden(adata, seed=seed, restrict_to=('cell_type_round4', ['trans CD8']), resolution=0.5, key_added=key)
sc.tl.leiden(adata, seed=seed, restrict_to=('cell_type_round4', ['trans T cycling']), resolution=0.2, key_added=key)
sc.tl.leiden(adata, seed=seed, restrict_to=(key, ['trans T exhausted']), resolution=0.2, key_added=key)
sc.tl.leiden(adata, seed=seed, restrict_to=(key, ['cDC1']), resolution=0.2, key_added=key)
sc.tl.leiden(adata, seed=seed, restrict_to=(key, ['cDC2']), resolution=0.2, key_added=key)
# dotplot and umap for the fifth round of clustering
sc.pl.dotplot(adata, var_names=remaining_markers, groupby=key, show=False)
plt.savefig(f'{fig_dir}/dotplot_round5v3.pdf')
plt.close()
sc.pl.umap(adata, color=key, return_fig=True)
plt.savefig(f'{fig_dir}/umap_round5v3.pdf', bbox_inches='tight')
plt.close()
# %% dict
round_5_dict = {'B-cell': 'B-cell',
                'APC Neutrophil': 'APC Neutrophil',
                'CD4': 'CD4',
                'CD169+ Macrophage': 'CD169+ Macrophage',
                'FDC': 'FDC',
                'FRC': 'FRC',
                'MSM': 'MSM',
                'Mature B-cell': 'Mature B-cell',
                'NK': 'NK',
                'Neutrophil': 'Neutrophil',
                'RPM': 'RPM',
                'Treg': 'Treg',
                'cDC1,0': 'recluster',
                'cDC1,1': 'recluster',
                'cDC1,2': 'cDC1',
                'cDC1,3': 'cDC1',
                'cDC1,4': 'cDC1',
                'cDC2,0': 'cDC2',
                'cDC2,1': 'cDC2',
                'cDC2,2': 'cDC2',
                'endo CD8': "endo CD8",
                'endo T cycling': "endo T cycling",
                'endo T exhausted': 'endo T exhausted',
                'mature NK': 'mature NK',
                'Monocyte': 'Monocyte',
                'Unclassified': 'Unclassified',
                'trans CD8': 'trans CD8',
                'trans T cycling,0': 'trans T cycling',
                'trans T cycling,1': 'trans T cycling',
                'trans T cycling,2': 'trans T cycling',  # considered making it endo, but that didn't help the noise
                # (now 22% of trans in tb compared to d14 instead of 20%). the dynamics also make more sense this way.
                'trans T cycling,3': 'trans T cycling',
                'trans T exhausted,0': 'Unclassified',
                'trans T exhausted,1': 'Unclassified',
                'trans T exhausted,2': 'Unclassified',
                'trans T exhausted,3': 'Unclassified',
                }
adata.obs['cell_type_round5'] = (adata.obs[key].map(round_5_dict).astype('category'))
print(f'amount of cells in Noise: {len(adata[adata.obs["cell_type_round5"] == "Noise"])}')
adata = adata[adata.obs['cell_type_round5'] != 'Noise']
# %%dotplot and umap
fig = sc.pl.dotplot(adata, var_names=remaining_markers, groupby='cell_type_round5', show=False)
plt.savefig(f'{fig_dir}/dotplot_round5_finished.pdf')
plt.close()
sc.pl.umap(adata, color='cell_type_round5', return_fig=True)
plt.savefig(f'{fig_dir}/umap_round5_finished.pdf', bbox_inches='tight')
plt.close()

#%%
key = 'leiden_round6'
# remove key if it already exists
if key in adata.obs.keys():
    adata.obs = adata.obs.drop(key, axis=1)
sc.tl.leiden(adata, seed=seed, restrict_to=('cell_type_round5', ['recluster']), resolution=0.4, key_added=key)
# dotplot and umap for the sixth round of clustering
# sc.pl.dotplot(adata, var_names=remaining_markers, groupby=key, show=False)
# plt.savefig(f'{fig_dir}/dotplot_round6.pdf')
# plt.close()
# # umap of only clusters starting with recluster or cDC
# sc.pl.umap(adata[adata.obs[key].isin(['recluster', 'cDC1', 'cDC2'])], color=key, return_fig=True)
# plt.savefig(f'{fig_dir}/umap_round6.pdf', bbox_inches='tight')
# plt.close()
print(f'{key}: {adata.obs[key].value_counts()}')
# %% dict
round_6_dict = {'B-cell': 'B-cell',
                'APC Neutrophil': 'APC Neutrophil',
                'CD4': 'CD4',
                'CD169+ Macrophage': 'CD169+ Macrophage',
                'FDC': 'FDC',
                'FRC': 'FRC',
                'MSM': 'MSM',
                'Mature B-cell': 'Mature B-cell',
                'NK': 'NK',
                'Neutrophil': 'Neutrophil',
                'RPM': 'RPM',
                'Treg': 'Treg',
                'cDC1': 'cDC1',
                'cDC2': 'cDC2',
                'endo CD8': "endo CD8",
                'endo T cycling': "endo T cycling",
                'endo T exhausted': 'endo T exhausted',
                'mature NK': 'mature NK',
                'Monocyte': 'Monocyte',
                'Unclassified': 'Unclassified',
                'trans CD8': 'trans CD8',
                'trans T cycling': 'trans T cycling',
                'recluster,0': 'recluster',
                'recluster,1': 'B-cell',
                'recluster,2': 'recluster',
                'recluster,3': 'recluster',
                'recluster,4': 'B-cell',
                'recluster,5': 'cDC1',
                }
adata.obs['cell_type_round6'] = (adata.obs[key].map(round_6_dict).astype('category'))
print(f'cell_type_round6: {adata.obs["cell_type_round6"].value_counts()}')
print(f'amount of NA in cell_type_round6: {len(adata[adata.obs["cell_type_round6"].isna()])}')
# %%dotplot and umap
# fig = sc.pl.dotplot(adata, var_names=remaining_markers, groupby='cell_type_round6', show=False)
# plt.savefig(f'{fig_dir}/dotplot_round6_finished.pdf')
# plt.close()
# sc.pl.umap(adata, color='cell_type_round6', return_fig=True)
# plt.savefig(f'{fig_dir}/umap_round6_finished.pdf', bbox_inches='tight')
# plt.close()
#%% seventh round: only split the recluster compartments
key = 'leiden_round7'
# remove key if it already exists
if key in adata.obs.keys():
    adata.obs = adata.obs.drop(key, axis=1)
sc.tl.leiden(adata, seed=seed, restrict_to=('cell_type_round6', ['recluster']), resolution=0.3, key_added=key)
#dotplot and umap for the seventh round of clustering
sc.pl.dotplot(adata, var_names=remaining_markers, groupby=key, show=False)
plt.savefig(f'{fig_dir}/dotplot_round7.pdf')
plt.close()
sc.pl.umap(adata, color=key, return_fig=True)
plt.savefig(f'{fig_dir}/umap_round7.pdf', bbox_inches='tight')
plt.close()
print(f'{key}: {adata.obs[key].value_counts()}')
#%% dict
round_7_dict = {'B-cell': 'B-cell',
                'APC Neutrophil': 'APC Neutrophil',
                'CD4': 'CD4',
                'CD169+ Macrophage': 'CD169+ Macrophage',
                'FDC': 'FDC',
                'FRC': 'FRC',
                'MSM': 'MSM',
                'Mature B-cell': 'Mature B-cell',
                'NK': 'NK',
                'Neutrophil': 'Neutrophil',
                'RPM': 'RPM',
                'Treg': 'Treg',
                'cDC1': 'cDC1',
                'cDC2': 'cDC2',
                'endo CD8': "endo CD8",
                'endo T cycling': "endo T cycling",
                'endo T exhausted': 'endo T exhausted',
                'mature NK': 'mature NK',
                'Monocyte': 'Monocyte',
                'Unclassified': 'Unclassified',
                'trans CD8': 'trans CD8',
                'trans T cycling': 'trans T cycling',
                'recluster,0': 'recluster',
                'recluster,1': 'B-cell',
                'recluster,2': 'B-cell',
                'recluster,3': 'recluster'
                }
adata.obs['cell_type_round7'] = (adata.obs[key].map(round_7_dict).astype('category'))
# %%dotplot and umap
fig = sc.pl.dotplot(adata, var_names=remaining_markers, groupby='cell_type_round7', show=False)
plt.savefig(f'{fig_dir}/dotplot_round7_finished.pdf')
plt.close()
sc.pl.umap(adata, color='cell_type_round7', return_fig=True)
plt.savefig(f'{fig_dir}/umap_round7_finished.pdf', bbox_inches='tight')
plt.close()
#%% round 8: only split the recluster compartments
key = 'leiden_round8'
# remove key if it already exists
if key in adata.obs.keys():
    adata.obs = adata.obs.drop(key, axis=1)
sc.tl.leiden(adata, seed=seed, restrict_to=('cell_type_round7', ['recluster']), resolution=0.3, key_added=key)
#dotplot and umap for the eighth round of clustering
sc.pl.dotplot(adata, var_names=remaining_markers, groupby=key, show=False)
plt.savefig(f'{fig_dir}/dotplot_round8.pdf')
plt.close()
sc.pl.umap(adata, color=key, return_fig=True)
plt.savefig(f'{fig_dir}/umap_round8.pdf', bbox_inches='tight')
plt.close()
print(f'{key}: {adata.obs[key].value_counts()}')
# %%dict
round_8_dict = {'B-cell': 'B-cell',
                'APC Neutrophil': 'APC Neutrophil',
                'CD4': 'CD4',
                'CD169+ Macrophage': 'CD169+ Macrophage',
                'FDC': 'FDC',
                'FRC': 'FRC',
                'MSM': 'MSM',
                'Mature B-cell': 'Mature B-cell',
                'NK': 'NK',
                'Neutrophil': 'Neutrophil',
                'RPM': 'RPM',
                'Treg': 'Treg',
                'cDC1': 'cDC1',
                'cDC2': 'cDC2',
                'endo CD8': "endo CD8",
                'endo T cycling': "endo T cycling",
                'endo T exhausted': 'endo T exhausted',
                'mature NK': 'mature NK',
                'Monocyte': 'Monocyte',
                'Unclassified': 'Unclassified',
                'trans CD8': 'trans CD8',
                'trans T cycling': 'trans T cycling',
                'recluster,0': 'B-cell',
                'recluster,1': 'B-cell',
                'recluster,2': 'cDC1',
                }
adata.obs['cell_type_round8'] = (adata.obs[key].map(round_8_dict).astype('category'))
# %%dotplot and umap
fig = sc.pl.dotplot(adata, var_names=remaining_markers, groupby='cell_type_round8', show=False)
plt.savefig(f'{fig_dir}/dotplot_round8_finished.pdf')
plt.close()
sc.pl.umap(adata, color='cell_type_round8', return_fig=True)
plt.savefig(f'{fig_dir}/umap_round8_finished.pdf', bbox_inches='tight')
plt.close()
#%% check for typos etc:
key = 'cell_type_round8'
print(f'{key}: {adata.obs[key].value_counts()}')
print(f'amount of NA in {key}: {len(adata[adata.obs[key].isna()])}')
# %% save the adata to h5ad
adata_save = adata.copy()
adata_save.layers = {}
adata_save.uns = {}
adata_save.write(f"{processed_data_dir}/phenotyping_{modifier_base}_{modifier_finetuning}.h5ad")
