#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
# '''=================================================
# @Project -> File :cd8_mel_codex_phenotyping -> 0_initial_phenotyping.py
# @Author : Gemma van der Voort
# @Time   : 08.04.24 14:47
# @Desc   : draft phenotyping according to development version of the SPARQ-MI pipeline. Repository version used is
# available as a zip in the zenodo repository. branch: dev, commit: 6b1707ae. Phenotyping result heavily postprocessed,
# see 1_phenotyping_finetuning.py
# '''=================================================

import os
import sys
import warnings
import time

import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc
import random

PROJECT_ROOT = Path(__file__).resolve()
while PROJECT_ROOT.name != "ACT_melanoma_neutrophil" and PROJECT_ROOT.parent != PROJECT_ROOT:
    PROJECT_ROOT = PROJECT_ROOT.parent
if PROJECT_ROOT.name != "ACT_melanoma_neutrophil":
    raise RuntimeError("ACT_melanoma_neutrophil root not found")
sys.path.append(str(PROJECT_ROOT / "helper_files"))
PIPELINE_PATH = ""
sys.path.append(str(PIPELINE_PATH / "/codexpipe"))

from paths_parameters import processed_data_dir, fig_dir, parent_dir, data_fraction, modifier_base, modifier_finetuning
# %%
import codex_phenotyping.CodexPhenotyper as pheno
from phenotyping_tree_codex import tree
random.seed(42)
start_time = time.time()
warnings.filterwarnings("ignore")
modifier = f'{modifier}'
print(modifier)

exp_list = [os.path.join(parent_dir, f) for f in os.listdir(parent_dir) if os.path.isdir(os.path.join(parent_dir, f))]
exp_list = [f for f in exp_list if not f.endswith("_failed")]  # remove folders ending in _failed
exp_list = [f for f in exp_list if not f.endswith("_skip")]  # remove folders ending in _skip
# %% initialise
phenotyper = pheno.CodexPhenotyper(exp_list, compartment_tree=tree)

# %% parameters
failed_markers = ["Tim3", "CD137", "Granzyme B"]
phenotyper.filter_parameters |= {
    'frac': 1,  #  low for testing
    "a_max": 1200,
    "a_min": 100,
    "min_score": -np.inf,
    "empty_datapoint_threshold": 1,  # what counts as an empty datapoint
    "empty_datapoint_quantile": 0.99,  # cells without any signal
    "full_datapoint_threshold": 1,
    "full_datapoint_quantile": 1,
    "normalized_empty_cell_threshold": -2,
    "normalized_empty_cell_quantile": 0.99,
    "max_artifact_score": 1225,  # no removal
    "max_artifact_overlap": 1,
    "min_log_pval": 0,  # large can make sense, lots of pixels so p-vals can go low.
    "min_fc": 0,
    "measures_to_normalize": ["mean", "intensity_mean", "foreground-mean", "background-mean"],
    "adata_main_measure": "fc",
    "failed_markers": failed_markers,
}
phenotyper.clustering_parameters |= {
    "resolution": 30,  # yields 1100 clusters
    "subset_fraction": 1,  # 0.1 for full dataset, 0.5 for 5% subset. 1 default
}

phenotyper.expression_group_parameters |= {
            "min_diff": 0.5,  # default 1. mixed cd90.1 expression too much 0.1 couldn't find endo cd8. 0.5 try now
            "n_processes": 1  # 1 works, higher crashes. run overnight in 12h
            }

# %%
phenotyper.compartment_tree = tree
phenotyper.run_clustering()
phenotyper.save(f"{processed_data_dir}/phenotyping_{modifier}_intermediate.pheno", adata=True, filtered_dfs=False,
                raw_dfs=False)
print(f'Time taken: {(time.time() - start_time)/60:.2f} minutes')

phenotyper.run_cell_type_assignment()
phenotyper.save(f"{processed_data_dir}/phenotyping_{modifier}.pheno", adata=True, filtered_dfs=False,
                raw_dfs=False)
# filtered_dfs and raw_dfs bool arguments to cut down saving time
#%% to file
adata = phenotyper.adata
adata.write(f"{processed_data_dir}/phenotyping_{modifier}.h5ad")
end_time = time.time()
print(f'Time taken: {(end_time - start_time)/60:.2f} minutes')
# %%var names
print(adata.var_names)
# %%visualize umap in adata.obsm['X_umap']
os.makedirs(f'{fig_dir}/raw/{modifier}', exist_ok=True)  # for figs prior to postprocessing
sc.pl.umap(adata, color='exp_name', return_fig=True)
plt.savefig(f"{fig_dir}/{modifier}/raw/umap_{modifier}_exp_name.png", bbox_inches='tight', dpi=300)

sc.pl.umap(adata, color='Compartment', return_fig=True)
plt.savefig(f"{fig_dir}/{modifier}/raw/umap_{modifier}_compartment.png", bbox_inches='tight', dpi=300)

#%% visualise marker intensity across cell types in heatmap
drop_markers = ['NPR', 'NCR', 'CPR', 'nuclear_concavity', 'area']
all_markers = list(adata.var_names)
remaining_markers = [m for m in all_markers if m not in drop_markers]
sc.pl.heatmap(adata, var_names=remaining_markers, groupby='Compartment', cmap='viridis', show=False)
plt.savefig(f"{fig_dir}/{modifier}/raw/heatmap_{modifier}.png", bbox_inches='tight', dpi=300)
# %%same but dotplot
sc.pl.dotplot(adata, var_names=remaining_markers, groupby='Compartment', show=False)
plt.savefig(f"{fig_dir}/{modifier}/raw/dotplot_{modifier}_compartment.png", bbox_inches='tight', dpi=300)
# %%same but dotplot
sc.pl.dotplot(adata, var_names=remaining_markers, groupby='exp_name', show=False)
plt.savefig(f"{fig_dir}/{modifier}/raw/dotplot_{modifier}_exp_name.png", bbox_inches='tight', dpi=300)
