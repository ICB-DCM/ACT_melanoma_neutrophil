#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
# '''=================================================
# @Project: ACT_melanoma_neutrophil (data: act_neutrophil_data_repository)
# @Author : Gemma van der Voort
# @Desc   : Main CODEx neighborhood-based region classification workflow.
# this is the main analysis script for fig 4: nbh classification, assignment to regions, and filtering of regions.
# It exports 2 dfs with the region and filtered region abundances, which are used in the visualisation scripts.
# 4A: clustermap icon from here.
# note: switch to venv-spacec env for this script (python 3.9)
# '''=================================================
# %% imports
import numpy as np
import anndata as ad
import pandas as pd
from tqdm import tqdm
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.patches import Patch
import seaborn as sns
# import spacec as sp
import os
import warnings
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve()
while PROJECT_ROOT.name != "ACT_melanoma_neutrophil" and PROJECT_ROOT.parent != PROJECT_ROOT:
    PROJECT_ROOT = PROJECT_ROOT.parent
if PROJECT_ROOT.name != "ACT_melanoma_neutrophil":
    raise RuntimeError("ACT_melanoma_neutrophil root not found")
sys.path.append(str(PROJECT_ROOT / "helper_files"))

from paths_parameters import (processed_data_dir, fig_dir, modifier, modifier_base, parent_dir,
                              processed_spatial_data_dir, data_repo_path)
from utils_codex import map_cell_type, spatial_metacluster_filtering

# %% load data
# import rbg images (needed if visualise=True)
rgb_images_folder = (f'{data_repo_path}/figures/codex/overview_intensity_images/full_ln_architecture/'
                     'channels_of_interest_highQ')  # input figure_annotation_dir to images
rgb_modifier = 'DAPI_B220_aSMA_CD3_ERTR7'
warnings.filterwarnings("ignore")

visualise = True  # whether to visualise the neighbourhoods and metaclusters on the RGB images
img_to_plot = 6  # number of images to plot if visualise=True
white_background = True  # whether to plot the images with a white background or rgb background

img_modifier = 'v4_corrected'  # for the data corrections done at the phenotyping level.
cell_type_col = f'High level cell type {img_modifier}'  # cell type column in adata_mod
cells_of_interest = ['Neutrophil', 'APC Neutrophil', 'cDC1', 'cDC2']


k = 50  # amount of neighbours. set high for large spatial context
n = 20  # amount of neighbourhoods
cn_col = f"CN_k{k}_n{n}"
min_zone_size = 10


fig_dir = f'{fig_dir}/spatial_analysis/'
os.makedirs(fig_dir, exist_ok=True)
# import adata object with cell types
exp_list = [os.path.join(parent_dir, f) for f in os.listdir(parent_dir) if os.path.isdir(os.path.join(parent_dir, f))]
# remove entries in exp_list that don't start with 'data-reg' (since we added a patches directory called selected...)
exp_list = [f for f in exp_list if os.path.basename(f).startswith('data-reg')]
# drop mice HOE-2257 from exp list
exp_list = [f for f in exp_list if 'HOE-2257' not in f]

#  # %% import adata
# adata_mod = ad.read_h5ad(f"{processed_spatial_data_dir}/sa_{modifier}_filtered.h5ad")  # check if it needs the regions from v2
#
# cd45_cells = set(adata_mod.obs['Cell type'].unique()) - {'FRC', 'FDC', 'Unclassified'}
# # # #%%
# # adata_mod = adata.copy()
# #%%
# high_level_phenotypes_mapping = {'B cell': 'B-cell',
#                                  'Mature B-cell': 'Mature B-cell',
#                                  'CD8 T cell': ['endo CD8', 'endo T cycling', 'endo T exhausted', 'trans CD8',
#                                                 'trans T cycling'],
#                                  # other the same for now:
#                                  'CD169+ M': ['CD169+ M'],
#                                  'MSM': ['MSM'],
#                                  'mature NK': ['mature NK'],
#                                  'FRC': ['FRC'],
#                                  'cDC1': ['cDC1'],
#                                  'cDC2': ['cDC2'],
#                                  'NK': ['NK'],
#                                  'CD4 T cell': ['CD4', 'Treg'],
#                                  'FDC': ['FDC'],
#                                  'Monocyte': ['Monocyte'],
#                                  'Neutrophil': ['Neutrophil', 'APC Neutrophil'],
#                                  # 'Treg': ['Treg'],
#                                  'RPM': ['RPM'],
#                                  # 'APC Neutrophil': ['APC Neutrophil']
#                                  }
# adata_mod.obs['X'] = adata_mod.obsm['position'][:, 0]  # for spacec
# adata_mod.obs['Y'] = adata_mod.obsm['position'][:, 1]  # for spacec
# # %%set high level cell type
# adata_mod.obs[cell_type_col] = adata_mod.obs['Cell type'].apply(map_cell_type,
#                                                                          args=(high_level_phenotypes_mapping,))
# # set unclassified to nan so they are not considered in the neighbourhood analysis
# # adata_mod.obs[cell_type_col] = adata_mod.obs[cell_type_col].replace("Unclassified", np.nan)
# adata_mod.obs[cell_type_col] = adata_mod.obs[cell_type_col].replace(np.nan, "Unclassified")
# #%% with n and k set above
# sp.tl.neighborhood_analysis(
#     adata_mod,
#     unique_region="dataset_name",
#     cluster_col=cell_type_col,
#     X='X',
#     Y='Y',
#     k=k,  # k nearest neighbors
#     n_neighborhoods=n,  # number of CNs
#     elbow=False,
# )
# # %% plot and save the heatmap
# if n<= 20:
#     palette = sns.color_palette("tab20", n)
# else:
#     palette = sns.color_palette("tab20b") + sns.color_palette("tab20c")
# adata_mod.uns[cn_col + '_colors'] = [mcolors.rgb2hex(c) for c in palette[:n]]
# # adata_mod.uns[cn_col + '_colors'] = [mcolors.rgb2hex(rgb) for rgb in get_hls_colors(30)]
# adata_mod.uns[cn_col + '_order'] = adata_mod.obs[cn_col].unique()
# # %% rectangular dataset for heatmap and metacluster finding. x: cell types, y: neighbourhoods
# # get the neighbourhoods
# neighbourhoods = adata_mod.obs[cn_col].unique()
# # get the cell types
# cell_types = adata_mod.obs[cell_type_col].unique()
# # create a dataframe with cell types as columns and neighbourhoods as rows
# df_nbh_ct_counts = pd.DataFrame(index=neighbourhoods, columns=cell_types)
# # fill the dataframe with the neighbourhood enrichment values
# for cell_type in cell_types:
#     for neighbourhood in neighbourhoods:
#         df_nbh_ct_counts.loc[neighbourhood, cell_type] = adata_mod.obs[(adata_mod.obs[cell_type_col] ==
#                                                                         cell_type) & (adata_mod.obs[cn_col] ==
#                                                                                       neighbourhood)].shape[0]
# # convert values to floats
# df_nbh_ct_counts = df_nbh_ct_counts.astype(float)
# # %%counts to log-two-fold change in a 2.5 to -2.5 range over tissue average
# df_nbh_ct_log2fc = np.log2(df_nbh_ct_counts / df_nbh_ct_counts.mean(axis=0))
# assert not df_nbh_ct_log2fc.isnull().values.any()  # check if it contains nans
# # set all values below -2.5 to -2.5
# df_nbh_ct_log2fc[df_nbh_ct_log2fc < -2.5] = -2.5
# # %% seaborn clustermap
# if visualise:
#     fig, ax = plt.subplots(1, 1, figsize=(15, round(n/3)))
#     # cluster both axes and center colormap on 0, increase text size
#     clus_map = sns.clustermap(df_nbh_ct_log2fc, cmap='bwr', figsize=(30, round(n/2)), row_cluster=True, col_cluster=True,
#                               center=0, vmin=-3, vmax=3, row_colors=adata_mod.uns[cn_col + '_colors'],
#                               cbar_kws={'label': 'log2 fold change'})
#     # Set colorbar label size
#     clus_map.ax_cbar.set_ylabel('log2 fold change', size=16)  # Adjust size
#
#     # Set tick label sizes for axes
#     for ax in [clus_map.ax_row_dendrogram, clus_map.ax_col_dendrogram, clus_map.ax_heatmap]:
#         ax.set_xticklabels(ax.get_xticklabels(), fontsize=16)
#         ax.set_yticklabels(ax.get_yticklabels(), fontsize=16)
#
#     # Set tick label sizes for colorbar
#     clus_map.ax_cbar.tick_params(labelsize=16)  # Adjust colorbar tick labels
#     # title
#     plt.suptitle(f"Cellular neighbourhood enrichment log2 fold change for all lymph nodes", fontsize=30)
#     # pad space between title and plot
#     plt.tight_layout(rect=[0, 0.03, 1, 0.95])
#
#     plt.savefig(f"{fig_dir}/cn_exp_heatmap_all_lymph_nodes_{cn_col}_{cell_type_col}_{img_modifier}"
#                 f"_log2fc.pdf", bbox_inches='tight',
#                 dpi=300)
#     plt.close()
# # %% plot and save the scatterplot with the nbhs
# metacluster_key = cn_col  # 'CN_k4_n10'
# if visualise:
#     n_cols = 3
#     n_rows = int(np.ceil(img_to_plot / n_cols))
#     fig, axes = plt.subplots(n_rows, n_cols, figsize=(5 * n_cols, 5 * n_rows))
#     axes = axes.flatten()
#     for index, exp_name in enumerate(exp_list[:img_to_plot]):
#         if index ==0:
#             continue  # somehow index 0 is truncated. Unknown why, will need to rerun the overview images script.
#         ax = axes[index]
#         base_name = os.path.basename(exp_name)
#         # load rgb image
#         rgb_image = plt.imread(f"{rgb_images_folder}/{base_name}_{rgb_modifier}.png")
#         adata_exp = adata_mod[adata_mod.obs['dataset_name'] == base_name].copy()
#         if white_background:
#             rgb_image = np.zeros_like(rgb_image)
#         ax.imshow(rgb_image)
#         ax.grid(False)  # remove seaborn grid
#         # adjust point size based on the size of the rgb image. 8000x8000=0.5, adjust the rest accordingly
#         point_size = 0.5 * ((8000*8000) / (rgb_image.shape[0]*rgb_image.shape[1]))
#         sns.scatterplot(data=adata_exp.obs, x='X', y='Y', hue=metacluster_key, ax=ax, s=point_size, alpha=1,  #0.7,
#                         palette=adata_exp.uns[f"{metacluster_key}_colors"], hue_order=adata_exp.uns[f"{metacluster_key}_order"])
#         ax.set_title(f"{base_name}")
#         # remove legend
#         ax.get_legend().remove()
#         ax.set_xticks([])  # remove x and y-axis and ticks
#         ax.set_yticks([])
#         ax.set_xlabel('')  # remove axis labels
#         ax.set_ylabel('')
#         if white_background:
#             ax.spines['top'].set_visible(False)
#             ax.spines['right'].set_visible(False)
#             ax.spines['bottom'].set_visible(False)
#             ax.spines['left'].set_visible(False)
#     # legend outside plot and dot size bigger
#     plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', markerscale=5)
#     if white_background:
#         plt.savefig(f"{fig_dir}/{img_to_plot}grid_{cn_col}_{img_modifier}_{metacluster_key}.pdf",
#                     bbox_inches='tight', dpi=500)
#     else:
#         plt.savefig(f"{fig_dir}/{img_to_plot}grid_{cn_col}_{img_modifier}_rgb_{metacluster_key}_overlay.pdf",
#                     bbox_inches='tight', dpi=500)
#     plt.close()
#
# # %% cluster the neighbourhoods into metaclusters
# tissue_zones = ['B-cell follicle', 'T-cell zone', 'Medulla-Interfollicular zone-SCS']
# nbh_zone_dict = {'B-cell follicle': [5, 11, 12, 15], 'T-cell zone': [1, 4, 6,9,13,19, 10],  # moved 10 to T cell zone
#                      'Medulla-Interfollicular zone-SCS': [0,2,3,7,8,14,16,17, 18]}
# version_key = 'v4_8'  # version of the metacluster assignment
# metacluster_key_orig = f'Metacluster {version_key}'
# for zone in tissue_zones:
#     adata_mod.obs.loc[adata_mod.obs[cn_col].isin(nbh_zone_dict[zone]), metacluster_key_orig] = zone
# # add metacluster colors to adata
# adata_mod.uns[f'{metacluster_key_orig}_colors'] = [mcolors.rgb2hex(color) for color in sns.color_palette("tab10", len(tissue_zones))]
# adata_mod.uns[f'{metacluster_key_orig}_order'] = tissue_zones
#
# #%% save adata_mod and df_nbh_ct_log2fc
# adata_mod.write(f"{processed_spatial_data_dir}/sa_{modifier_base}_inclregions_{version_key}.h5ad")
# df_nbh_ct_log2fc.to_csv(f"{processed_spatial_data_dir}/df_nbh_ct_log2fc_{version_key}.csv", index=False)
# #%% alt: read from here. comment out lines 69-224 and uncomment this block (lines 230-237):
version_key = 'v4_8'
metacluster_key_orig = f'Metacluster {version_key}'
nbh_zone_dict = {'B-cell follicle': [5, 11, 12, 15], 'T-cell zone': [1, 4, 6,9,13,19, 10],  # moved 10 to T cell zone
                      'Medulla-Interfollicular zone-SCS': [0,2,3,7,8,14,16,17, 18]}
tissue_zones = ['B-cell follicle', 'T-cell zone', 'Medulla-Interfollicular zone-SCS']
adata_mod = ad.read_h5ad(f"{processed_spatial_data_dir}/sa_{modifier_base}_inclregions_{version_key}.h5ad")
df_nbh_ct_log2fc = pd.read_csv(f"{processed_spatial_data_dir}/df_nbh_ct_log2fc_{version_key}.csv")
cd45_cells = set(adata_mod.obs['Cell type'].unique()) - {'Endothelial', 'FRC', 'FDC', 'Unclassified'}

#%% filter out small regions inside the big ones

# majority_thresholds = [0.5] #[0.3, 0.5, 0.6, 0.7]
#%% revision: run with max_iterations set to 25 to achieve real stability. Update metacluster_key to  reflect the
# change and rerun downstream visualisations with it.
max_iterations = 25
majority_threshold = 0.5
metacluster_key = f"{metacluster_key_orig}_filtered_{majority_threshold}_maxit_{max_iterations}_{version_key}"
if metacluster_key not in adata_mod.obs.columns: # skip if possible, since it takes about 15 min.
    print(metacluster_key)
    spatial_metacluster_filtering(adata_mod,  metacluster_key=metacluster_key_orig,
                                      spatial_key='position', n_neighbors=10,
                                      majority_threshold=majority_threshold, max_iterations=max_iterations,
                                  key_added=metacluster_key)
# relabel key_added to metacluster_key (since I acidentally added the old one)
# adata_mod.obs[metacluster_key] = adata_mod.obs['Metacluster v4_8_filtered_0.5_v4_8']
# adata_mod.obs.drop(columns=['Metacluster v4_8_filtered_0.5_v4_8'], inplace=True) # remove old one
#  seaborn clustermap with annotation of metacluster assignment
if visualise:
    # --- Assign colors based on cn_col (neighbourhoods) and metacluster
    neighbourhoods = df_nbh_ct_log2fc.index

    # Outer layer: metacluster colors
    # Map from neighbourhood label to metacluster
    neigh_to_metacluster = {
        neigh: zone for zone, neigh_list in nbh_zone_dict.items() for neigh in neigh_list
    }
    metacluster_colors = adata_mod.uns[f'{metacluster_key_orig}_colors']  # length 4
    cluster_colors = adata_mod.uns[f'{cn_col}_colors']  # length 10
    # map metacluster colours to the neighbourhoods
    metacluster_name_to_idx = {name: idx for idx, name in enumerate(tissue_zones)}
    metacluster_colors_dict = {
        neigh: metacluster_colors[metacluster_name_to_idx[neigh_to_metacluster[neigh]]]
        for neigh in neighbourhoods
    }
    metacluster_colors_mapped = [metacluster_colors[metacluster_name_to_idx[neigh_to_metacluster[neigh]]]
        for neigh in neighbourhoods
    ]
    # Combine into a DataFrame for row_colors
    row_colors = pd.DataFrame({
        'Tissue zone': metacluster_colors_mapped,
        'Neighbourhood': cluster_colors
    }, index=df_nbh_ct_log2fc.index)

    legend_handles = [Patch(facecolor=color, label=label) for label, color in metacluster_colors_dict.items()]
    # --- Plot
    fig, ax = plt.subplots(1, 1, figsize=(15, 5))
    clus_map = sns.clustermap(df_nbh_ct_log2fc,
                              cmap='bwr',
                              figsize=(30, round(n/2)),
                              row_cluster=True,
                              col_cluster=True,
                              center=0, vmin=-3, vmax=3,
                              row_colors=row_colors,
                              cbar_kws={'label': 'log2 fold change'})

    # Set labels and title
    clus_map.ax_cbar.set_ylabel('log2 fold change', size=16)
    for ax in [clus_map.ax_row_dendrogram, clus_map.ax_col_dendrogram, clus_map.ax_heatmap]:
        ax.set_xticklabels(ax.get_xticklabels(), fontsize=16)
        ax.set_yticklabels(ax.get_yticklabels(), fontsize=16)
    clus_map.ax_cbar.tick_params(labelsize=16)

    # Add the legend to the right of the clustermap
    clus_map.ax_heatmap.legend(
        handles=legend_handles,
        title='Tissue zone',
        loc='center left',
        bbox_to_anchor=(1.05, 0.5),
        frameon=False
    )
    plt.suptitle(f"Cellular neighbourhood enrichment log2 fold change for all lymph nodes, tissue zones annotated",
                 fontsize=30)
    # Adjust layout to accommodate the legend
    plt.tight_layout(rect=[0, 0, 0.85, 1])
    plt.savefig(f"{fig_dir}/cn_exp_heatmap_all_lymph_nodes_{cn_col}_{cell_type_col}_{metacluster_key}_"
                f"{img_modifier}_log2fc_tissueregions_{version_key}.pdf", bbox_inches='tight', dpi=300)
    plt.close()

# %% plot and save the scatterplot with the 3 meta-neighbourhoods with optional overlay across the RGB image
if visualise:
    n_cols = 3
    n_rows = int(np.ceil(img_to_plot / n_cols))
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(5 * n_cols, 5 * n_rows))
    axes = axes.flatten()
    for index, exp_name in enumerate(exp_list[:img_to_plot]):
        if index == 0:
            continue  # somehow index 0 is truncated. Unknown why, will need to rerun the overview images script.
        ax = axes[index]
        base_name = os.path.basename(exp_name)
        # load rgb image
        rgb_image = plt.imread(f"{rgb_images_folder}/{base_name}_{rgb_modifier}.png")
        adata_exp = adata_mod[adata_mod.obs['dataset_name'] == base_name].copy()
        if white_background:
            rgb_image = np.zeros_like(rgb_image)
        ax.imshow(rgb_image)
        ax.grid(False)  # remove seaborn grid
        # adjust point size based on the size of the rgb image. 8000x8000=0.5, adjust the rest accordingly
        point_size = 0.5 * ((8000*8000) / (rgb_image.shape[0]*rgb_image.shape[1]))
        sns.scatterplot(data=adata_exp.obs, x='X', y='Y', hue=metacluster_key, ax=ax, s=point_size, alpha=0.7,
                        palette=adata_exp.uns[f'{metacluster_key_orig}_colors'], hue_order=adata_exp.uns[f'{metacluster_key_orig}_order'])
        ax.set_title(f"{base_name}")
        # remove legend
        ax.get_legend().remove()
        ax.set_xticks([])  # remove x and y-axis and ticks
        ax.set_yticks([])
        ax.set_xlabel('')  # remove axis labels
        ax.set_ylabel('')
        if white_background:
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.spines['bottom'].set_visible(False)
            ax.spines['left'].set_visible(False)
    # legend outside plot and dot size bigger
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', markerscale=5)
    if white_background:
        plt.savefig(f"{fig_dir}/{img_to_plot}grid_{cn_col}_{img_modifier}_{metacluster_key}.pdf",
                    bbox_inches='tight', dpi=300)
    else:
        plt.savefig(f"{fig_dir}/{img_to_plot}grid_{cn_col}_{img_modifier}_rgb_{metacluster_key}_overlay.pdf",
                    bbox_inches='tight', dpi=300)
    plt.close()
    # %% count abundances of cells of interest in different metadata columns
    # dict with metacluster as key and neutrophil count as value
    for cell_of_interest in cells_of_interest:
        coi_nbh_df = pd.DataFrame(columns=['exp_name', 'Organ', 'Treatment', 'mouse_id'])
        for index_exp, exp_name in enumerate(tqdm(exp_list)):
            adata_exp = adata_mod[adata_mod.obs['dataset_name'] == os.path.basename(exp_name)].copy()
            if adata_exp is None or len(adata_exp.obs) == 0:
                print(f"Skipping empty adata for {exp_name}")  # correctly skips the removed mouse HOE-2257.
                continue
            organ = adata_exp.obs['Organ'].values[0]
            treatment = adata_exp.obs['Treatment'].values[0]
            mouse_id = adata_exp.obs['mouse_id'].values[0]
            meta_cl_abundance = {'exp_name': exp_name, 'Organ': organ, 'Treatment': treatment, 'mouse_id': mouse_id}
            # sum of all neutrophils
            all_neut = (adata_exp.obs[cell_type_col] == cell_of_interest).sum()
            # sum all cells in adata_exp
            all_cells = len(adata_exp.obs)
            # find rows where neutrophil is true and medulla (metacluster) is true
            for metacluster in tissue_zones:
                neut_meta = ((adata_exp.obs[cell_type_col] == cell_of_interest) &
                             (adata_exp.obs[metacluster_key] == metacluster)).sum()
                # normalise by total neutrophils
                neut_meta_norm = neut_meta / all_neut
                neut_all_cell_norm = neut_meta / all_cells
                # normalise by total cells in that metacluster
                all_cells_meta = len(adata_exp.obs[adata_exp.obs[metacluster_key] == metacluster])
                neut_meta_norm_cells = neut_meta / all_cells_meta
                # normalise by cd45+ cells in that metacluster.
                cd45_cells_meta = len(adata_exp.obs[(adata_exp.obs[metacluster_key] == metacluster) &
                                                    (adata_exp.obs['Cell type'].isin(cd45_cells))])
                neut_meta_norm_cd45 = neut_meta / cd45_cells_meta
                meta_cl_abundance[f'{metacluster}_percent_n'] = neut_meta_norm * 100  # norm by total neutrophils
                meta_cl_abundance[f'{metacluster}_percent_c'] = neut_all_cell_norm * 100  # norm by total cells
                meta_cl_abundance[f'{metacluster}_percent_z'] = neut_meta_norm_cells * 100
                meta_cl_abundance[f'{metacluster}_percent_cd45_z'] = neut_meta_norm_cd45 * 100  # norm by cd45+ cells in that metacluster
                meta_cl_abundance[f'{metacluster}_count'] = neut_meta
            coi_nbh_df = pd.concat([coi_nbh_df, pd.DataFrame(meta_cl_abundance, index=[0])], ignore_index=True)
        # %% drop double rows from coi_nbh_df (in case I restarted the loop)
        coi_nbh_df.drop_duplicates(inplace=True)
        # check if df now has the same length as exp_list
        assert len(coi_nbh_df) == len(exp_list)
        # %%save the dataframe
        coi_nbh_df.to_csv(f"{processed_spatial_data_dir}/{cell_of_interest}_abundance_3regions_"
                          f"{cn_col}_{metacluster_key}_filtered.csv",
                          index=False)
        del coi_nbh_df

# %%save the adata
adata_mod.write(f"{processed_spatial_data_dir}/sa_{modifier_base}_inclregions_{img_modifier}_{metacluster_key}.h5ad")
print('Saved adata_mod with metacluster assignments and filtering applied.')

#%% plot subset as a smaller clustermap to use as an example in the workflow (fig 4A)
# --- Define subset of neighbourhoods and cell types to include ---
selected_neighbourhoods = [0, 2, 3, 8, 10]
selected_cells = ['B cell', 'CD4 T cell', 'CD8 T cell', 'Neutrophil', 'CD169+ M']

# --- Subset data ---
df_subset = df_nbh_ct_log2fc.loc[selected_neighbourhoods, selected_cells]

# --- Tissue zone color mapping ---
# Assign each neighbourhood to its metacluster name
selected_zone_names = [neigh_to_metacluster[neigh] for neigh in selected_neighbourhoods]
# Map zone name to index for metacluster color lookup
metacluster_name_to_idx = {name: idx for idx, name in enumerate(tissue_zones)}
tissue_zone_colors = [metacluster_colors[metacluster_name_to_idx[zone]] for zone in selected_zone_names]

# --- Neighbourhood-specific color mapping ---
cluster_colors_subset = [cluster_colors[neigh] for neigh in selected_neighbourhoods]

# --- Row colors DataFrame ---
row_colors = pd.DataFrame({
    'Tissue zone': tissue_zone_colors,
    'Neighbourhood': cluster_colors_subset
}, index=selected_neighbourhoods)

# --- Legend: one handle per tissue zone (not per neighbourhood) ---
# To avoid repeated labels, build a mapping from zone name to color
unique_zone_color_dict = {
    zone: metacluster_colors[metacluster_name_to_idx[zone]]
    for zone in set(selected_zone_names)
}
legend_handles = [
    Patch(facecolor=color, label=zone) for zone, color in unique_zone_color_dict.items()
]

# --- Plot ---
fig, ax = plt.subplots(1, 1, figsize=(15, 5))
clus_map = sns.clustermap(df_subset,
                          cmap='bwr',
                          figsize=(len(selected_cells)*2, round(len(selected_neighbourhoods))),
                          row_cluster=True,
                          col_cluster=True,
                          center=0, vmin=-3, vmax=3,
                          row_colors=row_colors,
                          cbar_kws={'label': 'log2 fold change'})

# Clean ticks
for ax in [clus_map.ax_row_dendrogram, clus_map.ax_col_dendrogram, clus_map.ax_heatmap]:
    ax.set_xticklabels(ax.get_xticklabels())
    ax.set_yticklabels(ax.get_yticklabels())

plt.tight_layout(rect=[0, 0, 0.85, 1])
plt.savefig(f"{fig_dir}/cn_exp_heatmap_subset_{cn_col}_{cell_type_col}_"
            f"{img_modifier}_log2fc_tissueregions_{metacluster_key}.pdf", bbox_inches='tight', dpi=300)
plt.close()
