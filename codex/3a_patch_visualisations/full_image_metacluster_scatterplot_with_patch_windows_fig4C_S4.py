#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
# '''=================================================
# @Project: ACT_melanoma_neutrophil (data: act_neutrophil_data_repository)
# @Author : Gemma van der Voort
# @Desc   : adapted from inpection_neut_infl_selected_tissues.py to plot the metacluster distribution with
# @Desc updated: Patch visualization with metacluster overlays.
# neutrophils and APC neutrophils in the lymph nodes for the two selected tissues.
# update 20260109: remove APC neutrophils, led to confusion from reviewers. Also add two more boxes for S5 (supplement
# on spatial things)
# '''=================================================
# %% imports
import sys
from pathlib import Path

import numpy as np
import anndata as ad
from tqdm import tqdm
import matplotlib.pyplot as plt
import seaborn as sns
import os
import warnings
import matplotlib.colors as mcolors

PROJECT_ROOT = Path(__file__).resolve()
while PROJECT_ROOT.name != "ACT_melanoma_neutrophil" and PROJECT_ROOT.parent != PROJECT_ROOT:
    PROJECT_ROOT = PROJECT_ROOT.parent
if PROJECT_ROOT.name != "ACT_melanoma_neutrophil":
    raise RuntimeError("ACT_melanoma_neutrophil root not found")
sys.path.append(str(PROJECT_ROOT / "helper_files"))

from paths_parameters import (fig_dir, modifier_base, processed_spatial_data_dir, cell_colours_hex, data_repo_path)

# import rbg images (needed if visualise=True)
rgb_images_folder = (f'{data_repo_path}/figures/codex/overview_intensity_images/full_ln_architecture/'
                     'channels_of_interest_highQ')  # input figure_annotation_dir to images
rgb_modifier = 'DAPI_B220_aSMA_CD3_ERTR7'
warnings.filterwarnings("ignore")
cell_type_col = 'Cell type'
cells_of_interest = ['Neutrophil']


fig_dir = f'{fig_dir}/patch_visualisations'
os.makedirs(fig_dir, exist_ok=True)
# get the colour of neutrophils from cell_colours_hex and transform to rgb in 0-1 range
neutrophil_rgb = np.array(mcolors.hex2color(cell_colours_hex['Neutrophil']))
# apc_neutrophil_rgb = np.array(mcolors.hex2color(cell_colours_hex['APC Neutrophil']))
#%%
img_modifier = 'v4_corrected'  # for the cn clustering
version_key = 'v4_8'  # version key for the metacluster annotations
metacluster_key_orig = f'Metacluster {version_key}'  # original metacluster key # for large regions
# metacluster_key_orig = f'CN_k50_n20'  # for 20 nbhs. update modifier!
majority_threshold = 0.5
modifier = '3_windows_3regions_no_cells'
adata = ad.read_h5ad(f"{processed_spatial_data_dir}/sa_{modifier_base}_inclregions_{img_modifier}_{version_key}.h5ad")
#%%
base_names = ['data-reg_5x2_HOE-2367_inLNr', 'data-reg_3x3_HOE-2367_inLNl']
zoom_windows = [[(8200, 9000, 1100, 1700), (2200, 3000, 900, 1500), (3200, 4000, 2200, 2800)],
                [(3200, 4000, 4100, 4700), (2800,3600, 2000,2600), (2400,3200,3100,3700)]]
metacluster_col_name = f"{metacluster_key_orig}_filtered_{majority_threshold}_{version_key}"
# metacluster_col_name = metacluster_key_orig  # for 20 nbhs. (CNs)
# %%make 2 plots of one image each.
# plotting decision logic:
white_background = False
black_background = True  # only one of these can be True
add_legend = False
#%%
for image_idx, base_name in tqdm(enumerate(base_names)):
    fig, ax = plt.subplots(1, 1)
    if black_background:
        fig.patch.set_facecolor('black')
        ax.set_facecolor('black')

    # load rgb
    rgb_image = plt.imread(f"{rgb_images_folder}/{base_name}_{rgb_modifier}.png")
    adata_exp = adata[adata.obs['dataset_name'] == base_name].copy()
    # remove the rgb image by setting it to 0 to plot the metacluster colours on a white background while
    # preserving the shape/orientation of the image
    if white_background:
        rgb_image = np.ones_like(rgb_image)
        zoom_window_color = 'black'
    elif black_background:  # remove image
        rgb_image = np.zeros_like(rgb_image)
        zoom_window_color = 'white'
    ax.imshow(rgb_image)
    ax.grid(False)  # remove seaborn grid
    # adjust point size based on the size of the rgb image. 8000x8000=0.5, adjust the rest accordingly
    point_size = 0.5 * ((8000*8000) / (rgb_image.shape[0]*rgb_image.shape[1]))
    # plot metacluster using metacluster colours
    sns.scatterplot(data=adata_exp.obs, x='X', y='Y', hue=metacluster_col_name, ax=ax, s=point_size, alpha=1,  #0.7,
                    palette=adata_exp.uns[f'{metacluster_key_orig}_colors'], hue_order=adata_exp.uns[f'{metacluster_key_orig}_order'])
    # remove legend
    ax.get_legend().remove()
    # plot neutrophils from Cell type column in red
    # neut = adata_exp[adata_exp.obs[cell_type_col] == 'Neutrophil']
    # ax.scatter(neut.obs['X'], neut.obs['Y'], s=point_size*2, c=neutrophil_rgb)
    # add zoom window in black
    zoom_windows_1_image = zoom_windows[image_idx]
    for zoom_window in zoom_windows_1_image:
        x1, x2, y1, y2 = zoom_window
        ax.add_patch(plt.Rectangle((x1, y1), x2-x1, y2-y1, fill=False, edgecolor=zoom_window_color, lw=1))

    ax.set_xticks([])  # remove x and y-axis and ticks
    ax.set_yticks([])
    ax.set_xlabel('')  # remove axis labels
    ax.set_ylabel('') # remove borders
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    if add_legend:
        # modify legend to include the scatterplot and the metacluster legend
        handles, labels = ax.get_legend_handles_labels()
        # add neutrophil legend
        handles.append(plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=neutrophil_rgb, markersize=1,
                                  label='Neutrophil'))
        labels.append('Neutrophil')
        plt.legend(handles=handles, labels=labels, bbox_to_anchor=(1.05, 1), loc='upper left', markerscale=5)
        # save legend separately as svg
        plt.savefig(f"{fig_dir}/{base_name}_{modifier}_neut_region_overlay_legend_{version_key}.png", bbox_inches='tight', dpi=300)
    if white_background:
        plt.savefig(f"{fig_dir}/{base_name}_{modifier}_neut_region_overlay_white_{version_key}.png",
                    bbox_inches='tight', pad_inches=0,
                    dpi=300)
    elif black_background:
        plt.savefig(f"{fig_dir}/{base_name}_{modifier}_neut_region_overlay_black_{version_key}.png",
                    bbox_inches='tight', pad_inches=0,
                    dpi=300)
    else:
        plt.savefig(f"{fig_dir}/{base_name}_{modifier}_neut_region_overlay_{version_key}.png",
                    bbox_inches='tight',  pad_inches=0, dpi=300)
    plt.close()
