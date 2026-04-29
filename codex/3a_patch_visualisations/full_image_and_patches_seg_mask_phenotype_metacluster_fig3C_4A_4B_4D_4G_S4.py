#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
# '''=================================================
# @Project: ACT_melanoma_neutrophil (data: act_neutrophil_data_repository)
# @Author : Gemma van der Voort
# @Desc   : script to create nice snapshots for the paper. produce a full image with zoom window, the zoomed-in portion
# @Desc updated: Patch visualizations with segmentation and phenotyping overlays.
# at high quality, a segmentation mask overlaying the zoom,and the phenotyping information overlaying the zoom.
# '''=================================================
# %% imports
import sys
from pathlib import Path

import numpy as np
import anndata as ad
import skimage
import matplotlib.pyplot as plt

import matplotlib.colors as mcolors
import os
import warnings

PROJECT_ROOT = Path(__file__).resolve()
while PROJECT_ROOT.name != "ACT_melanoma_neutrophil" and PROJECT_ROOT.parent != PROJECT_ROOT:
    PROJECT_ROOT = PROJECT_ROOT.parent
if PROJECT_ROOT.name != "ACT_melanoma_neutrophil":
    raise RuntimeError("ACT_melanoma_neutrophil root not found")
sys.path.append(str(PROJECT_ROOT / "helper_files"))

from paths_parameters import fig_dir, modifier_base, parent_dir, processed_spatial_data_dir, cell_colours_hex, data_repo_path


# Helper function to find nearest nonzero cell ID in a local window
def find_nearest_nonzero(matrix, x, y, window_size=3):
    """
    Find the nearest nonzero value in a local window around the given coordinates.
    :param matrix:
    :param x:
    :param y:
    :param window_size:
    :return:
    """
    half_window = window_size // 2
    x_min, x_max = max(0, x - half_window), min(matrix.shape[1], x + half_window + 1)
    y_min, y_max = max(0, y - half_window), min(matrix.shape[0], y + half_window + 1)

    window = matrix[y_min:y_max, x_min:x_max]
    nonzero_values = window[window > 0]
    return nonzero_values[0] if len(nonzero_values) > 0 else 0


# %%import rbg images (needed if visualise=True)
rgb_images_folder = (f'{data_repo_path}/figures/codex/overview_intensity_images/full_ln_architecture/'
                     'channels_of_interest_highQ')  # input figure_annotation_dir to images
rgb_modifier = 'DAPI_B220_aSMA_CD3_ERTR7'
warnings.filterwarnings("ignore")
cell_type_col = 'Cell type'
cell_type_col_file_name = 'cell_type'  # for use in file names, no spaces or caps. also poss shorter
version_key = 'v4_8'  # version of the metacluster key, large scale classifier.
metacluster_key_orig = f'Metacluster {version_key}'
cn_col = 'CN_k50_n20'  # column name for the CN in adata

majority_threshold = 0.5  # threshold for filtering metaclusters
metacluster_col_name = f"{metacluster_key_orig}_filtered_{majority_threshold}_{version_key}"
img_modifier = 'v4_corrected'
modifier = f'{img_modifier}_{version_key}_no_APC'

fig_dir = f'{fig_dir}/patch_visualisations'
os.makedirs(fig_dir, exist_ok=True)
# data we need:
# - adata with cell type column
# - rgb images
exp_list = [os.path.join(parent_dir, f) for f in os.listdir(parent_dir) if os.path.isdir(os.path.join(parent_dir, f))]
adata = ad.read_h5ad(f"{processed_spatial_data_dir}/sa_{modifier_base}_inclregions_{img_modifier}_{version_key}.h5ad")

#%% # for data-reg_5x2_HOE-2367_inLNr, the zooms window are:
zoom_windows = [[(8200, 9000, 1100, 1700), (2200, 3000, 900, 1500), (3200, 4000, 2200, 2800)],
                [(3200, 4000, 4100, 4700), (2800,3600, 2000,2600), (2400,3200,3100,3700)]]
# image_path = f'{data_repo_path}/raw/codex/data-reg_5x2_HOE-2367_inLNr'  # d14, inLNr.
# zoom_windows_1_image = zoom_windows[0]

# ALTERNATIVE: D14, inLNl
image_path = f'{data_repo_path}/raw/codex/data-reg_3x3_HOE-2367_inLNl'  # d14, inLNl.
# for data-reg_3x3_HOE-2367_inLNl, the zoom window is: (for fig 4D etc)
zoom_windows_1_image = zoom_windows[1]
# for fig 4F: CHECK MODIFIER!
base_name = os.path.basename(image_path)
# load rgb & subset adata
rgb_image = plt.imread(f"{rgb_images_folder}/{base_name}_{rgb_modifier}.png")
adata_exp = adata[adata.obs['dataset_name'] == base_name].copy()

# plot to explore
_, ax = plt.subplots(1, 1)
ax.imshow(rgb_image)  # imshow with 0,0 in the top left corner
ax.grid(True)  # show seaborn grid for finding a good zoom window
plt.savefig(f'{fig_dir}/full_image_{base_name}.png', dpi=300)
plt.close()
# load segmentation mask
seg_mask_path = f'{parent_dir}/{base_name}/mask/R1_membrane.tiff'
seg_mask = plt.imread(seg_mask_path)
# %%draw rectangles on full image
_, ax = plt.subplots(1, 1)
ax.imshow(rgb_image)
ax.grid(False)  # now remove grid
# remove axis labels and ticks
ax.set_xticks([])
ax.set_yticks([])
for zoom_window in zoom_windows_1_image:
    x1, x2, y1, y2 = zoom_window
    rect = plt.Rectangle((x1, y1), x2 - x1, y2 - y1, edgecolor='w', facecolor='none')
    ax.add_patch(rect)
plt.savefig(f'{fig_dir}/full_image_with_zoom_{base_name}_{len(zoom_windows_1_image)}zooms.png', dpi=300,
            bbox_inches='tight', pad_inches=0)  # removes white space around the image
plt.close()
# %%zoom proportions width 80, height 60 looks OK in fig outline.
# zoom in
for zoom_index, zoom_window in enumerate(zoom_windows_1_image):
    x1, x2, y1, y2 = zoom_window
    _, ax = plt.subplots(1, 1)
    ax.imshow(rgb_image[y1:y2, x1:x2])
    # remove axis labels and ticks
    ax.set_xticks([])
    ax.set_yticks([])
    plt.savefig(f'{fig_dir}/zoomed_image_{base_name}_{zoom_index}.png', dpi=300,
                bbox_inches='tight', pad_inches=0)  # removes white space around the image
    plt.close()
    # %%plot the zoomed image with the segmentation mask overlay. only plot the outlines of the mask.
    rgb_zoom = rgb_image[y1:y2, x1:x2]  # zoomed rgb image
    seg_mask_zoom = seg_mask[y1:y2, x1:x2]  # zoomed seg mask
    seg_mask_outline = skimage.segmentation.find_boundaries(seg_mask_zoom, mode='thick')  # thick outlines
    overlay = rgb_zoom.copy()
    # Apply red color to boundary pixels only in RGB channels (ignore alpha channel)
    overlay[seg_mask_outline, :3] = [1, 1, 1]  # Set boundaries to red (R=1, G=0, B=0), keeping alpha unchanged
    # plot figure
    _, ax = plt.subplots(1, 1)
    ax.imshow(overlay)
    # turn of axis labels and ticks
    ax.set_xticks([])
    ax.set_yticks([])
    # remove the white background around the image
    plt.savefig(f'{fig_dir}/zoomed_image_seg_mask_overlay_{base_name}_{modifier}_zoom_{zoom_index}.png', dpi=300,
                bbox_inches='tight', pad_inches=0)  # removes white space around the image
    plt.close()


    # %% now add the phenotype information overlaying the zoom. Don't do this for the full image, as it is too intensive.
    cell_colours_rgb = {cell_type: mcolors.hex2color(cell_colours_hex[cell_type]) for cell_type in cell_colours_hex}
    # prepare for plotting the phenotype overlay
    labeled_mask_zoom = skimage.measure.label(seg_mask_zoom)  # Each cell gets a unique ID
    # num_cells_zoom = labeled_mask_zoom.max()  # Number of detected cells, should be => than rows in adata_exp subset.
    # Generate a colormap for unique cell types. Replace with new colourmap from umap once designed.
    unique_cell_types = adata_exp.obs[cell_type_col].unique()
    num_cell_types = len(unique_cell_types)
    # Create an RGB version of the mask
    seg_colored_zoom = np.zeros((*seg_mask_zoom.shape, 3))  # Shape: (H, W, 3)
    # subset df to only contain the zoomed in cells
    adata_exp_zoom = adata_exp[(adata_exp.obs['X'] >= x1) & (adata_exp.obs['X'] <= x2) &
                                (adata_exp.obs['Y'] >= y1) & (adata_exp.obs['Y'] <= y2)].copy()
    # map labeled regions to cell types
    for index, row in adata_exp_zoom.obs.iterrows():
        cell_x, cell_y = int(row["X"]), int(row["Y"])
        cell_type = row[cell_type_col]
        # Clip coordinates to ensure they are within valid bounds
        cell_x_clipped = np.clip(cell_x - x1, 0, labeled_mask_zoom.shape[1] - 1)
        cell_y_clipped = np.clip(cell_y - y1, 0, labeled_mask_zoom.shape[0] - 1)
        # Get initial cell ID at centroid location
        cell_id = labeled_mask_zoom[cell_y_clipped, cell_x_clipped]
        # If cell_id is 0, search for nearest labeled pixel
        if cell_id == 0:
            cell_id = find_nearest_nonzero(labeled_mask_zoom, cell_x_clipped, cell_y_clipped, window_size=3)
        # Assign color only if the cell_id is valid
        if cell_id > 0:
            seg_colored_zoom[labeled_mask_zoom == cell_id] = cell_colours_rgb[cell_type]
        else:
            print(f"Warning: No valid cell_id found for cell {index} at ({cell_x}, {cell_y})")

    # Compute segmentation boundaries
    seg_boundaries_zoom = skimage.segmentation.find_boundaries(labeled_mask_zoom, mode='thick')
    seg_colored_zoom[seg_boundaries_zoom] = [0, 0, 0]  # Overlay boundaries in black
    # plot figure
    plt.figure()
    plt.imshow(seg_colored_zoom)
    plt.axis("off")
    plt.xticks([])
    plt.yticks([])
    plt.savefig(f'{fig_dir}/zoomed_image_{cell_type_col_file_name}_overlay_zoom_{zoom_index}_{base_name}_'
                f'{version_key}_{modifier}.png', dpi=300,
                bbox_inches='tight', pad_inches=0)  # removes white space around the image
    plt.close()
    #%% do the same but for the metaclusters (with associated colours in adata_exp.uns)

    # load metacluster colours
    metacluster_colors = adata_exp.uns[f'{metacluster_key_orig}_colors']
    # convert metacluster colors from hex to rgb in 0-1 range
    metacluster_colors_rgb = np.array([np.array(mcolors.hex2color(metacluster_colors[i])) for i, metacluster in enumerate(metacluster_colors)])
    metacluster_order = adata_exp.uns[f'{metacluster_key_orig}_order']
    # both are arrays of length 3. make a dict of metacluster to colour
    metacluster_color_dict = {metacluster: metacluster_colors_rgb[i] for i, metacluster in enumerate(metacluster_order)}
    # Create an RGB version of the mask
    seg_colored_zoom = np.zeros((*seg_mask_zoom.shape, 3))  # Shape: (H, W, 3)
    # map labeled regions to metacluster id
    for index, row in adata_exp_zoom.obs.iterrows():
        cell_x, cell_y = int(row["X"]), int(row["Y"])
        metacluster = row[metacluster_col_name]  # print metacluster_key_orig for unfiltered
        # Clip coordinates to ensure they are within valid bounds
        cell_x_clipped = np.clip(cell_x - x1, 0, labeled_mask_zoom.shape[1] - 1)
        cell_y_clipped = np.clip(cell_y - y1, 0, labeled_mask_zoom.shape[0] - 1)
        # Get initial cell ID at centroid location
        cell_id = labeled_mask_zoom[cell_y_clipped, cell_x_clipped]
        # If cell_id is 0, search for the nearest labeled pixel
        if cell_id == 0:
            cell_id = find_nearest_nonzero(labeled_mask_zoom, cell_x_clipped, cell_y_clipped, window_size=3)
        # Assign color only if the cell_id is valid
        if cell_id > 0:
            seg_colored_zoom[labeled_mask_zoom == cell_id] = metacluster_color_dict[metacluster]
        else:
            print(f"Warning: No valid cell_id found for cell {index} at ({cell_x}, {cell_y})")
    # Compute segmentation boundaries
    seg_boundaries_zoom = skimage.segmentation.find_boundaries(labeled_mask_zoom, mode='thick')
    seg_colored_zoom[seg_boundaries_zoom] = [0, 0, 0]  # Overlay boundaries in black
    # plot figure
    plt.figure()
    plt.imshow(seg_colored_zoom)
    plt.axis("off")
    plt.xticks([])
    plt.yticks([])
    plt.savefig(f'{fig_dir}/zoomed_image_{metacluster_col_name}_overlay_zoom_{zoom_index}_{base_name}_{version_key}_{modifier}.png', dpi=300,
                bbox_inches='tight', pad_inches=0)  # removes white space around the image
    plt.close()
    #%% small change of the above with the unfiltered metacluster assignments, for the workflow image
    metacluster_col_name_unfiltered = f"{metacluster_key_orig}"
    # load metacluster colours
    metacluster_colors = adata_exp.uns[f'{metacluster_key_orig}_colors']
    # convert metacluster colors from hex to rgb in 0-1 range
    metacluster_colors_rgb = np.array([np.array(mcolors.hex2color(metacluster_colors[i])) for i, metacluster in enumerate(metacluster_colors)])
    metacluster_order = adata_exp.uns[f'{metacluster_key_orig}_order']
    # both are arrays of length 3. make a dict of metacluster to colour
    metacluster_color_dict = {metacluster: metacluster_colors_rgb[i] for i, metacluster in enumerate(metacluster_order)}
    # Create an RGB version of the mask
    seg_colored_zoom = np.zeros((*seg_mask_zoom.shape, 3))  # Shape: (H, W, 3)
    # map labeled regions to metacluster id
    for index, row in adata_exp_zoom.obs.iterrows():
        cell_x, cell_y = int(row["X"]), int(row["Y"])
        metacluster = row[metacluster_col_name_unfiltered]  # print metacluster_key_orig for unfiltered
        # Clip coordinates to ensure they are within valid bounds
        cell_x_clipped = np.clip(cell_x - x1, 0, labeled_mask_zoom.shape[1] - 1)
        cell_y_clipped = np.clip(cell_y - y1, 0, labeled_mask_zoom.shape[0] - 1)
        # Get initial cell ID at centroid location
        cell_id = labeled_mask_zoom[cell_y_clipped, cell_x_clipped]
        # If cell_id is 0, search for nearest labeled pixel
        if cell_id == 0:
            cell_id = find_nearest_nonzero(labeled_mask_zoom, cell_x_clipped, cell_y_clipped, window_size=3)
        # Assign color only if the cell_id is valid
        if cell_id > 0:
            seg_colored_zoom[labeled_mask_zoom == cell_id] = metacluster_color_dict[metacluster]
        else:
            print(f"Warning: No valid cell_id found for cell {index} at ({cell_x}, {cell_y})")
    # Compute segmentation boundaries
    seg_boundaries_zoom = skimage.segmentation.find_boundaries(labeled_mask_zoom, mode='thick')
    seg_colored_zoom[seg_boundaries_zoom] = [0, 0, 0]  # Overlay boundaries in black
    # plot figure
    plt.figure()
    plt.imshow(seg_colored_zoom)
    plt.axis("off")
    plt.xticks([])
    plt.yticks([])
    plt.savefig(f'{fig_dir}/zoomed_image_{metacluster_col_name_unfiltered}_overlay_zoom_{zoom_index}_{base_name}_{version_key}_{modifier}.png', dpi=300,
                bbox_inches='tight', pad_inches=0)  # removes white space around the image
    plt.close()
    # %%adapt the above a bit: in addition to the metacluster colours, also plot the neutrophils
    neutrophil_rgb = np.array([0.83921569, 0.15294118, 0.15686275])  # RGB for Neutrophil, red
    # Create an RGB version of the mask
    seg_colored_zoom = np.zeros((*seg_mask_zoom.shape, 3))  # Shape: (H, W, 3)
    # map labeled regions to metacluster id
    for index, row in adata_exp_zoom.obs.iterrows():
        cell_x, cell_y = int(row["X"]), int(row["Y"])
        metacluster = row[metacluster_col_name]
        # Clip coordinates to ensure they are within valid bounds
        cell_x_clipped = np.clip(cell_x - x1, 0, labeled_mask_zoom.shape[1] - 1)
        cell_y_clipped = np.clip(cell_y - y1, 0, labeled_mask_zoom.shape[0] - 1)
        # Get initial cell ID at centroid location
        cell_id = labeled_mask_zoom[cell_y_clipped, cell_x_clipped]
        # If cell_id is 0, search for nearest labeled pixel
        if cell_id == 0:
            cell_id = find_nearest_nonzero(labeled_mask_zoom, cell_x_clipped, cell_y_clipped, window_size=3)
        # Assign color only if the cell_id is valid
        if cell_id > 0:
            seg_colored_zoom[labeled_mask_zoom == cell_id] = metacluster_color_dict[metacluster]
        else:
            print(f"Warning: No valid cell_id found for cell {index} at ({cell_x}, {cell_y})")
    # plot neutrophils from Cell type column in red
    neut = adata_exp_zoom[adata_exp_zoom.obs[cell_type_col] == 'Neutrophil']
    for index, row in neut.obs.iterrows():
        cell_x, cell_y = int(row["X"]), int(row["Y"])
        # Clip coordinates to ensure they are within valid bounds
        cell_x_clipped = np.clip(cell_x - x1, 0, labeled_mask_zoom.shape[1] - 1)
        cell_y_clipped = np.clip(cell_y - y1, 0, labeled_mask_zoom.shape[0] - 1)
        # Get initial cell ID at centroid location
        cell_id = labeled_mask_zoom[cell_y_clipped, cell_x_clipped]
        # If cell_id is 0, search for nearest labeled pixel
        if cell_id == 0:
            cell_id = find_nearest_nonzero(labeled_mask_zoom, cell_x_clipped, cell_y_clipped, window_size=3)
        # Assign color only if the cell_id is valid
        if cell_id > 0:
            seg_colored_zoom[labeled_mask_zoom == cell_id] = neutrophil_rgb
        else:
            print(f"Warning: No valid cell_id found for cell {index} at ({cell_x}, {cell_y})")
    # Compute segmentation boundaries
    seg_boundaries_zoom = skimage.segmentation.find_boundaries(labeled_mask_zoom, mode='thick')
    seg_colored_zoom[seg_boundaries_zoom] = [0, 0, 0]  # Overlay boundaries in black
    # plot figure
    plt.figure()
    plt.imshow(seg_colored_zoom)
    plt.axis("off")
    plt.xticks([])
    plt.yticks([])
    plt.savefig(f'{fig_dir}/zoomed_image_{cell_type_col_file_name}_{metacluster_col_name}_overlay_zoom_'
                f'{zoom_index}_{base_name}_{version_key}_{modifier}_metacluster_and_neutrophils_new_colours.png', dpi=300,
                bbox_inches='tight', pad_inches=0)  # removes white space around the image
    plt.close()

    #%% repeat the above but for the CN_k50_n20 column
    cn_colours = adata_exp.uns[f'{cn_col}_colors']  # CN colours from adata
    cn_colours_rgb = np.array([np.array(mcolors.hex2color(cn_colours[i])) for i in range(len(cn_colours))])  # convert to rgb
    cn_order = adata_exp.uns[f'{cn_col}_order']  # CN order from adata
    # Create an RGB version of the mask
    seg_colored_zoom = np.zeros((*seg_mask_zoom.shape, 3))  # Shape: (H, W, 3)
    # map labeled regions to CN id
    for index, row in adata_exp_zoom.obs.iterrows():
        cell_x, cell_y = int(row["X"]), int(row["Y"])
        cn = row[cn_col]  # print metacluster_key_orig for unfiltered
        # Clip coordinates to ensure they are within valid bounds
        cell_x_clipped = np.clip(cell_x - x1, 0, labeled_mask_zoom.shape[1] - 1)
        cell_y_clipped = np.clip(cell_y - y1, 0, labeled_mask_zoom.shape[0] - 1)
        # Get initial cell ID at centroid location
        cell_id = labeled_mask_zoom[cell_y_clipped, cell_x_clipped]
        # If cell_id is 0, search for nearest labeled pixel
        if cell_id == 0:
            cell_id = find_nearest_nonzero(labeled_mask_zoom, cell_x_clipped, cell_y_clipped, window_size=3)
        # Assign color only if the cell_id is valid
        if cell_id > 0:
            seg_colored_zoom[labeled_mask_zoom == cell_id] = cn_colours_rgb[cn]
        else:
            print(f"Warning: No valid cell_id found for cell {index} at ({cell_x}, {cell_y})")
    # Compute segmentation boundaries
    seg_boundaries_zoom = skimage.segmentation.find_boundaries(labeled_mask_zoom, mode='thick')
    seg_colored_zoom[seg_boundaries_zoom] = [0, 0, 0]  # Overlay boundaries in black
    # plot figure
    plt.figure()
    plt.imshow(seg_colored_zoom)
    plt.axis("off")
    plt.xticks([])
    plt.yticks([])
    plt.savefig(f'{fig_dir}/zoomed_image_{cn_col}_overlay_zoom_{zoom_index}_{base_name}_{version_key}_{modifier}.png', dpi=300,
                bbox_inches='tight', pad_inches=0)  # removes white space around the image
    plt.close()
