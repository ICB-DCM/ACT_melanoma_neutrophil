#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
# '''=================================================
# @Project: ACT_melanoma_neutrophil (data: act_neutrophil_data_repository)
# @Author : Gemma van der Voort
# @Desc   : for the size matters fig, create one visualisation per mouse with the lymph nodes side by side in their
# @Desc updated: Side-by-side whole-LN visualizations per mouse.
# original size. so suptitle: mouse ID, treatment condition. then three subplots: inLNr, inLNl, and brLNr. pad the
# images with a black background.
# '''=================================================
# %%imports
import os
import sys
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import matplotlib.gridspec as gridspec # for subplots with correct (real) proportions
from matplotlib.patches import Rectangle
import numpy as np

PROJECT_ROOT = Path(__file__).resolve()
while PROJECT_ROOT.name != "ACT_melanoma_neutrophil" and PROJECT_ROOT.parent != PROJECT_ROOT:
    PROJECT_ROOT = PROJECT_ROOT.parent
if PROJECT_ROOT.name != "ACT_melanoma_neutrophil":
    raise RuntimeError("ACT_melanoma_neutrophil root not found")
sys.path.append(str(PROJECT_ROOT / "helper_files"))

from paths_parameters import data_repo_path, modifier

images_dir = f'{data_repo_path}/figures/codex/overview_intensity_images/full_ln_architecture/channels_of_interest_lowQ'
processed_data_dir = f'{data_repo_path}/processed/codex'
images_list = os.listdir(images_dir)
fig_dir = f'{data_repo_path}/figures/codex/whole_LN_visualisations/side_by_side'
os.makedirs(fig_dir, exist_ok=True)

# %%load adata for mapping
adata = ad.read_h5ad(f'{processed_data_dir}/phenotyping_{modifier}.h5ad')

# %%create mapping: adata.obs['mouse_ID'] and adata.obs['Treatment']. find 1st occurence of each mouse ID and save the
# treatment condition. save as a dictionary with mouse ID as key and treatment condition as value.
mouse_mapping = {}
for mouse_id in adata.obs['mouse_id'].unique():
    treatment = adata.obs.loc[adata.obs['mouse_id'] == mouse_id, 'Treatment'].values[0]
    mouse_mapping[mouse_id] = treatment

# %%save and export mapping
mapping_file = f'{processed_data_dir}/mouse_treatment_mapping.txt'
with open(mapping_file, 'w') as f:
    for mouse_id, treatment in mouse_mapping.items():
        f.write(f"{mouse_id}: {treatment}\n")
# %%check: how many mice per condition in the mapping
mice_per_condition = {}
for mouse_id, treatment in mouse_mapping.items():
    if treatment not in mice_per_condition:
        mice_per_condition[treatment] = []
    mice_per_condition[treatment].append(mouse_id)  # OK

# %%now per mouse_id, get the images by searching for the mouse_id in the image name
# create a dictionary with the images per mouse ID
images_per_mouse = {}
for image in images_list:
    for mouse_id in mouse_mapping.keys():
        if mouse_id in image:
            if mouse_id not in images_per_mouse:
                images_per_mouse[mouse_id] = []
            images_per_mouse[mouse_id].append(image)
            break  # stop searching for this image, move to the next one
# %%order the images by inLNr, inLNl, brLNr in images_per_mouse.
ln_order = {'inLNr': 0, 'inLNl': 1, 'brLNr': 2}
ln_order_ls = ['inLNr', 'inLNl', 'brLNr']
# Sort the image filenames per mouse according to the defined priority
for mouse_id, images in images_per_mouse.items():
    images_per_mouse[mouse_id] = sorted(
        images,
        key=lambda img: next((ln_order[key] for key in ln_order if key in img), float('inf'))
    )

#%% generate a plot per mouse_id, with the images side by side.
dpi = 50
pixels_per_mm = 3077  # 325 nm/px resolution
scale = 0.5

for mouse_id, images in images_per_mouse.items():
    print(mouse_id)
    # if mouse_id != "HOE-2371":
    #     continue

    treatment = mouse_mapping[mouse_id]
    loaded_images = [mpimg.imread(os.path.join(images_dir, image)) for image in images]
    img_shapes = [img.shape for img in loaded_images]

    max_height = max(shape[0] for shape in img_shapes)
    max_width = max(shape[1] for shape in img_shapes)

    # Pad images to the same size (black background)
    padded_images = []
    for img in loaded_images:
        h, w = img.shape[:2]
        pad_top = (max_height - h) // 2
        pad_bottom = max_height - h - pad_top
        pad_left = (max_width - w) // 2
        pad_right = max_width - w - pad_left

        # Pad with black (0) for RGB images
        padded_img = np.pad(
            img,
            pad_width=((pad_top, pad_bottom), (pad_left, pad_right), (0, 0)),
            mode='constant',
            constant_values=0
        )
        padded_images.append(padded_img)

    total_width_px = len(padded_images) * max_width
    figsize = (scale * total_width_px / dpi, scale * max_height / dpi)  # scale for smaller images

    fig = plt.figure(figsize=figsize, dpi=dpi)
    gs = gridspec.GridSpec(1, len(images), width_ratios=[1] * len(images))

    fig.suptitle(f'Mouse ID: {mouse_id}, Treatment: {treatment}')

    for i, (img, image_name) in enumerate(zip(padded_images, images)):
        ax = fig.add_subplot(gs[0, i])
        ax.imshow(img)
        ax.set_title(f'{ln_order_ls[i]}_{image_name}', fontsize=80)
        ax.axis('off')

        # Add scale bar to each image
        bar_length = int(pixels_per_mm)
        bar_height = int(0.01 * max_height)
        margin = int(0.03 * max_height)

        rect = Rectangle(
            (max_width - margin - bar_length, max_height - margin - bar_height),
            bar_length, bar_height,
            linewidth=0, edgecolor=None, facecolor='white'
        )
        ax.add_patch(rect)

        ax.text(
            max_width - margin - bar_length / 2,
            max_height - margin - bar_height - 10,
            '1 mm',
            color='white', ha='center', va='bottom',
            fontsize=50, weight='bold'
        )

    plt.savefig(os.path.join(fig_dir, f'{mouse_id}_{treatment}_inclbox_scale_{scale}.png'), bbox_inches='tight', dpi=dpi)
    plt.close()
