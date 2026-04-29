#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
# '''=================================================
# @Project: ACT_melanoma_neutrophil (data: act_neutrophil_data_repository)
# @Author : Gemma van der Voort
# @Desc   : uses ome.tiff and marker csv to create an overview image of the channels of interest
# @Desc updated: Generate CODEx overview intensity images for region detection.
# '''=================================================

# %%imports
import os
import sys
from pathlib import Path
import numpy as np
import tifffile as tiff
from tqdm import tqdm
import pathlib
import matplotlib.pyplot as plt
from skimage import exposure, util

PROJECT_ROOT = Path(__file__).resolve()
while PROJECT_ROOT.name != "ACT_melanoma_neutrophil" and PROJECT_ROOT.parent != PROJECT_ROOT:
    PROJECT_ROOT = PROJECT_ROOT.parent
if PROJECT_ROOT.name != "ACT_melanoma_neutrophil":
    raise RuntimeError("ACT_melanoma_neutrophil root not found")
sys.path.append(str(PROJECT_ROOT / "helper_files"))

from utils_codex import prune_csv_to_markers
from paths_parameters import data_repo_path

# locations
processed_image_folder = f'{data_repo_path}/raw/codex'
figure_dir = f'{data_repo_path}/figures/codex/overview_intensity_images'

# marker_list = ["DAPI", "Tox", "Tcf1", "NKp46", "FoxP3", "Eomes", "CD31", "Bcl2", "Ki67", "CD21-35", "Tbet",
#                "Granzyme B", "CD11b", "CD137", "Ly108", "CD4", "Tim3", "Lag3", "B220", "CD90.1", "CD11c", "MHCII",
#                "CD8", "CD45", "F4-80", "CD103", "CD169", "aSMA", "PD1", "CD69", "CD3", "Ly6C", "ERTR7", "Vimentin",
#                "Ly6G"]

# %%MODIFY THIS PART TO YOUR NEEDS
# img_folder_name ='b_t_cell_zones'
# markers_of_interest = ["DAPI", "B220", "CD3"]
# colors = {
#     'DAPI': (192/255, 192/255, 192/255),        # Grey
#     'B220': (55/255, 126/255, 184/255),         # Blue
#     'CD3': (255/255, 255/255, 0)              # Yellow
# }
# img_folder_name ='membrane_markers'
# markers_of_interest = ["DAPI", "CD45", "aSMA", "CD21-35", "CD31"]
# colors = {
#     'DAPI': (192/255, 192/255, 192/255),        # Grey
#     'CD45': (55/255, 126/255, 184/255),         # Blue
#     'aSMA': (152/255, 78/255, 163/255),         # Purple
#     'CD21-35': (77/255, 175/255, 74/255),       # Green
#     'CD31': (255/255, 127/255, 0)               # Orange
# }
img_folder_name ='full_ln_architecture'
markers_of_interest = ["DAPI", "B220", "aSMA", "CD3", "ERTR7"]  # grey, red, blue, green and yellow, respectively
colors = {
    'DAPI': (192/255, 192/255, 192/255),        # Grey, all nuclei
    'B220': (255/255, 0, 0),                    # Red, B cells
    'aSMA': (0, 0, 255/255),                    # Blue, smooth muscle actin
    'CD3': (0, 255/255, 0),                     # Green, T cells
    'ERTR7': (255/255, 255/255, 0)               # Yellow, reticular fibroblastic network (surrounds vessels)
}
dapi_images = False  # whether to make additional images of only the dapi channel
# %%standard from here again:
# make the figure dir if it does not exist
full_fig_dir = f'{figure_dir}/{img_folder_name}'
pathlib.Path(figure_dir).mkdir(parents=True, exist_ok=True)
# also make subdirectories for figures
if dapi_images:
    subdirs = ['dapi_channel_lowQ', 'dapi_channel_highQ', 'inverted_dapi_lowQ', 'channels_of_interest_lowQ',
               'channels_of_interest_highQ']
else:
    subdirs = ['channels_of_interest_lowQ', 'channels_of_interest_highQ']
for subdir in subdirs:
    pathlib.Path(f'{full_fig_dir}/{subdir}').mkdir(parents=True, exist_ok=True)
# list of files
folder_list = os.listdir(processed_image_folder)
# %% for each file (ome.tiff)
for folder in tqdm(folder_list):
    # open .ome.tiff in folder. name can vary, but only one with this extension
    marker_csv = f'{processed_image_folder}/{folder}/csv_unfiltered_new_filter/R1_marker_headers.csv'
    # read the first line of csv
    with open(marker_csv, 'r') as file:
        first_line = file.readline()
        indices, marker_list_csv = prune_csv_to_markers(first_line, markers_of_interest)
    print(indices)

    file_list = os.listdir(f'{processed_image_folder}/{folder}')
    for file in file_list:
        if file.endswith('.ome.tif'):
            image = tiff.imread(f'{processed_image_folder}/{folder}/{file}')
            print(folder)
            print(image.shape)
            assert image.shape[0] == len(marker_list_csv), (f'Image shape 0 (number of channels) does not match marker '
                                                            f'list in csv file. \n Number of channels in image: '
                                                            f'{image.shape[0]}. \n Number of markers in csv file: '
                                                            f'{len(marker_list_csv)}')
            # extract the channels of interest by index
            image_channels = image[indices]
            # reorder dimensions to height, width, channels
            image_channels = np.moveaxis(image_channels, 0, -1)
            # Create an empty RGB image
            rgb_image = np.zeros((image.shape[1], image.shape[2], 3), dtype=np.float32)
            # Map each channel to its corresponding color
            for i in range(len(markers_of_interest)):
                channel = image_channels[:, :, i]
                # rescale the intensity of the channel by removing peaks and normalizing to 0-1
                vmin, vmax = np.percentile(channel, q=(0.5, 99.5))
                channel_normalized = exposure.rescale_intensity(
                    channel, in_range=(vmin, vmax), out_range=np.float32
                )
                channel_normalized = ((channel_normalized - channel_normalized.min()) /
                                      (channel_normalized.max() - channel_normalized.min()))
                # check: plot a histogram of the channel
                # plt.hist(channel_normalized.flatten(), bins=100)
                # plt.title(f'{folder} {markers_of_interest[i]}')
                # plt.show()
                # plt.close()
                # if channel is dapi, make an image of just the dapi channel
                if dapi_images and markers_of_interest[i] == 'DAPI':
                    plt.figure(figsize=(10, 10))
                    plt.axis('off')
                    plt.imshow(channel_normalized, cmap='gray')
                    plt.imsave(f'{full_fig_dir}/dapi_channel_lowQ/{folder}_dapi_channel_lowQ.jpg', channel_normalized, cmap='gray',
                               dpi=150)
                    plt.imsave(f'{full_fig_dir}/dapi_channel_highQ/{folder}_dapi_channel.png', channel_normalized, cmap='gray',
                               dpi=300)
                    plt.close()
                    # flip the pixels so the background is white and the channel is black
                    channel_normalized_inv = 1 - channel_normalized
                    plt.figure(figsize=(10, 10))
                    plt.axis('off')
                    plt.imsave(f'{full_fig_dir}/inverted_dapi_lowQ/{folder}_dapi_channel_inv_lowQ.jpg',
                               channel_normalized_inv,
                               cmap='gray',
                               dpi=150)
                    plt.close()

                # for the multichannel image, assign colors based on the channel index
                rgb_image[..., 0] += channel_normalized * colors[list(colors.keys())[i]][0]  # Red channel
                rgb_image[..., 1] += channel_normalized * colors[list(colors.keys())[i]][1]  # Green channel
                rgb_image[..., 2] += channel_normalized * colors[list(colors.keys())[i]][2]  # Blue channel

            # Clip values to ensure they are valid for an image
            rgb_image = np.clip(rgb_image, 0, 1)

            # Display the image using imshow
            plt.figure(figsize=(10, 10))
            plt.axis('off')  # Turn off axis labels
            # markers_of_interest as a string
            markers_of_interest_str = '_'.join(markers_of_interest)
            plt.imsave(f'{full_fig_dir}/channels_of_interest_highQ/{folder}_{markers_of_interest_str}.png', rgb_image, dpi=300)
            plt.imsave(f'{full_fig_dir}/channels_of_interest_lowQ/{folder}_{markers_of_interest_str}_lowQ.jpg', rgb_image, dpi=150)
            plt.close()

    # later: edge detection and calculate features of artefacts
#%% dapi in print-friendly format
# 71 images. make 12 images per page. 6 pages. 3 columns, 4 rows. A4 format (for printing): dapi inv
if dapi_images:
    folder_list.sort()  # order folder list alphabetically
    nrows = 4
    ncols = 3
    # infer amount of pages
    pages = len(folder_list) // (nrows*ncols) + 1
    for page_index in tqdm(range(pages)):
        fig, axs = plt.subplots(nrows, ncols, figsize=(20, 20))
        for ax_index, ax in enumerate(axs.flatten()):
            if page_index * nrows * ncols + ax_index < len(folder_list):
                folder = folder_list[page_index * nrows * ncols + ax_index]
                print(ax_index, folder)
                image = plt.imread(f'{full_fig_dir}/inverted_dapi_lowQ/{folder}_dapi_channel_inv_lowQ.jpg')
                ax.imshow(image, cmap='gray')
                ax.axis('off')
                ax.title.set_text(folder)
        plt.tight_layout()
        plt.savefig(f'{full_fig_dir}/dapi_inv_lowQ_page_{page_index}.png', dpi=300)
        plt.close()
