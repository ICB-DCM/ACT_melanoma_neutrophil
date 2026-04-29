#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
# '''=================================================
# @Project: ACT_melanoma_neutrophil (data: act_neutrophil_data_repository)
# @Author : Gemma van der Voort
# @Desc   : helper script that outputs only legends to add to several figures.
# @Desc updated: Generate standalone legend figures for flow/codex plots.
# ================================================='''


from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
from paths_parameters import bar_colours_hex_codex, marker_per_organ_dict, data_repo_path

plt.style.use(str(Path(__file__).resolve().parent / "journal_style.mplstyle"))
fig_dir = f'{data_repo_path}/figures/flow_cytometry/preprocessing'
os.makedirs(fig_dir, exist_ok=True)

# %%use marker-per-organ-dict to plot legend in three colors: grey, and the final 2 colors in bar_colours_hex_codex. omit
# the 'spleen' entry from the marker-per-organ-dict always. Generate three separate images. Save to pdf
markers_to_plot = {k: v for k, v in marker_per_organ_dict.items() if k != 'Spleen'}
# make a dict with informative names for the colors
colors_to_plot = {'grey': '#808080ff',
                  'yellow': bar_colours_hex_codex[-2],
                  'orange': bar_colours_hex_codex[-1]}
for color_name, color_hex in colors_to_plot.items():
    fig, ax = plt.subplots(figsize=(2, 2))
    for organ, marker in markers_to_plot.items():
        ax.plot([], [], marker=marker, color=color_hex, markeredgecolor='black', markeredgewidth=0.5,
                linestyle='None', label=organ)
    ax.legend(frameon=False, loc='center')
    ax.axis('off')
    plt.savefig(f'{fig_dir}/legend_all_markers_{color_name}.pdf', dpi=300, bbox_inches='tight')
    plt.close()
#%% additionally, plot only the tdLN and clLN entries in the three colors. Save to individual images in dpi=300 .png format.
markers_to_plot_subset = {k: v for k, v in marker_per_organ_dict.items() if k in ['tdLN', 'clLN']}
for color_name, color_hex in colors_to_plot.items():
    fig, ax = plt.subplots(figsize=(2, 2))
    for organ, marker in markers_to_plot_subset.items():
        ax.plot([], [], marker=marker, color=color_hex, markeredgecolor='black', markeredgewidth=0.5,
                linestyle='None', label=organ)
    ax.legend(frameon=False, loc='center')
    ax.axis('off')
    plt.savefig(f'{fig_dir}/legend_tdLN_clLN_{color_name}.pdf', dpi=300, bbox_inches='tight')
    plt.close()
# %% plot all LNs in all three colors. Save to individual images in dpi=300 .png format.
markers_to_plot_LNs = {k: v for k, v in marker_per_organ_dict.items() if k in ['tdLN', 'clLN', 'intLN']}
for color_name, color_hex in colors_to_plot.items():
    fig, ax = plt.subplots(figsize=(2, 2))
    for organ, marker in markers_to_plot_LNs.items():
        ax.plot([], [], marker=marker, color=color_hex, markeredgecolor='black', markeredgewidth=0.5,
                linestyle='None', label=organ)
    ax.legend(frameon=False, loc='center')
    ax.axis('off')
    plt.savefig(f'{fig_dir}/legend_all_LNs_{color_name}.pdf', dpi=300, bbox_inches='tight')
    plt.close()
# %% plot all LNs and tumor
markers_to_plot_LNs = {k: v for k, v in marker_per_organ_dict.items() if k in ['tdLN', 'clLN', 'intLN', 'Tumor']}
for color_name, color_hex in colors_to_plot.items():
    fig, ax = plt.subplots(figsize=(2, 2))
    for organ, marker in markers_to_plot_LNs.items():
        ax.plot([], [], marker=marker, color=color_hex, markeredgecolor='black', markeredgewidth=0.5,
                linestyle='None', label=organ)
    ax.legend(frameon=False, loc='center')
    ax.axis('off')
    plt.savefig(f'{fig_dir}/legend_LN_tumor_{color_name}.pdf', dpi=300, bbox_inches='tight')
    plt.close()
#%% add an additional grey legend containing all LNs, Tumor, and the ratio tdLN/clLN with 'd' the marker
markers_to_plot_incl_ratio = {k: v for k, v in marker_per_organ_dict.items() if k in ['tdLN', 'clLN', 'intLN', 'Tumor']}
markers_to_plot_incl_ratio['tdLN/clLN'] = 'd'
color_hex = '#808080ff'
fig, ax = plt.subplots(figsize=(2, 2))
for organ, marker in markers_to_plot_incl_ratio.items():
    ax.plot([], [], marker=marker, color=color_hex, markeredgecolor='black', markeredgewidth=0.5,
            linestyle='None', label=organ)
ax.legend(frameon=False, loc='center')
ax.axis('off')
plt.savefig(f'{fig_dir}/legend_LN_tumor_ratio_grey.pdf', dpi=300, bbox_inches='tight')
plt.close()
