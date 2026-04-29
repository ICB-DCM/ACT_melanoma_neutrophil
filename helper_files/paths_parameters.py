#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
# '''=================================================
# @Project: ACT_melanoma_neutrophil (data: act_neutrophil_data_repository)
# @Author : Gemma van der Voort
# @Desc   : paths and parameters for this project, to be imported in other scripts.
# @Desc updated: Centralized data-repository paths and shared plotting parameters.
# '''=================================================
# common path
data_repo_path = "/media/gemma/X10 Pro/act_neutrophil_data_repository"
#%% codex paths etc
parent_dir = f"{data_repo_path}/raw/codex"
processed_data_dir = f'{data_repo_path}/processed/codex'
processed_spatial_data_dir = f'{data_repo_path}/processed/codex'
fig_dir = f'{data_repo_path}/figures/codex'
# fig_dir = f'{data_repo_path}/figures/codex/spatial_analysis' # should not need after renaming
modifier_base = f'20241021_all_ln_tree6_rep2'
modifier = f'{modifier_base}_layer_uns_rem_finetuned_DCs_fixed_transcleaned_postprocessed'
modifier_finetuning = 'layer_uns_rem_finetuned_DCs_fixed'
modifier_postprocessing = 'transcleaned_postprocessed'
data_fraction = 1

# %%plotting conventions, all modalities
x_axis_tick_labels_blood = ["Na\u00EFve",
                            "Tumor\nbearing", "Cy", "ACT\nday 3 \n(\u00B11)",
                      "ACT\nday 7\n(\u00B12)", "ACT\nday 14\n(\u00B12)", "Relapse"]
bar_colours_hex_blood_flow = ["#8f3ca6ff",  # lighter purple, naive
                         "#5a5eb8ff",  # lighter dark blue, tumour bearing
                         "#4c8ac1ff",  # lighter blue, cyclo
                         "#1f968bff",  # turquoise, ACT day 3
                         "#fde725ff",  # yellow, ACT day 7
                         "#fa815fff",  # orange, ACT day 14
                         "#f71480ff"   # pink, relapse
                        ]
experiments_flow_noCpG = ['69', '70', '76', '86', '96', '97', '87']
conditions_flow_no_CpG = ['Naive', 'Tumour-bearing 3-5 mm', 'Cy only', 'ACT-d3', 'ACT-d7',
       'ACT-d14', 'ACT_Relapse']
x_axis_tick_labels_flow = ["Na\u00EFve", "Tumor\nbearing", "Cy", "ACT\nday 3",
                      "ACT\nday 7", "ACT\nday 14", "Relapse"]  # us english!
titles_tumour_growth = ["Tumor bearing", "Cy", "ACT day 3",
                      "ACT day 7", "ACT day 14", "Relapse", "ACT day 7, no CpG", "ACT day 14, no CpG"]
colours_tumour_growth = bar_colours_hex_blood_flow[1:]
x_axis_tick_labels_codex = x_axis_tick_labels_flow[1:-1]
bar_colours_hex_codex = bar_colours_hex_blood_flow[1:-1]
lymph_nodes_label = ['inLNr\n(tumor-draining)', 'inLNl\n(contralateral)', 'brLNr\n(intermediate)']
organ_names_short = ['inLNr', 'inLNl', 'brLNr', 'Spleen', 'Blood', 'Tumor']
organ_names_long = ['inLNr\n(tumor-draining)', 'inLNl\n(contralateral)', 'brLNr\n(intermediate)', 'spleen', 'blood', 'tumor']
marker_per_organ = ['o', 's', '^', 'D', 'X', 'P']
marker_per_organ_dict = dict(zip(organ_names_short, marker_per_organ))
#%% ln names updated
ln_names_data = ['inLNr', 'inLNl', 'brLNr']
ln_labels_updated = ['tdLN', 'clLN', 'intLN']
updated_ln_label_dict = dict(zip(ln_names_data, ln_labels_updated))

# naming conventions flow cytometry
neutrophil_label = 'Neutrophils\n[% of CD45.2$^{+}$ cells]'
teff_label = 'Pmel-1 $\mathrm{T_{EFF}}$\n[% of Pmel-1 T cells]'
tcm_label = 'Pmel-1 $\mathrm{T_{CM}}$\n[% of Pmel-1 T cells]'
#%% all possible cell types for codex phenotyping (including higher level phenotypes for uncertain clusters)
possible_cell_types_codex = ["Non Immune", "Myeloid", "Lymphoid", "Endothelial", "FRC", "FDC", "Monocyte", "cDC2", "RPM",
                  "CD169+ Macrophage", "MSM", "Neutrophil", "NK", "B-cell", "T-cell", "cDC1", "Mature NK",
                  "Mature B-cell", "CD4", "endo CD8", 'trans CD8', "APC Neutrophil", "Treg", "endo T cycling",
                  "endo T exhausted", 'trans T cycling', 'trans T exhausted']
# %%Dictionary with hex codes for codex cell types (umap, patches, etc
cell_colours_hex = {
    "Mature B-cell": "#8E44AD",  # Rich purple
    "B-cell": "#A569BD",  # Lighter purple for contrast

    "CD169+ M": "#8C564B",
    "Monocyte": "#C49C94",
    "MSM": "#008080",  # Dark teal
    "RPM": "#B8860B",  # Dark orange

    "cDC1": "#2CA02C",
    "cDC2": "#98DF8A",

    "Neutrophil": "#D62728",  # Bright red
    "APC Neutrophil": "#FF9896",  # Lighter red

    "FRC": "#C5B0D5",
    "FDC": "#8C6D31",

    "NK": "#17BECF",  # Cyan
    "mature NK": "#9EDAE5",  # Light cyan

    "CD4": "#F7B6D2",  # Pink
    "Treg": "#E377C2",  # Magenta

    "endo CD8": "#BCBD22",  # Olive green
    "endo T cycling": "#DBDB8D",  # Light olive
    "endo T exhausted": "#556B2F",  # Dark olive

    "trans CD8": "#1F78B4",  # Deep blue
    "trans T cycling": "#A6CEE3",  # Light blue

    "Unclassified": "#BDBDBD"  # Neutral gray
}
