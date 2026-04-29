#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
# '''=================================================
# @Project: ACT_melanoma_neutrophil (data: act_neutrophil_data_repository)
# @Author : Gemma van der Voort
# @Desc   : Patch proximity analysis setup for CODEx data. revision: rerun with smoothing max_iteration 25
# note: switch to spacec venv for this script
# '''=================================================
# imports
import sys
from pathlib import Path
import anndata as ad
import spacec as sp
import warnings
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve()
while PROJECT_ROOT.name != "ACT_melanoma_neutrophil" and PROJECT_ROOT.parent != PROJECT_ROOT:
    PROJECT_ROOT = PROJECT_ROOT.parent
if PROJECT_ROOT.name != "ACT_melanoma_neutrophil":
    raise RuntimeError("ACT_melanoma_neutrophil root not found")
sys.path.append(str(PROJECT_ROOT / "helper_files"))

from paths_parameters import fig_dir, modifier_base, parent_dir, processed_spatial_data_dir, cell_colours_hex
# %%custom paths and parameters
cell_type_col = 'Cell type'
cell_type_col_file_name = 'cell_type'  # for use in file names, no spaces or caps. also poss shorter
version_key = 'v4_8'  # version of the metacluster key, large scale classifier.
cn_col = 'CN_k50_n20'  # column name for the CN in adata

metacluster_key = 'Metacluster v4_8_filtered_0.5_maxit_25_v4_8'
img_modifier = 'v4_corrected'
adata = ad.read_h5ad(f"{processed_spatial_data_dir}/sa_{modifier_base}_inclregions_{img_modifier}_"
                     f"{metacluster_key}.h5ad")  # includes regions
#%%
results, outlines_results = sp.tl.patch_proximity_analysis(
    adata,
    region_column = "dataset_name",
    patch_column = metacluster_key,
    group='T-cell zone',
    min_cluster_size=30,
    x_column='X', y_column='Y',
    radius = [20, 60, 100], # radii to calculate
    edge_neighbours = 1,
    plot = False,
    original_unit_scale = round(1/0.325), #  325 nm/px
    method= "border_cell_radius",
    key_name = "ppa_result_20_60_100_border_cell_radius",
    save_geojson = False,
    )

#%% save adata
adata.write(f"{processed_spatial_data_dir}/sa_{modifier_base}_inclregions_{img_modifier}_"
            f"{metacluster_key}_ppa.h5ad")
# #%% donut plot (unused)
# for treatment in adata.obs['Treatment'].unique():
#     print(treatment)
#     adata_treatment = adata[adata.obs['Treatment'] == treatment].copy()
#     sp.pl.ppa_res_donut(
#         adata_treatment,
#         cat_col = 'Cell type',
#         key_name="ppa_result_20_60_100_border_cell_radius",
#         palette=None,
#         distance_mode="between",  # "within" or "between"
#         unit="µm",
#         figsize=(10, 10),
#         add_guides=True,
#         text="Cell types around T-cell zones, within",
#         label_color="black",
#         group_by= 'Organ',
#         title="PPA",
#     )


#%% for all treatments, count the abundance of neutrophils within 20, 60, 100 µm of T cell zones
# and normalize to the total number of neutrophils in that sample. make a df with columns  ['Treatment', 'Organ',
# 'mouse_id', abundance_below_20, abundance_above_20]

df_t_zone_proximity = pd.DataFrame() # columns=columns)
for treatment in adata.obs['Treatment'].unique():
    print(treatment)
    adata_treatment = adata[adata.obs['Treatment'] == treatment].copy()
    results_treatment = results[results['Treatment'] == treatment].copy()
    for organ in adata_treatment.obs['Organ'].unique():
        print(organ)
        adata_treatment_organ = adata_treatment[adata_treatment.obs['Organ'] == organ].copy()
        results_treatment_organ = results_treatment[results_treatment['Organ'] == organ].copy()
        for mouse_id in adata_treatment_organ.obs['mouse_id'].unique():
            print(mouse_id)
            adata_treatment_organ_mouse = adata_treatment_organ[adata_treatment_organ.obs['mouse_id'] == mouse_id].copy()
            results_treatment_organ_mouse = results_treatment_organ[results_treatment_organ['mouse_id'] == mouse_id].copy()
            total_neutrophils = len(adata_treatment_organ_mouse[adata_treatment_organ_mouse.obs[cell_type_col] == 'Neutrophil'])
            # calculate the number of neutrophils within the T cell zone (so we can correct for it)
            neut_in_t_cell_zone = len(adata_treatment_organ_mouse[
                                          (adata_treatment_organ_mouse.obs[cell_type_col] == 'Neutrophil') &
                                          (adata_treatment_organ_mouse.obs[metacluster_key] == 'T-cell zone')])
            total_external_neutrophils = total_neutrophils - neut_in_t_cell_zone
            results_treatment_organ_mouse_neutrophils = results_treatment_organ_mouse[results_treatment_organ_mouse[cell_type_col] == 'Neutrophil'].copy()
            total_neut_in_radii = len(results_treatment_organ_mouse_neutrophils)
            print('Total neutrophils:', total_neutrophils)
            print('Neutrophils in all radii:', total_neut_in_radii)
            print('Neutrophils in T cell zone:', neut_in_t_cell_zone)
            print('Neutrophils outside T cell zone:', total_external_neutrophils)
            if total_neut_in_radii> total_external_neutrophils:
                # this can happen if a radius overlaps with the next t cell zone, counting cells double.
                # we take the total external neutrophils as the maximum to correct for this.
                total_neut_in_radii = total_external_neutrophils
            if total_neutrophils == 0:
                warnings.warn(f"No neutrophils found in {treatment}, {organ}, {mouse_id}. Skipping.")
                continue
            abundances = {}
            for distance in [0, 60, 180, 300]:  # !!! the distances saved in the results dataframe are in px!
                if distance == 0:
                    # abundance within T cell zone
                    abundances[f'abundance_below_{distance}_px'] = neut_in_t_cell_zone / total_neutrophils
                    print(f"Neutrophils within T cell zone:", neut_in_t_cell_zone)
                    abundances[f'abundance_above_{distance}_px'] = total_external_neutrophils / total_neutrophils
                    print(f"Neutrophils outside T cell zone:", total_external_neutrophils)
                    continue
                else:
                    abundances[f'abundance_below_{distance}_px'] = len(results_treatment_organ_mouse_neutrophils[
                        results_treatment_organ_mouse_neutrophils[f'distance_from_patch'] <= distance]) / total_external_neutrophils
                    if abundances[f'abundance_below_{distance}_px'] > 1:
                        abundances[f'abundance_below_{distance}_px'] = 1  # can happen due to overlapping radii
                    abundances[f'abundance_above_{distance}_px'] = 1-abundances[f'abundance_below_{distance}_px']
            print(f"Abundances for {treatment}, {organ}, {mouse_id}:")
            print(abundances)
            row = {
                'Treatment': treatment,
                'Organ': organ,
                'mouse_id': mouse_id,
            }
            row.update(abundances)
            df_t_zone_proximity = pd.concat([df_t_zone_proximity, pd.DataFrame([row])], ignore_index=True)
#%% save df
df_t_zone_proximity.to_csv(f"{processed_spatial_data_dir}/df_t_zone_proximity_neutrophils_{img_modifier}_{metacluster_key}.csv")
