#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
# '''=================================================
# @Project: ACT_melanoma_neutrophil (data: act_neutrophil_data_repository)
# @Author : Gemma van der Voort
# @Desc   : 20250526: updated without brLNr, with grouped pmel counts.
# @Desc updated: Relative abundance paired analysis for CODEx Pmel CD8 T cells.
# update 20260108: make compositional by dividing by total cd45+ cells in that LN, then do a ratio paired T-test.
# '''=================================================
import os
import sys
from pathlib import Path

import anndata as ad
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

PROJECT_ROOT = Path(__file__).resolve()
while PROJECT_ROOT.name != "ACT_melanoma_neutrophil" and PROJECT_ROOT.parent != PROJECT_ROOT:
    PROJECT_ROOT = PROJECT_ROOT.parent
if PROJECT_ROOT.name != "ACT_melanoma_neutrophil":
    raise RuntimeError("ACT_melanoma_neutrophil root not found")
sys.path.append(str(PROJECT_ROOT / "helper_files"))

from functions_stat_test import paired_test, ratio_paired_test
from functions_vis import boxplot_paired, single_paired_plot, boxplot_paired_no_brLNr
from paths_parameters import x_axis_tick_labels_codex, bar_colours_hex_codex, marker_per_organ_dict, data_repo_path

plt.style.use(str(PROJECT_ROOT / "helper_files" / "journal_style.mplstyle"))

processed_codex_data_dir = f'{data_repo_path}/processed/codex'
processed_spatial_data_dir = f'{data_repo_path}/processed/codex'
stats_dir = f'{processed_spatial_data_dir}/statistics'
os.makedirs(stats_dir, exist_ok=True)
modifier_base = f'20241021_all_ln_tree6_rep2'
modifier_finetuning = 'layer_uns_rem_finetuned_DCs_fixed'
modifier_postprocessing = f'transcleaned_postprocessed'
file_name = f'{modifier_base}_{modifier_finetuning}_{modifier_postprocessing}'
df_grouped_all_cells = pd.read_csv(f'{processed_codex_data_dir}/df_grouped_coi_{file_name}_all_cells.csv', index_col=0)
fig_dir = f'{data_repo_path}/figures/codex/manuscipt_plots/abundances'
os.makedirs(fig_dir, exist_ok=True)
# drop 'Population' column and then remove duplicate rows
df_grouped_all_cells = df_grouped_all_cells.drop(columns='Population').drop_duplicates()
group_col = 'Organ'
condition_col = 'Treatment'
id_col = 'mouse_id'
df_condition_names = df_grouped_all_cells[condition_col].unique()
#%% update 20260108: generate a column with total CD45+ cells per mouse_id and Organ, omitting the non-immune cells.
non_immune_cell_types = ['FRC', 'FDC', 'Endothelial']
df_cd45_counts = (
    df_grouped_all_cells[~df_grouped_all_cells['Cell type'].isin(non_immune_cell_types)]
    .groupby([id_col, group_col, condition_col], as_index=False)
    .agg({'count': 'sum'})
    .rename(columns={'count': 'CD45+ cell count'})
)
# %%merge the total cd45+ counts and back to the original df
df_grouped_all_cells = df_grouped_all_cells.merge(df_cd45_counts, on=[id_col, group_col, condition_col], how='left')
# create a new column with the relative abundance: cell type count / CD45+ cell count
df_grouped_all_cells['relative_abundance_cd45'] = (df_grouped_all_cells['count'] / df_grouped_all_cells['CD45+ cell count']) * 100
#%% end of update 20260108
#%% add 'all Pmel CD8 T' to the df_grouped_all_cells, counting all Pmel CD8 T cells, regardless of cycling or not
df_grouped_all_cells['Higher level cell type'] = None
df_grouped_all_cells.loc[df_grouped_all_cells['Cell type'].str.contains('trans'), 'Higher level cell type'] = 'Pmel CD8 T'
# fill all others as in the original column (to prevent pivot issue)
df_grouped_all_cells['Higher level cell type'].fillna(df_grouped_all_cells['Cell type'], inplace=True)
#%% start update: also count the fraction of Pmel CD8 T cells over total CD45+ cells in that LN
# first, aggregate the counts for Pmel CD8 T cells
df_pmel_agg = (
    df_grouped_all_cells[df_grouped_all_cells['Higher level cell type'] == 'Pmel CD8 T']
    .groupby([id_col, group_col, condition_col], as_index=False)
    .agg({'count': 'sum', 'CD45+ cell count': 'first'})  # take first CD45+ cell count, should be the same for all rows
)
# calculate the fraction
df_pmel_agg['percent_pmel_cd45'] = (df_pmel_agg['count'] / df_pmel_agg['CD45+ cell count']) * 100
# merge back to the main df
df_grouped_all_cells = df_grouped_all_cells.merge(df_pmel_agg[[id_col, group_col, condition_col, 'percent_pmel_cd45']],
                                                  on=[id_col, group_col, condition_col], how='left')

# %% version 1
# phenotyping_column = 'Cell type'
# value_cols = df_grouped_all_cells[phenotyping_column].unique()
# value_cols_names = ['APC Neutrophils', 'B-cells', 'CD169$^+$ Macrophages', 'CD4$^+$ T cells', 'FDC', 'FRC', 'MSM',
#        'Mature B-cells', 'Monocytes', 'NKs', 'Neutrophils', 'RPM', 'Tregs',
#        'Unclassified', 'cDC1', 'cDC2', 'Endogenous CD8$^+$ T other', 'Endogenous CD8$^+$ T cycling',
#        'Endogenous CD8$^+$ T exhausted', 'Mature NK', 'Pmel-1 T other',
#        'Pmel-1 T cycling']
# abundance_column_name = 'relative_abundance_cd45'
# version 2
phenotyping_column = 'Higher level cell type'
value_cols = ['Pmel CD8 T']
value_cols_names = ['Pmel-1 T cells']
abundance_column_name = 'percent_pmel_cd45'
# general
group_order = ['inLNr', 'inLNl']  # 'brLNr'  # order of the lymph nodes in the plot
ln_to_compare = [('inLNr', 'inLNl')]  # ('inLNr', 'brLNr'),  ('brLNr', 'inLNl')
# %%reshape the df to have every unique entry in the "CN" column as a new column
# add a column to the df for each value in the CN column
if phenotyping_column == 'Higher level cell type':
    df_grouped_agg = (
        df_grouped_all_cells[df_grouped_all_cells[phenotyping_column].isin(value_cols)]
        .groupby([id_col, group_col, condition_col, phenotyping_column], as_index=False)
        .agg({abundance_column_name: 'sum'})
    )
    df_grouped_pivotted = df_grouped_agg.pivot(index=[id_col, group_col, condition_col],
                                               columns=phenotyping_column,
                                               values=abundance_column_name).reset_index()
else:
    df_grouped_pivotted = df_grouped_all_cells.pivot(index=[id_col, group_col, condition_col], columns=phenotyping_column,
                                               values=abundance_column_name).reset_index()
df_grouped_coi = df_grouped_pivotted.copy()
# %%stats_df
stats_df = pd.DataFrame(columns=[phenotyping_column, 'Condition', 'LN 1', 'LN 2', 'p-value', 'test_used'])
for value_col in value_cols:
    for condition in df_condition_names:
        for LN1, LN2 in ln_to_compare:
            # for some mice not both lymph nodes were processed. drop these mice
            ln_df_LN1 = df_grouped_coi[(df_grouped_coi[group_col] == LN1) & (df_grouped_coi[condition_col] == condition)]
            ln_df_LN2 = df_grouped_coi[(df_grouped_coi[group_col] == LN2) & (df_grouped_coi[condition_col] == condition)]
            # check if the mouse_ids are the same and in the same order
            if not  ln_df_LN1[id_col].tolist() == ln_df_LN2[id_col].tolist():
                # drop the mouse ID that is not in both conditions
                only_ln1 = set(ln_df_LN1[id_col].tolist()) - set(ln_df_LN2[id_col].tolist())
                only_ln2 = set(ln_df_LN2[id_col].tolist()) - set(ln_df_LN1[id_col].tolist())
                print(f'Only in {LN1}: {only_ln1}, only in {LN2}: {only_ln2}. Dropped single mice')
                ln_df_LN1 = ln_df_LN1[~ln_df_LN1[id_col].isin(only_ln1)]
                ln_df_LN2 = ln_df_LN2[~ln_df_LN2[id_col].isin(only_ln2)]
            else:
                # print(f'Mice in {LN1} and {LN2} are the same')
                pass

            t_stat, p_val, test_used = ratio_paired_test(ln_df_LN1, ln_df_LN2, value_col, non_parametric=False)
            stats_df = pd.concat([stats_df, pd.DataFrame([[value_col, condition, LN1, LN2, p_val,
                                                           test_used]],
                                                         columns=[phenotyping_column, 'Condition', 'LN 1', 'LN 2', 'p-value',
                                                                  'test_used'])])
# write to file
stats_df.to_csv(f'{stats_dir}/codex_relative_abundance_paired_stats_no_brLNr.csv')

# %% plot the data
n_categories = len(df_condition_names)
ncols = 1
fig_width = ((n_categories*1.5)-1)*ncols  # a bit broader because we have 3 entries per condition
markers_LN = [marker_per_organ_dict[ln] for ln in group_order]
for cell_type_i, cell_type in enumerate(value_cols):
    print('Working on', cell_type)
    for LN1, LN2 in ln_to_compare:
        try:
            fig, ax = plt.subplots(1, ncols)  # , figsize=(fig_width, 6))
            data = df_grouped_coi
            # get stats df for cell type and LN combination
            stats_df_cell = stats_df[(stats_df[phenotyping_column] == cell_type) & (stats_df['LN 1'] == LN1) & (stats_df['LN 2']
                                                                                                     == LN2)]
            stats_df_cell = stats_df_cell.drop([phenotyping_column, 'Condition'], axis=1)  # drop for compatibility with plot
            boxplot_paired_no_brLNr(data, condition_col, cell_type, x_group_col=group_col, stats_df=stats_df_cell, ax=ax,
                           x_group_order=group_order, markers=markers_LN, condition_colours=bar_colours_hex_codex,
                                    omit_ns=True)
            # add title with the LN combination
            # ax.set_title(f'Lymph nodes, paired test {LN1} and {LN2}', y=1.05)
            ax.set_ylabel(f'{value_cols_names[cell_type_i]}\n[% of CD45+ cells]')
            # ax.set_xlabel('Condition')
            ax.set_xlabel('')
            ax.set_xticklabels(x_axis_tick_labels_codex)
            plt.tight_layout()
            plt.savefig(f'{fig_dir}/paired_{cell_type}_{LN1}_{LN2}_no_brLNr_omit_ns_narrow_relative_abundance.pdf',
                        bbox_inches='tight', dpi=300)
            plt.close()
        except TypeError:
            print('No data for', cell_type, 'in', LN1, 'and', LN2)
            continue
