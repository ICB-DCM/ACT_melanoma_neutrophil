#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
# '''=================================================
# @Project: ACT_melanoma_neutrophil (data: act_neutrophil_data_repository)
# @Author : Gemma van der Voort
# @Desc   : to compare the different codex phenotyping outcomes to the already established flow pan-immune panels.
# @Desc updated: Flow vs CODEx abundance comparison for manuscript plots.
# '''=================================================

import sys
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import anndata as ad
import scipy.stats as stats
from tqdm import tqdm
from matplotlib.patches import Ellipse

PROJECT_ROOT = Path(__file__).resolve()
while PROJECT_ROOT.name != "ACT_melanoma_neutrophil" and PROJECT_ROOT.parent != PROJECT_ROOT:
    PROJECT_ROOT = PROJECT_ROOT.parent
if PROJECT_ROOT.name != "ACT_melanoma_neutrophil":
    raise RuntimeError("ACT_melanoma_neutrophil root not found")
sys.path.append(str(PROJECT_ROOT / "helper_files"))

from utils_codex import adata_to_df
from paths_parameters import processed_data_dir, modifier_base, modifier_postprocessing, modifier_finetuning, data_repo_path

plt.style.use(str(PROJECT_ROOT / "helper_files" / "journal_style.mplstyle"))
modifier_base = '20241021_all_ln_tree6_rep2'  # because v7 is not done yet.
# %% load the flow data
flow_data_dir = f'{data_repo_path}/raw/flow_cytometry/pan_immune'
# load dataset, skipping the first 5 rows and setting header to 6th row. 1st column is the highest level index, 2nd
# column is the second level index, 3rd column is the index per row.
flow_data_raw = pd.read_excel(f'{flow_data_dir}/Cell Counts_Pan-Immune_20240219.xlsx', sheet_name='LN_Pan-Immune',
                              skiprows=5, index_col=[0, 1])
treatments = ['tumour_bearing', 'Cyclo_d1', 'ACT_d3', 'ACT_d7', 'ACT_d14']
organs = ['inLNl', 'brLNr', 'inLNr']
# %% load example codex data
codex_modifiers = [f'{modifier_base}_{modifier_finetuning}_{modifier_postprocessing}']
fig_dir = f'{data_repo_path}/figures/flow_codex_comparison'
os.makedirs(fig_dir, exist_ok=True)
# %% clean flow
# %%add numeric index, change the two level index to two extra regular columns
flow_data_2 = flow_data_raw.reset_index()
# rename 1st 3 columns to 'experiment', 'treatment' and 'file name'
flow_data_2 = flow_data_2.rename(columns={'level_0': 'experiment', 'level_1': 'Treatment',
                                          'Unnamed: 2': 'file_name'})
# drop columns that don't contain 'Count normalized' (except first 3 columns)
flow_data_2 = flow_data_2.loc[:, flow_data_2.columns.str.contains('Count normalised|experiment|Treatment|file_name')]
# we need to split column file_name into two columns: 'Organ' and 'mouse_id'
# split by _ and add the first part to 'Organ' and the last part to 'mouse_id' (plus remove extension in mouse
flow_data_2['Organ'] = flow_data_2['file_name'].str.split('_').str[0]
flow_data_2['mouse_id'] = flow_data_2['file_name'].str.split('_').str[-1]
flow_data_2['mouse_id'] = flow_data_2['mouse_id'].str.replace('.fcs', '', regex=False)
for alt_spelling in ['brLNri', 'brLN-r', 'brLN-right', 'brLNright']:
    flow_data_2['Organ'] = flow_data_2['Organ'].str.replace(alt_spelling, 'brLNr')
for alt_spelling in ['inLNle', 'inLN-l', 'inLN-left', 'inLNleft']:
    flow_data_2['Organ'] = flow_data_2['Organ'].str.replace(alt_spelling, 'inLNl')
for alt_spelling in ['inLNri', 'inLN-r', 'inLN-right', 'inLNright']:
    flow_data_2['Organ'] = flow_data_2['Organ'].str.replace(alt_spelling, 'inLNr')
# remove 'Count normalized' from column names
flow_data_2.columns = flow_data_2.columns.str.replace('Count normalised', '')
flow_data_2.columns = flow_data_2.columns.str.replace('|', '', regex=False)
# remove trailing spaces from column names
flow_data_2.columns = flow_data_2.columns.str.strip()
# drop treatments that don't exist and rename the rest to match codex
flow_data_2 = flow_data_2.loc[
    ~flow_data_2['Treatment'].isin(['Naive', 'ACT-d7_no CpG', 'ACT-d14_no CpG', 'ACT_Relapse'])]
flow_data_2['Treatment'] = flow_data_2['Treatment'].replace('Tumour-bearing 3-5 mm', 'tumour_bearing')
flow_data_2['Treatment'] = flow_data_2['Treatment'].replace('Cyclo only', 'Cyclo_d1')
flow_data_2['Treatment'] = flow_data_2['Treatment'].replace('ACT-d3', 'ACT_d3')
flow_data_2['Treatment'] = flow_data_2['Treatment'].replace('ACT-d7', 'ACT_d7')
flow_data_2['Treatment'] = flow_data_2['Treatment'].replace('ACT-d14', 'ACT_d14')
# %%
for modifier in tqdm(codex_modifiers):
    file_name = (f'phenotyping_{modifier}')
    adata = ad.read_h5ad(f"{processed_data_dir}/{file_name}.h5ad")
    # %% create grouped df from codex data in same format as flow data
    try:
        phenotyping_levels = ['name', 'Cell type']  # adata to df renames compartment to name
        codex_df = adata_to_df(adata, additional_columns=['Organ', 'Treatment', 'mouse_id', 'Cell type'])
    except KeyError:  # for non-fine tuned data
        phenotyping_levels = ['Compartment']
        codex_df = adata_to_df(adata, additional_columns=['Organ', 'Treatment', 'mouse_id'])
    codex_df = codex_df.drop(columns=['X', 'Y'])  # drop x and y columns
    # %%
    phenotyping_levels = ['Cell type']
    for phenotyping_level in tqdm(phenotyping_levels):
        # drop the other phenotyping columns
        other_phenotyping_levels = [level for level in phenotyping_levels if level != phenotyping_level]
        codex_df_2 = codex_df.drop(columns=other_phenotyping_levels)
        codex_df_3 = codex_df_2.groupby(['mouse_id', 'Treatment', 'Organ']).value_counts().unstack().reset_index()
        # remove all rows where all values are 0 for all columns except the first 3
        codex_df_3 = codex_df_3.loc[codex_df_3.iloc[:, 4:].sum(axis=1) != 0]

        # %% create a dict tom map the flow and the codex cell types to consensus cell types
        cell_type_map_consensus_flow = {'DC': ['DCs'],
                                        # 'CD11b+ Dcs', 'Xcr1+ Dcs': don't include these since we would be
                                        # counting them twice
                                        'Endo CD8 T': ['Endo CD8 T'],
                                        'Trans CD8 T': ['Pmel T'],
                                        'CD4 T': ['CD4 Tconv', 'CD25+ CD4+ T'],
                                        'Neutrophil': ['Neutrophils'],
                                        'NK': ['Nk cells'],
                                        'B-cell': ['B cells']}
        cell_type_map_consensus_codex = {'DC': ['cDC1', 'cDC2'],
                                         'Endo CD8 T': ['endo CD8', 'endo T cycling', 'endo T exhausted'],
                                         'Trans CD8 T': ['trans CD8', 'trans T cycling', 'trans T exhausted'],
                                         'CD4 T': ['CD4', 'Treg'],
                                         'Neutrophil': ['Neutrophil', 'APC Neutrophil'],
                                         'NK': ['NK', 'Mature NK'],
                                         'B-cell': ['B-cell', 'Mature B-cell']}
        # modify codex and flow columns to match consensus
        consensus_cell_types = list(cell_type_map_consensus_flow.keys())
        flow_data_3 = flow_data_2.loc[:, ['Treatment', 'Organ', 'mouse_id']]
        codex_df_4 = codex_df_3.loc[:, ['Treatment', 'Organ', 'mouse_id']]
        for consensus_cell_type, cell_types in cell_type_map_consensus_flow.items():
            flow_data_3[consensus_cell_type] = flow_data_2[cell_types].sum(axis=1)
        for consensus_cell_type, cell_types in cell_type_map_consensus_codex.items():
            try:
                codex_df_4[consensus_cell_type] = codex_df_3[cell_types].sum(axis=1)
            except KeyError:
                # try with each cell type separately
                for cell_type in cell_types:
                    try:
                        codex_df_4[consensus_cell_type] = codex_df_3[cell_type]
                    except KeyError:
                        print(f'Cell type {cell_type} not found in codex data')
        # we now have a 1383 by 10 df for codex (65 samples) and a 140 by 10 df for flow (140 samples) in the same
        # format.
        # group each by treatment and organ and plot the dataframes as a scatterplot side by side (should have same
        # shape after grouping). Keep the variance resulting from the different mouse samples in the plot.
        # %% group by treatment and organ and mouse
        codex_df_5 = codex_df_4.groupby(['Treatment', 'Organ', 'mouse_id']).sum().reset_index()  # totals per mouse
        flow_data_4 = flow_data_3.groupby(['Treatment', 'Organ', 'mouse_id']).sum().reset_index()
        # %%remove all rows where all values are 0 for all columns except the first 3
        codex_df_5 = codex_df_5.loc[codex_df_5.iloc[:, 3:].sum(axis=1) != 0]
        # fraction of cells of interest
        codex_df_5['total'] = codex_df_5.iloc[:, 3:].sum(axis=1)
        flow_data_4 = flow_data_4.loc[flow_data_4.iloc[:, 3:].sum(axis=1) != 0]
        # fraction of cells of interest
        flow_data_4['total'] = flow_data_4.iloc[:, 3:].sum(axis=1)
        for cell_type in consensus_cell_types:
            codex_df_5[f'{cell_type}_fraction'] = codex_df_5[cell_type] / codex_df_5['total']
            flow_data_4[f'{cell_type}_fraction'] = flow_data_4[cell_type] / flow_data_4['total']

        # group the dataframes by treatment and organ, calculate the mean and standard deviation of the fractions.
        # %%there is still a categorical column mouse_id left, correct for that
        # codex_frac = codex_df_5.groupby(['Treatment', 'Organ'], observed=True).agg(['mean', 'std'], numeric_only=True)
        categorical_cols = ['Treatment', 'Organ', 'mouse_id']
        numerical_cols_codex = codex_df_5.select_dtypes(include=['number']).columns
        agg_dict_codex = {col: ['mean', 'std'] for col in numerical_cols_codex}
        codex_frac = codex_df_5.groupby(['Treatment', 'Organ'], observed=True).agg(agg_dict_codex)
        # %%
        # flow_frac = flow_data_4.groupby(['Treatment', 'Organ'], observed=True).agg(['mean', 'std'], numeric_only=True)
        numerical_cols_flow = flow_data_4.select_dtypes(include=['number']).columns
        agg_dict_flow = {col: ['mean', 'std'] for col in numerical_cols_flow}
        flow_frac = flow_data_4.groupby(['Treatment', 'Organ'], observed=True).agg(agg_dict_flow)

        # %%set each fraction in codex against the corresponding fraction in flow
        plot_df = pd.DataFrame(columns=['Treatment', 'Organ', 'cell_type', 'mean codex', 'mean flow', 'std codex',
                                        'std flow'])
        for treatment in treatments:
            for organ in organs:
                for cell_type in consensus_cell_types:
                    new_row = pd.DataFrame({'Treatment': [treatment], 'Organ': [organ], 'cell_type': [cell_type],
                                            'mean codex': [codex_frac.loc[(treatment, organ), (f'{cell_type}_fraction',
                                                                                               'mean')]],
                                            'mean flow': [flow_frac.loc[(treatment, organ), (f'{cell_type}_fraction',
                                                                                             'mean')]],
                                            'std codex': [codex_frac.loc[(treatment, organ), (f'{cell_type}_fraction',
                                                                                              'std')]],
                                            'std flow': [flow_frac.loc[(treatment, organ), (f'{cell_type}_fraction',
                                                                                            'std')]]
                                            })
                    plot_df = pd.concat([plot_df, new_row], ignore_index=True)
        # remove the rows where both fractions are 0
        plot_df = plot_df.loc[(plot_df['mean codex'] != 0) & (plot_df['mean flow'] != 0)]
        # %% make a scatterplot with mean_flow on x-axis and mean_codex on y-axis, with a circle with the width of the
        # standard deviation in both directions around each point. match circle colours to cell type
        fig, ax = plt.subplots(1, 1, figsize=(5, 5))
        sns.scatterplot(data=plot_df, x='mean flow', y='mean codex', hue='cell_type', # style='Treatment',
                        # size='Organ',
                        s=10,
                        hue_order=consensus_cell_types, ax=ax)
        # get legend colours
        handles, labels = ax.get_legend_handles_labels()
        colours = [handle.get_color() for handle in handles]
        # dict of the cell type and the corresponding colour
        cell_type_colour_dict = {cell_type: color for cell_type, color in zip(consensus_cell_types, colours)}  # [1:8]
         # add a diagonal line
        x = np.linspace(0, 0.5, 100)
        y = x
        plt.plot(x, y, 'k--')
        handles.append(plt.Line2D([0], [0], color='black', linestyle='--'))
        labels.append('diagnonal')
        ax.legend(handles=handles, labels=labels, loc='upper left',  markerscale=3, title='Consensus cell types')
        legend = ax.legend_
        legend.get_title().set_fontweight('bold')
        # bbox_to_anchor=(1.5, 1),and upper right for outside plot

        # plot an oval with the flow and codex standard deviation as the width and height
        for i, row in plot_df.iterrows():
            ax.add_patch(Ellipse((row['mean flow'], row['mean codex']), row['std flow'], row['std codex'],
                                 facecolor=cell_type_colour_dict[row['cell_type']],
                                 edgecolor=cell_type_colour_dict[row['cell_type']], lw=0.5, alpha=0.1))
        sns.scatterplot(data=plot_df, x='mean flow', y='mean codex', hue='cell_type',
                        # style='Treatment', # size='Organ',
                        s=40,
                        hue_order=consensus_cell_types, ax=ax,
                        legend=False)  # again so the points are clearer

        sns.despine()
        plt.xlabel('Flow mean fraction')
        plt.ylabel('Codex mean fraction')
        spearman_corr = stats.spearmanr(plot_df['mean flow'], plot_df['mean codex'])
        plt.savefig(f'{fig_dir}/{phenotyping_level}_{modifier}_crcl_dia1std_sp_{spearman_corr[0]:.2f}_new_plotting_style_v2.pdf')
        plt.close()
        # %% same on log scale
        fig, ax = plt.subplots(1, 1)
        sns.scatterplot(data=plot_df, x='mean flow', y='mean codex', hue='cell_type', style='Treatment',
                        size='Organ', hue_order=consensus_cell_types, ax=ax)
        # add a diagonal line
        x = np.linspace(0, 0.5, 100)
        y = x
        plt.plot(x, y, 'k--')
        # plot an oval with the flow and codex standard deviation as the width and height of the oval around each point
        for i, row in plot_df.iterrows():
            # assign colour to cell type and plot oval
            ax.add_patch(Ellipse((row['mean flow'], row['mean codex']), row['std flow'], row['std codex'],
                                 facecolor=cell_type_colour_dict[row['cell_type']],
                                 edgecolor=cell_type_colour_dict[row['cell_type']], lw=0.5, alpha=0.1))
        sns.scatterplot(data=plot_df, x='mean flow', y='mean codex', hue='cell_type', style='Treatment',
                        size='Organ', hue_order=consensus_cell_types, ax=ax, legend=False)


        plt.yscale('log')
        plt.xscale('log')
        # plt.title(f'{phenotyping_level} {modifier}')
        plt.xlabel('log(Flow mean fraction)')
        plt.ylabel('log(Codex mean fraction)')
        # add log pearson correlation
        # pearson_corr_log = stats.pearsonr(np.log(plot_df['mean flow']), np.log(plot_df['mean codex']))
        spearman_corr_log = stats.spearmanr(np.log(plot_df['mean flow']), np.log(plot_df['mean codex']))
        # add pearson to plot
        # plt.text(0.001, 0.2, f'Pearson correlation: {pearson_corr_log[0]:.2f}', fontsize=12)
        # plt.text(0.001, 0.2, f'Spearman correlation: {spearman_corr_log[0]:.2f}')
        sns.despine()
        # add legend outside, next to plot (on the right)
        ax.legend(handles=handles, labels=labels, loc='upper right', bbox_to_anchor=(1.4, 1))
        plt.savefig(f'{fig_dir}/{phenotyping_level}_{modifier}_crcl_dia1std_log_sp_{spearman_corr_log[0]:.2f}_new_plotting_style.pdf')
        # plt.show()
