#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
# '''=================================================
# @Project: ACT_melanoma_neutrophil (data: act_neutrophil_data_repository)
# @Author : Gemma van der Voort
# @Desc   : adds relevant metadata and filters out redundant entries/columns from the phenotyping data.
# '''=================================================
import anndata as ad
import pandas as pd
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve()
while PROJECT_ROOT.name != "ACT_melanoma_neutrophil" and PROJECT_ROOT.parent != PROJECT_ROOT:
    PROJECT_ROOT = PROJECT_ROOT.parent
if PROJECT_ROOT.name != "ACT_melanoma_neutrophil":
    raise RuntimeError("ACT_melanoma_neutrophil root not found")
sys.path.append(str(PROJECT_ROOT / "helper_files"))

from paths_parameters import processed_data_dir, modifier_base, modifier_finetuning, modifier_postprocessing, possible_cell_types_codex


# %%
modifier = f'{modifier_base}_{modifier_finetuning}'
phenotyping_column = 'cell_type_round8'
full_adata = ad.read_h5ad(f"{processed_data_dir}/phenotyping_{modifier}.h5ad")
copy_adata = full_adata.copy()


metadata_sheet = pd.read_csv(f'{processed_data_dir}/codex_simplified_metadata.csv')
# replace the 4th character with a dot and the 6th with a dash, remove _CIM from experiment names
metadata_sheet['experiment'] = (
            metadata_sheet['experiment'].str[:3] + '.' + metadata_sheet['experiment'].str[4] + '-'
            + metadata_sheet['experiment'].str[6:-4])

# post-processing
markers = list(full_adata.var_names)[:-1]

# add a column 'Organ' to the adata.obs, consisting of the word after the last _ in the 'dataset_name' column
full_adata.obs['Organ'] = full_adata.obs['dataset_name'].str.split('_').str[-1]
# renaming for consistency
full_adata.obs['Organ'] = full_adata.obs['Organ'].replace('inLNR', 'inLNr')
full_adata.obs['Organ'] = full_adata.obs['Organ'].replace('S', 'Spleen')
# remove entries with organ 'skip' from the adata
full_adata = full_adata[full_adata.obs['Organ'] != 'skip']

# add a column 'Treatment' to adata.obs, which finds the cohort from the metadata sheet that matches the exp_name in
# the 'cohort' column from metadata sheet
if full_adata.obs['exp_name'].str.startswith('ME').any():
    full_adata.obs['exp_name'] = full_adata.obs['exp_name'].str[2:]
# %% replace '100.5-2_rpt' (in brLNr) with '100.5-2' in exp_name
full_adata.obs['exp_name'] = full_adata.obs['exp_name'].str.replace('100.5-2_rpt', '100.5-2')
full_adata.obs['Treatment'] = full_adata.obs['exp_name'].map(metadata_sheet.set_index('experiment')['Cohort'])
# print the 'exp_name' columns where the 'Treatment' column is NaN
print(full_adata.obs[full_adata.obs['Treatment'].isna()]['exp_name'].unique())

# add column 'population' to adata.obs, which is 'endogenous' if 'phenotyping_column' starts with 'endo', 'transferred' if
# 'phenotyping_column' starts with 'trans' and 'other' if neither.
full_adata.obs['Population'] = 'other'
full_adata.obs.loc[full_adata.obs[phenotyping_column].str.startswith('endo'), 'Population'] = 'endogenous'
full_adata.obs.loc[full_adata.obs[phenotyping_column].str.startswith('trans'), 'Population'] = 'transferred'

# %% set new columns to category type
new_columns = ['Treatment', 'Organ', 'Population']
for column in new_columns:
    full_adata.obs[column] = full_adata.obs[column].astype('category')

# rename the too long phenotyping_column names to shorter versions
full_adata.obs[phenotyping_column] = full_adata.obs[phenotyping_column].replace('CD169+ Macrophage', 'CD169+ M')

# rename phenotyping_column to 'cell_type'
full_adata.obs['Cell type'] = full_adata.obs[phenotyping_column]
full_adata.obs = full_adata.obs.drop(phenotyping_column, axis=1)

# rename '100_9_1_CIM-reg_3x3_HOE_2275_brLNr' to 100_9_1_CIM-reg_3x3_HOE-2275_brLNr (and others)
full_adata.obs['dataset_name'] = full_adata.obs['dataset_name'].str.replace('HOE_', 'HOE-')
full_adata.obs['dataset_name'] = full_adata.obs['dataset_name'].str.replace('BAL_', 'BAL-')
# create columns mouse_id in adata.obs. from '100.2-1_CIM-reg_2x2_BAL-3614_brLNr', extract the BAL-3614 part.
full_adata.obs['mouse_id'] = full_adata.obs['dataset_name'].str.split('_').str[-2].astype('category')

# check if all data was added correctly:
print(full_adata)
# %% print all unique values of the obs variables in the full_adata
for variable in list(full_adata.obs.columns):
    print(f'{variable}: {full_adata.obs[variable].unique()}')

# %% save post-processed data without further filtering
# full_adata.write(f"{processed_data_dir}/phenotyping_{modifier}_incl_metadata.h5ad")

# # subsample 0.5% of the data for testing
# full_adata = sc.pp.subsample(full_adata, fraction=0.005, copy=True)
# full_adata.write(f"{processed_data_dir}/phenotyping_{modifier}_subsampled_0005.h5ad")
#
# # subsample 5% of the data for testing
# full_adata = sc.pp.subsample(full_adata, fraction=0.05, copy=True)
# full_adata.write(f"{processed_data_dir}/phenotyping_{modifier}_subsampled_005.h5ad")

# %%list the var names
print(full_adata.var_names)

# comparison
# rename 'CD169+ Macrophage' to 'CD169+ M' in all_cell_types
all_cell_types = [cell_type.replace('CD169+ Macrophage', 'CD169+ M') for cell_type in possible_cell_types_codex]
print(set(all_cell_types) - set(full_adata.obs['Cell type'].unique()))
# print amounts of each phenotyping_column
print(full_adata.obs['Cell type'].value_counts())

# %%check on amount of trans (and subtype) T cells in the control samples and remove them from the adata
# per treatment, print percentage of population that is trans vs all cells in treatment
trans_frac_dict = {}
for treatment in full_adata.obs['Treatment'].unique():
    print(treatment)
    print(
        f"""{full_adata.obs[(full_adata.obs['Treatment'] == treatment) & (full_adata.obs['Population'] == 'transferred')]
             ['Population'].count() / full_adata.obs[full_adata.obs['Treatment'] == treatment]['Population'].count():5f}""")
    trans_frac_dict[treatment] = \
    full_adata.obs[(full_adata.obs['Treatment'] == treatment) & (full_adata.obs['Population'] ==
                                                                 'transferred')][
        'Population'].count() / full_adata.obs[full_adata.obs['Treatment'] == treatment][
        'Population'].count()

# %% now check how many of the 3 transferred cell types (trans T cycling , trans T exhausted, trans CD8) are in the
# control treatments
print('transferred populations in control treatments')
print(full_adata.obs[full_adata.obs['Population'] == 'transferred']['Treatment'].value_counts())
print('trans T cycling')
print(full_adata.obs[full_adata.obs['Cell type'] == 'trans T cycling']['Treatment'].value_counts())
print('trans T exhausted')
print(full_adata.obs[full_adata.obs['Cell type'] == 'trans T exhausted']['Treatment'].value_counts())
print('trans CD8')
print(full_adata.obs[full_adata.obs['Cell type'] == 'trans CD8']['Treatment'].value_counts())

# these are all acceptable background levels.
# %%remove the transferred cells from the treatments tumour_bearing and Cyclo_d1
full_adata = full_adata[~((full_adata.obs['Treatment'] == 'tumour_bearing') & (full_adata.obs['Population'] ==
                                                                               'transferred'))]
full_adata = full_adata[~((full_adata.obs['Treatment'] == 'Cyclo_d1') & (full_adata.obs['Population'] ==
                                                                         'transferred'))]
# %% prevent typeError in saving due to 'root' column in adata.uns['codex_xml'], which is an 'Element Study' object
# from lxml.etree._Element, I think the id and location of the first 2 columns should be enough.
# remove the 'root' column from adata.uns['codex_xml']:
# full_adata.uns['codex_xml'] = full_adata.uns['codex_xml'][['exp_id', 'path']]
# this seems to be cropping up in other .uns columns as well. take the chance to slim down the adata. I will remove
# all 'uns' and all 'layers'. I will keep the 'obsm' and 'var' properties, as they are useful for down stream analysis.
full_adata.uns = {}
full_adata.layers = {}
# in case there are cells with umap coordinates beyond -10 and 10, remove them
# this is also done in during fine-tuning, so the umap is more informative there. Keep in case the data is not
# fine-tuned.
print(f'length of adata before umap-based filtering: {len(full_adata)}')
full_adata = full_adata[full_adata.obsm['X_umap'][:, 0] < 10]
full_adata = full_adata[full_adata.obsm['X_umap'][:, 0] > -10]
full_adata = full_adata[full_adata.obsm['X_umap'][:, 1] < 10]
full_adata = full_adata[full_adata.obsm['X_umap'][:, 1] > -10]
print(f'length of adata after umap-based filtering: {len(full_adata)}')
# %%save
full_adata.write(f"{processed_data_dir}/phenotyping_{modifier_base}_{modifier_finetuning}_{modifier_postprocessing}.h5ad")
# from 4gb to 900MB, great!
# %% make a summary document of this data. include the following:
# - amount of files and their names (dataset name unique)
# - number of mice per treatment
# - amount and percentage of trans cells in each Treatment

# map mouse_id to treatment
mouse_treatment_mapping = full_adata.obs[['mouse_id', 'Treatment']].drop_duplicates()
# how many mice per treatment
mice_per_treatment = mouse_treatment_mapping['Treatment'].value_counts()

with open(f'{processed_data_dir}/phenotyping_{modifier_base}_{modifier_finetuning}_{modifier_postprocessing}'
          f'_summary.txt', 'w') as f:
    f.write(f'Amount of files: {len(full_adata.obs["dataset_name"].unique())}\n\n')
    f.write(f'Names of files: {list(full_adata.obs["dataset_name"].unique())}\n\n')
    f.write(f'Number of segmented cells per treatment: \n{full_adata.obs["Treatment"].value_counts()}\n\n')
    f.write(f'Number of mice per treatment: \n{mice_per_treatment}\n\n')
    f.write(f'Percentage of transferred cells in each treatment:\n')
    for treatment in full_adata.obs['Treatment'].unique():
        f.write(f"""{treatment}: {full_adata.obs[(full_adata.obs["Treatment"] == treatment) &
                                                 (full_adata.obs["Population"] == "transferred")]["Population"].count()
                                  / full_adata.obs[full_adata.obs["Treatment"] == treatment]["Population"].count():5f}\n""")
    f.write(f'\nAmount of each cell type:\n{full_adata.obs["Cell type"].value_counts()}')
    # add missing cell types
    f.write(f'\nMissing cell types: {set(all_cell_types) - set(full_adata.obs["Cell type"].unique())}')

# # %% we might have some mixed cell-types. set all 'Cell types' that contain a '/' to 'mixed'
# full_adata.obs['Cell type'] = full_adata.obs['Cell type'].str.replace('.*/.*', 'mixed')
