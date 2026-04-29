#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
# '''=================================================
# @Project: ACT_melanoma_neutrophil (data: act_neutrophil_data_repository)
# @Author : Gemma van der Voort
# @Desc   : paired plot but without the brLNr for less crowding.
# @Desc updated: clLN/tdLN relative no-CpG analysis and plots.
# update 8.8: us english, drop useless x axis titles, add single paired d7 plot
# update 2026.01.12: reviewer comment
# '''=================================================
import math
import os
import sys
from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
from matplotlib.transforms import ScaledTranslation
import numpy as np

PROJECT_ROOT = Path(__file__).resolve()
while PROJECT_ROOT.name != "ACT_melanoma_neutrophil" and PROJECT_ROOT.parent != PROJECT_ROOT:
    PROJECT_ROOT = PROJECT_ROOT.parent
if PROJECT_ROOT.name != "ACT_melanoma_neutrophil":
    raise RuntimeError("ACT_melanoma_neutrophil root not found")
sys.path.append(str(PROJECT_ROOT / "helper_files"))

from paths_parameters import x_axis_tick_labels_flow, bar_colours_hex_blood_flow, lymph_nodes_label, \
    marker_per_organ_dict, neutrophil_label, data_repo_path
from functions_stat_test import ratio_paired_test, unpaired_multiple_conditions, unpaired_test
from functions_vis import boxplot_paired_no_brLNr, single_paired_plot, boxplot_unpaired

plt.style.use(str(PROJECT_ROOT / "helper_files" / "journal_style.mplstyle"))

# paths etc
processed_data_path = f'{data_repo_path}/processed/flow_cytometry'
stats_dir = f'{processed_data_path}/statistics'
os.makedirs(stats_dir, exist_ok=True)
fig_path = f'{data_repo_path}/figures/flow_cytometry/no_cpg_conditions/relative'
ln_df_raw = pd.read_csv(f'{processed_data_path}/ln_pan_immune_normalised_v2_cd45_corrected.csv')
modifier = 'cd45_corrected_cl_div_td'  # divide cl by td.
mtc = False # single test per cell type, so no multiple testing correction needed
force_parametric = False
non_parametric = False  # set after initial normalicy testing
os.makedirs(fig_path, exist_ok=True)

#%%
ln_df = ln_df_raw.copy()
df_condition_names = ln_df['Condition'].unique().tolist()
#%% only doing the CpG vs no CpG comparison for ACT-d14 here
condition_1 = 'ACT-d14'
condition_2 = 'ACT-d14_no CpG'
colours_conditions = ["#fa815fff","#fa815fff"] # orange, orange
x_axis_tick_labels = ['ACT\nday 14',  'ACT\nday 14\nno CpG']  # preserve order!

value_cols = ['DCs percent', 'CD11b+ Dcs percent', 'Xcr1+ Dcs percent',
              'CD8 T percent', 'Endo CD8 T percent', 'Pmel T percent',
              'CD4 Tconv percent', 'CD25+ CD4+ T percent',
              'Neutrophils percent', 'Nk cells percent', 'B cells percent']
value_cols_y_label = ['DCs\n% of CD45.2$^+$ cells', 'CD11b$^+$ DCs\n% of CD45.2$^+$ cells', 'Xcr1$^+$ DCs\n% of CD45.2$^+$ cells',
                        'CD8 T\n% of CD45.2$^+$ cells', 'Endo CD8 T\n% of CD45.2$^+$ cells', 'Pmel CD8$^+$ T cells\n% of CD45.2$^+$ cells',
                        'Conventional CD4 T\n% of CD45.2$^+$ cells', 'CD25$^+$ CD4$^+$ T\n% of CD45.2$^+$ cells',
                        'clLN/tdLN Neutrophil ratio', 'NK cells\n% of CD45.2$^+$ cells', 'B cells\n% of CD45.2$^+$ cells'
                      ]  # plural except when ending in T (is short for T cells, but that takes too much space)
modifier_2 = 'percent'
# alternative: run with counts (remember to rename downstream)
# value_cols = ['DCs | Count normalised', 'CD11b+ Dcs | Count normalised', 'Xcr1+ Dcs | Count normalised',
#                 'CD8 T | Count normalised', 'Endo CD8 T | Count normalised', 'Pmel T  | Count normalised',
#                 'CD4 Tconv | Count normalised', 'CD25+ CD4+ T | Count normalised',
#                 'Neutrophils | Count normalised', 'Nk cells | Count normalised', 'B cells | Count normalised']
# value_cols_y_label = ['DCs\nCount normalised', 'CD11b$^+$ DCs\nCount normalised', 'Xcr1$^+$ DCs\nCount normalised',
#                         'CD8 T\nCount normalised', 'Endogenous CD8 T\nCount normalised', 'Pmel CD8$^+$ T cells\nCount normalised',
#                         'Conventional CD4 T\nCount normalised', 'CD25$^+$ CD4$^+$ T\nCount normalised',
#                         'Neutrophils\nCount normalised', 'NK cells\nCount normalised', 'B cells\nCount normalised']
# modifier_2 = 'counts'
# %% update: in a new df, create relative changes per value col by dividing clLN by tdLN per mouse and condition.
ln_df_relative = ln_df.copy()
for value_col in value_cols:
    ln_df_relative[f'{value_col}_relative'] = np.nan  # initialize column
    ln_df_relative[f'{value_col}_log_relative'] = np.nan  # initialize column
    for condition in df_condition_names:
        # get tdLN and clLN data for this condition
        tdLN_data = ln_df[(ln_df['LN_type'] == 'inLNr') & (ln_df['Condition'] == condition)]
        clLN_data = ln_df[(ln_df['LN_type'] == 'inLNl') & (ln_df['Condition'] == condition)]
        # for each mouse in tdLN_data, find the corresponding clLN value and divide
        for mouse_id in tdLN_data['Mouse_ID'].unique():
            tdLN_value = tdLN_data[tdLN_data['Mouse_ID'] == mouse_id][value_col].values
            clLN_value = clLN_data[clLN_data['Mouse_ID'] == mouse_id][value_col].values
            if len(tdLN_value) > 0 and len(clLN_value) > 0:
                relative_value = clLN_value[0] / tdLN_value[0] if tdLN_value[0] != 0 else np.nan
                log_relative_value = np.log(relative_value) if relative_value > 0 else np.nan
                ln_df_relative.loc[(ln_df_relative['Mouse_ID'] == mouse_id) &
                                   (ln_df_relative['Condition'] == condition) &
                                   (ln_df_relative['LN_type'] == 'inLNl'), f'{value_col}_relative'] = relative_value
                # also store log-relative value if needed
                ln_df_relative.loc[(ln_df_relative['Mouse_ID'] == mouse_id) &
                                   (ln_df_relative['Condition'] == condition) &
                                   (ln_df_relative['LN_type'] == 'inLNl'), f'{value_col}_log_relative'] = log_relative_value
            else:
                # if no corresponding value, leave as nan
                pass
# now drop all rows that are not inLNl, since only these have the relative values
ln_df_relative = ln_df_relative[ln_df_relative['LN_type'] == 'inLNl']
ln_df_relative['LN_type'] = 'clLN/tdLN' # rename the inLNl entries to clLN/tdLN for clarity
# drop rows with nan in any of the relative value columns
for value_col in value_cols:
    ln_df_relative = ln_df_relative[~ln_df_relative[f'{value_col}_relative'].isna()]
#%% statistical test of log-relative values
# value_col = 'Neutrophils percent_log_relative'
stats_df_columns = ['Cell type', 'LN_type', 'Condition 1', 'Condition 2', 'p-value', 'test_used']
stats_df = pd.DataFrame(columns=stats_df_columns)
conditions_to_compare = [(condition_1, condition_2)]
group_col = 'Condition'
for value_col in value_cols:
    value_col_log_rel = f'{value_col}_log_relative'
    print(f'Processing {value_col_log_rel}')
    # compare 2 conditions for the log relative values
    for condition1, condition2 in conditions_to_compare:
        # t_stat, p_val, test_used, mtc_res = unpaired_multiple_conditions(ln_df_relative, value_col_log_rel, group_col,
        #                                                                  (condition1, condition2),
        #                                                                     non_parametric=non_parametric,
        #                                                                     multiple_testing_correction=mtc,
        #                                                                     force_parametric=force_parametric
        #                                                                  )
        t_stat, p_val, test_used = unpaired_test(ln_df_relative, value_col_log_rel, group_col,
                                                 (condition1, condition2),
                                                 non_parametric=non_parametric,
                                                 multiple_testing_correction=mtc
                                                 )
        print(p_val)
        stats_df = pd.concat([stats_df, pd.DataFrame([[value_col_log_rel, 'clLN/tdLN', condition1, condition2, p_val,
                                                       test_used]],
                                                     columns=stats_df_columns)])
# save stats df
stats_df.to_csv(f'{stats_dir}/pan_immune_log_cl_div_td_unpaired_stats_{condition_1}_vs_{condition_2}.csv',
                index=False)
#%% plotting
LN = 'clLN/tdLN'
# limit to the two conditions of interest
data = ln_df_relative[ln_df_relative['Condition'].isin([condition_1, condition_2])]
for cell_type_i, cell_type in enumerate(value_cols):
    value_col_log_rel = f'{cell_type}_log_relative'
    print(f'Plotting unpaired clLN/tdLN log-relative for {cell_type}')
    fig, ax = plt.subplots(1, 1, figsize=(4, 6))
    # get stats df for cell type and LN
    stats_df_LN = stats_df[(stats_df['Cell type'] == value_col_log_rel) & (stats_df['LN_type'] == LN)]
    print('1')
    print(stats_df_LN)
    stats_df_LN = stats_df_LN.drop(['Cell type', 'LN_type'], axis=1)  # drop the LN_type column for the plot
    if mtc:
        p_val_col_name = 'mtc_res'
    else:
        p_val_col_name = 'p-value'
    boxplot_unpaired(data, group_col, value_col_log_rel, stats_df=stats_df_LN, ax=ax, p_val_col_name=p_val_col_name,
                     strip_kwargs = {'palette':colours_conditions, 'marker':'d'}, # y_lim=max_y,
                     # significance_bar_order=['low', 'low', 'high', 'higher']
                     )
    ax.set_title('')  # remove title
    ax.set_ylabel(f'{value_cols_y_label[cell_type_i]}')
    # ax.set_ylabel(y_label)
    ax.set_xlabel("")
    ax.set_xticklabels(x_axis_tick_labels)
    plt.tight_layout()
    plt.savefig(f'{fig_path}/flow_{value_col_log_rel}_unpaired_{condition_1}_vs_{condition_2}_{modifier}_plot.pdf',
                dpi=300, bbox_inches='tight')
    plt.close()
