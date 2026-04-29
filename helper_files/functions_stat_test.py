#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
# '''=================================================
# @Project: ACT_melanoma_neutrophil (data: act_neutrophil_data_repository)
# @Author : Gemma van der Voort
# @Desc   : script for testing significance. Make into functions later.
# @Desc updated: Statistical testing helper functions.
# '''=================================================
import pandas as pd
from scipy import stats
import warnings
import numpy as np
import pandas as pd
import numpy as np
from scipy import stats
from statsmodels.stats.multicomp import pairwise_tukeyhsd
from statsmodels.stats.multitest import multipletests
from scipy.stats import kruskal
from scikit_posthocs import posthoc_dunn
warnings.filterwarnings("ignore")


def paired_test(df_1, df_2, value_col, non_parametric=False):
    """
    Perform a paired t-test for the values in value_col for group1 and group2 in the df.
    :param df_1: dataframe with the data for group1
    :param df_2: dataframe with the data for group2
    :param value_col: column with the values
    :param non_parametric: if True, use a non-parametric test
    :return: t-statistic, p-value, test_used
    """
    group1_values = df_1[value_col]
    group2_values = df_2[value_col]
    # do shapiro test for normality
    stat1, p1 = stats.shapiro(group1_values)
    stat2, p2 = stats.shapiro(group2_values)
    # print(f'Normality test for {group1}: {p1}, for {group2}: {p2}')
    if p1 < 0.05 or p2 < 0.05 or non_parametric:
        # print(f'Normality test failed for {group1} or {group2}, use a non-parametric test')
        print(f'Normality test failed, use a non-parametric test')
        test_used = "Wilcoxon"
        result = stats.wilcoxon(group1_values, group2_values)
        t_stat, p_val = result.statistic, result.pvalue
    else:
        t_stat, p_val = stats.ttest_rel(group1_values, group2_values)
        test_used = "t-test"
    return t_stat, p_val, test_used

def ratio_paired_test(df_1, df_2, value_col, non_parametric=False, force_parametric=False):
    """
    Perform a ratio paired t-test for the values in value_col for group1 and group2 in the df.
    If force_parametric=True, always use the parametric test (t-test),
    even if normality fails.
    :param df_1:
    :param df_2:
    :param value_col:
    :param non_parametric:
    :param force_parametric:
    :return:
    """
    group1_values = np.array(df_1[value_col])
    group2_values = np.array(df_2[value_col])
    # assert that the values are not 0 or negative (since otherwise the test cannot be performed)
    if (group1_values <= 0).any() or (group2_values <= 0).any():
        print(f'Values in {value_col} for group1 or group2 contain 0 or negative values, cannot perform ratio test')
        return None, None, None

    ratios = group1_values / group2_values
    log_ratios = np.log(ratios)

    # Only test normality if not forced parametric
    if not force_parametric:
        stat, p = stats.shapiro(log_ratios)
        if p < 0.05 or non_parametric:
            print('Normality test failed or non_parametric=True: using Wilcoxon signed-rank test')
            test_used = "Wilcoxon signed-rank test"
            result = stats.wilcoxon(group1_values, group2_values)
            t_stat, p_val = result.statistic, result.pvalue
            return t_stat, p_val, test_used

    # Parametric path (either normal or forced)
    print('Using parametric ratio paired t-test')
    t_stat, p_val = stats.ttest_1samp(log_ratios, 0)
    test_used = "ratio paired t-test"
    return t_stat, p_val, test_used


def unpaired_test(df, value_col, group_col, groups, non_parametric=False, multiple_testing_correction=False):
    """
    Perform an unpaired t-test for the values in value_col for group1 and group2 in the df.
    :param df: dataframe with the data
    :param value_col: column with the values (continuous)
    :param group_col: column with the groups (categorical)
    :param groups: tuple with the two groups
    :param non_parametric: if True, use a non-parametric test
    :param multiple_testing_correction: if True, apply multiple testing correction
    :return: t-statistic, p-value
    """
    group1, group2 = groups
    group1_values = df[df[group_col] == group1][value_col]
    # print(f'Group1 values ({group1}): {group1_values}')
    group2_values = df[df[group_col] == group2][value_col]
    # print(f'Group2 values ({group2}): {group2_values}')
    # do shapiro test for normality
    stat1, p1 = stats.shapiro(group1_values)
    stat2, p2 = stats.shapiro(group2_values)
    if p1 < 0.05 or p2 < 0.05 or non_parametric:
        print(f'Normality test failed for {group1} or {group2}, use a non-parametric test')
        test_used = "Mann-Whitney"
        t_stat, p_val = stats.mannwhitneyu(group1_values, group2_values)
    else:
        t_stat, p_val = stats.ttest_ind(group1_values, group2_values)
        test_used = "t-test"
    if multiple_testing_correction:
        print('Multiple testing correction not implemented yet')
        pass
    return t_stat, p_val, test_used


def unpaired_multiple_conditions(df, value_col, group_col, groups, non_parametric=False, multiple_testing_correction=False,
                                 force_parametric=False
                                 ):
    """
    Perform a one-way ANOVA for the values in value_col for the groups in the df, or its non-parametric equivalent
    (Kruskal-Wallis).
    If multiple_testing_correction is True, apply post-hoc tests.
    :param df: dataframe with the data
    :param value_col: column with the values (continuous)
    :param group_col: column with the groups (categorical)
    :param groups: list with the groups # do all groups in the df
    :param force_parametric: force the use of parametric one-way anova. If biological assumptions permit.
    :param non_parametric: if True, use a non-parametric test. mutually exclusive with force_parametric.
    :param multiple_testing_correction: if True, apply multiple testing correction
    :return: f-statistic, p-value, and if MTC is applied, also the results of the post-hoc tests
    """
    group_values = [df[df[group_col] == group][value_col] for group in groups]

    # do shapiro test for normality
    p_values = [stats.shapiro(group_values[i])[1] for i in range(len(group_values))]

    if (any([p < 0.05 for p in p_values]) or non_parametric) and not force_parametric:
        # print(f'Normality test failed for one of the groups, use a non-parametric test')
        test_used = "Kruskal-Wallis"
        f_stat, p_val = stats.kruskal(*group_values)

        if multiple_testing_correction:
            # Apply Dunn's test for multiple comparisons
            mtc_dunn = 'holm'  # 'bonferroni', 'holm', 'fdr_bh', 'fdr_by'
            dunn_result = posthoc_dunn(df, val_col=value_col, group_col=group_col, p_adjust=mtc_dunn)
            col_group = groups[0]  # symmetrical matrix, so row/col order does not matter. but for clarity
            row_group = groups[1]
            dunn_result = dunn_result.loc[col_group, row_group]
            print(f'dunn result: {dunn_result}')
            test_used = f"Kruskal-Wallis, Dunn's test with {mtc_dunn} for mtc"

            return f_stat, p_val, test_used, dunn_result
        else:
            return f_stat, p_val, test_used, 'no mtc'
    else:
        f_stat, p_val = stats.f_oneway(*group_values)
        test_used = "ANOVA"
        if multiple_testing_correction:
            # Apply Tukey's HSD for multiple comparisons
            tukey_result = pairwise_tukeyhsd(endog=df[value_col], groups=df[group_col], alpha=0.05)
            group1 = groups[0]  # non-symmetrical
            group2 = groups[1]
            # Convert the SimpleTable to a DataFrame
            tukey_df = tukey_result._results_table.data
            test_used = f"ANOVA with Tukey's HSD for mtc"
            # Convert the DataFrame-like data into an actual DataFrame
            tukey_df = pd.DataFrame(tukey_df[1:], columns=tukey_df[0])
            for index, row in tukey_df.iterrows():
                if row['group1'] == group1 and row['group2'] == group2:
                    naive_tumour_diff = row['meandiff']
                    tukey_p_adj = row['p-adj']
                    break
                elif row['group1'] == group2 and row['group2'] == group1:
                    tukey_p_adj = row['p-adj']
                else:
                    pass
            return f_stat, p_val, test_used, tukey_p_adj
        else:
            return f_stat, p_val, test_used, 'no mtc'
