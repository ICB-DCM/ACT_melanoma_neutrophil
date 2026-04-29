#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
# '''=================================================
# @Project: ACT_melanoma_neutrophil (data: act_neutrophil_data_repository)
# @Author : Gemma van der Voort
# @Desc   : function to make a boxplot with unpaired data, in the style as developed in plotting_style_exploration_2.py
# @Desc updated: Visualization helper functions for flow/codex plots.
# '''=================================================
import seaborn as sns
import matplotlib.pyplot as plt
import math
import numpy as np
import pandas as pd
from scipy.special.cython_special import log_ndtr


def p_val_formatting(p_val: float):
    if p_val < 0.0001:
        stars = '****'
        text_fontsize = 20
        text_height_modifier = 0
    elif p_val < 0.001:
        stars = '***'
        text_fontsize = 20
        text_height_modifier = 0
    elif p_val < 0.01:
        stars = '**'
        text_fontsize = 20
        text_height_modifier = 0
    elif p_val < 0.05:
        stars = '*'
        text_fontsize = 20
        text_height_modifier = 0
    else:
        stars = 'ns'
        text_fontsize = 16
        text_height_modifier = 0.03
    return stars, text_fontsize, text_height_modifier

def boxplot_unpaired(data, x_col, y_col, stats_df=None, significance_bar_order=None, ax=None, y_lim=None,
                    p_val_col_name='p-value',
                     box_kwargs=None, strip_kwargs=None):
    """
    Make a boxplot + scatterplot with unpaired data, in the style as developed in plotting_style_exploration_2.py.
    If a statistical test is wanted, provide a stats_df of the form (condition 1, condition 2, p-value).
    Per default, the conditions for which test results are provided are plotted alternating low and high
    in the order they are provided in the stats_df. If a custom order is provided, the significance bars are plotted
    in that order.

    :param data: dataframe with the data
    :param x_col: column name of the x_col-axis
    :param y_col: column name of the y-axis
    :param stats_df: dataframe with the statistical test results
    :param significance_bar_order: list with the order in which to plot the significance bars
    :param ax: axis to plot on
    :param y_lim: limit for the y-axis, used for proportioning of the significance bars. If None, ceil of the
    highest value in the y-axis is used.
    :param p_val_col_name: name of the column in stats_df with the p-values
    :param box_kwargs: kwargs for the boxplot
    :param strip_kwargs: kwargs for the stripplot
    :return: a matplotlib plot

    Parameters
    ----------

    """
    if strip_kwargs is None:
        strip_kwargs = {}
    if box_kwargs is None:
        box_kwargs = {}
    if ax is None:
        ax = plt.gca()
    if y_lim is None:
        y_lim = math.ceil(data[y_col].max())
    if stats_df is not None and significance_bar_order is None:
        significance_bar_order = ['low', 'high'] * (len(stats_df) // 2) # even case
        significance_bar_order = significance_bar_order + ['low'] * (len(stats_df) % 2)  # for uneven case
    # boxplot
    sns.boxplot(x=x_col, y=y_col, data=data, ax=ax,
                fliersize=0,  # no outliers
                whis=0,  # no whiskers
                showcaps=False,  # no caps on the whiskers
                boxprops={'facecolor': 'none', 'linewidth': 0},  # no fill in the box
                medianprops={'color': 'black', 'linewidth': 2},  # median line emphasise
                # add a line for the 75 and 25 percentile

                width=0.5,  # adjust width of the box to be more narrow
                **box_kwargs
                )
    # Compute percentiles
    categories = data[x_col].unique()
    for i, category in enumerate(categories):
        subset = data[data[x_col] == category][y_col]

        # Compute Q1 (25th percentile) and Q3 (75th percentile)
        q1 = np.percentile(subset, 25)
        q3 = np.percentile(subset, 75)

        # Overlay thicker lines at Q1 and Q3
        ax.hlines(y=[q1, q3], xmin=i - 0.35 / 2, xmax=i + 0.35 / 2, color='black')
        # Add a vertical line between Q1 and Q3
        ax.vlines(x=i, ymin=q1, ymax=q3, color='black')

    sns.stripplot(x=x_col, y=y_col, data=data, ax=ax,
                  # palette=bar_colours_hex,  # should pass using the
                  jitter=0.2, size=10, edgecolor='black' , linewidth=1,
                  zorder=0,  # move to back (request Jan)
                  # marker='s',
                  **strip_kwargs
                  )
    sns.despine()
    if stats_df is not None:
        # get the 1st and second column and make a list of tuples
        conditions_to_compare = [(stats_df.iloc[i, 0], stats_df.iloc[i, 1]) for i in range(len(stats_df))]
        # p_vals_list = stats_df.iloc[:, 2].tolist()  # todo modify here
        p_vals_list = stats_df.loc[:, p_val_col_name].tolist()
        print(p_vals_list)
        for i_pair, cond_pair in enumerate(conditions_to_compare):
            p_val = p_vals_list[i_pair]
            print(p_val)
            # get different number of stars for different p values
            stars, text_fontsize, text_height_modifier = p_val_formatting(p_val)
            # get the y value of the highest point in the plot
            y_highest_point = data[y_col].max()  # for line placement correction
            h1 = data[data[x_col] == cond_pair[0]][y_col].max()
            h2 = data[data[x_col] == cond_pair[1]][y_col].max()
            y_horizontal_bar_low = y_highest_point + (y_lim * 0.07)  # a bit higher to be more visually appealing
            y_horizontal_bar_high = y_highest_point + (y_lim * 0.27)
            y_horizontal_bar_higher = y_highest_point + (y_lim * 0.47)  # a bit higher to be more visually appealing
            if significance_bar_order[i_pair] == 'low':
                h1 = h1 + (y_lim * 0.05)  # a bit higher to be more visually appealing
                h2 = h2 + (y_lim * 0.05)  # a bit higher to be more visually appealing
                y_max_text = y_horizontal_bar_low + (y_highest_point * text_height_modifier)  # position of the text
                y = y_horizontal_bar_low
            elif significance_bar_order[i_pair] == 'high':
                # set the same height to top of low condition for both
                h1 = y_horizontal_bar_low + (y_lim * 0.15)  # a bit higher to be more visually appealing
                h2 = y_horizontal_bar_low + (y_lim * 0.15)  # a bit higher to be more visually appealing
                y_max_text = y_horizontal_bar_high + (y_highest_point * text_height_modifier)  # position of the text
                y = y_horizontal_bar_high
            elif significance_bar_order[i_pair] == 'higher':
                # set the same height to top of low condition for both
                h1 = y_horizontal_bar_high + (y_lim * 0.15)  # a bit higher to be more visually appealing
                h2 = y_horizontal_bar_high + (y_lim * 0.15)  # a bit higher to be more visually appealing
                y_max_text = y_horizontal_bar_higher + (y_highest_point * text_height_modifier)  # position of the text
                y = y_horizontal_bar_higher
            else:
                raise ValueError(f"Unrecognized significance_bar_order value: {significance_bar_order[i_pair]}")

            # for function, add a conditional on if condition_plot_location is given, and if not, alternate.
            # unpaired: x1 is the x_col-tick location of cond_pair[0], x2 is the x_col-tick location of cond_pair[1]
            x1 = data[x_col].unique().tolist().index(cond_pair[0])  # + 0.25
            x2 = data[x_col].unique().tolist().index(cond_pair[1])  # - 0.25
            # get mid between x1 and x2 for the x_col location of the asterisk
            x_mid = (x1 + x2) / 2
            ax.text(x_mid, y_max_text, stars, fontsize=text_fontsize, ha='center', va='bottom')
            ax.hlines(y, x1, x2, color='black')
            # for the final vline to be drawn, pull it down to the 'low' condition h2
            if significance_bar_order[i_pair] == 'high' and i_pair == len(conditions_to_compare) - 1:
                h2 = data[data[x_col] == cond_pair[1]][y_col].max()
                h2 = h2 + (y_horizontal_bar_low * 0.05)
            # check that h1 and h2 are not higher than y. if they are, don't add that vline
            if h1 < y:
                ax.vlines(x1, y, h1, color='black')  # vertical line left
            if h2 < y:
                ax.vlines(x2, y, h2, color='black')  # vertical line right
    ax.set_title(y_col)  # , y=1.05)
    ax.set_ylabel(y_col) # , fontsize=16)
    return ax


def boxplot_paired(data, x_col, y_col, x_group_col, stats_df=None, ax=None, y_lim=None,
                     x_group_order=None, markers=None, condition_colours=None, omit_ns=False,
                   box_kwargs=None, scatter_kwargs=None):
    """
    Make a boxplot + scatterplot with paired data, in the style as developed in plotting_style_exploration_2.py
    If a statistical test is wanted, provide a stats_df of the form (paired condition 1, paired condition 2, p-value,
    with length of x-axis categories).

    :param data: dataframe with the data
    :param x_col: column name of the x_col-axis
    :param y_col: column name of the y-axis
    :param x_group_col: list of tuples with the pairs to compare
    :param stats_df: dataframe with the statistical test results. if none, no significance bars are plotted
    :param ax: axis to plot on
    :param y_lim: limit for the y-axis, used for proportioning of the significance bars. If None, ceil of the
    highest value in the y-axis is used.
    :param x_group_order: order in which to plot the x_group_col. If None, the order in which the data is provided
    in x_group_col is used.
    :param markers: custom markers for the scatterplot. If None, the markers are set to the first
    len(x_group_col.unique()) markers for the following: ['o', 's', '^', 'v', 'D', 'P', 'X'], then cycles.
    :param condition_colours: colours for the conditions, should be len(x_col.unique()).
    If None, the colours are default seaborn category colours.
    :param box_kwargs: kwargs for the boxplot
    :param scatter_kwargs: kwargs for the scatterplot
    :return: a matplotlib plot
    """
    np.random.seed(0)
    if scatter_kwargs is None:
         scatter_kwargs = {}
    if box_kwargs is None:
        box_kwargs = {}
    if ax is None:
        ax = plt.gca()
    if y_lim is None:
        y_lim = math.ceil(data[y_col].max())
    if x_group_order is None:  # if no order is given, use the order in which the data is provided
        x_group_order = data[x_group_col].unique()
        # print(x_group_order)
    else:
        data[x_group_col] = pd.Categorical(data[x_group_col], categories=x_group_order, ordered=True)
    if markers is None:
        markers = ['o', 's', '^', 'v', 'D', 'P', 'X']
        markers = markers * (len(data[x_group_col].unique()) // len(markers) + 1)
    if condition_colours is None:
        condition_colours = sns.color_palette('tab20', len(data[x_col].unique()))
    # boxplot
    median_width = 0.75
    percentile_width = 0.5
    sns.boxplot(x=x_col, y=y_col, data=data, ax=ax, hue=x_group_col,  # look at what this is called!
                fliersize=0,  # no outliers
                whis=0,  # no whiskers
                showcaps=False,  # no caps on the whiskers
                boxprops={'facecolor': 'none', 'linewidth': 0},  # no fill in the box
                medianprops={'color': 'black', 'linewidth': 2},  # median line emphasise
                width=median_width,  # adjust width of the box to be more narrow,
                # order=data[x_col].cat.categories, # ensures order by catergories
                **box_kwargs
                )
    # calculate offsets based on the width of the boxplot
    x_positions = np.arange(len(data[x_col].unique()))  # positions (0-6)
    offsets = [-0.25, 0, 0.25]  # calculate based on median width?


    # Compute percentiles
    categories = data[x_col].unique()
    y_coords_list = []
    for i_cond, category in enumerate(categories):
        for i_ln, ln in enumerate(x_group_order):
            subset = data[(data[x_col] == category) & (data[x_group_col] == ln)][y_col]
            # Compute Q1 (25th percentile) and Q3 (75th percentile)
            q1 = np.percentile(subset, 25)
            q3 = np.percentile(subset, 75)
            # find the middle of the x-axis for the boxplot for this condition and LN combo
            i_cond_ln = i_cond + (i_ln - 1) * 0.25
            xmin = i_cond_ln - (percentile_width / 6)
            xmax = i_cond_ln + (percentile_width / 6)
            # Overlay thicker lines at Q1 and Q3. The boxplot is 0.75 wide, so make the lines 0.5 wide
            ax.hlines(y=[q1, q3], xmin=xmin, xmax=xmax, color='black')
            # Add a vertical line between Q1 and Q3
            ax.vlines(x=i_cond_ln, ymin=q1, ymax=q3, color='black')

            # for the stripplot: (in new loop so it is plotted on top of the boxplot)
            x_coords = np.repeat(x_positions[i_cond] + offsets[i_ln], len(subset))
            x_coords_jittered = x_coords + np.random.uniform(-0.05, 0.05, len(subset))
            y_coords = subset
            y_coords_list.append(y_coords)
            ax.scatter(x_coords_jittered, y_coords, marker=markers[i_ln], color=condition_colours[i_cond],
                       zorder=0,  # move to back (request Jan)
                       edgecolor='black', **scatter_kwargs)
    sns.despine()
    ax.get_legend().remove()  # remove the legend
    if stats_df is not None:
        # get the 1st and second column and make a list of tuples
        ln_to_compare = [(stats_df.iloc[i, 0], stats_df.iloc[i, 1]) for i in range(len(stats_df))]  # compares ln1 vs ln2 7 times, once for each condition. OK
        i_ln1 = x_group_order.index(ln_to_compare[0][0])  # i_LN1 is index of ln in LN_order
        ln1 = ln_to_compare[0][0]
        i_ln2 = x_group_order.index(ln_to_compare[0][1])  # i_LN2 is index of ln in LN_order
        ln2 = ln_to_compare[0][1]
        p_vals_list = stats_df.iloc[:, 2].tolist()  # each condition has an associated p-value. OK
        # get the y value of the highest point in the plot
        y_highest_point = data[y_col].max()  # for line placement correction
        for i_cond, condition in enumerate(data[x_col].unique()):
            # df for condition
            data_condition = data[data[x_col] == condition]
            p_val = p_vals_list[i_cond]
            # get different amount of stars for different p values
            stars, text_fontsize, text_height_modifier = p_val_formatting(p_val)

            # skip if omit_ns=True and not significant
            if omit_ns and stars.lower() == "ns":
                continue
            data_cond_ln1 = data_condition[data_condition[x_group_col] == ln1]
            data_cond_ln2 = data_condition[data_condition[x_group_col] == ln2]
            h1 = data_cond_ln1[y_col].max()
            h2 = data_cond_ln2[y_col].max()
            y_horizontal_bar_low = y_highest_point + (y_lim * 0.07)  # a bit higher to be more visually appealing
            h1 = h1 + (y_lim * 0.05)  # a bit higher to be more visually appealing
            h2 = h2 + (y_lim * 0.05)  # a bit higher to be more visually appealing
            y_max_text = y_horizontal_bar_low + (y_highest_point * text_height_modifier)  # position of the text
            y = y_horizontal_bar_low
            x1 = x_positions[i_cond] + offsets[i_ln1]
            x2 = x_positions[i_cond] + offsets[i_ln2]
            # get mid between x1 and x2 for the x_col location of the asterisk
            x_mid = float((x1 + x2) / 2)
            ax.text(x_mid, y_max_text, stars, fontsize=text_fontsize, ha='center', va='bottom')
            ax.hlines(y, x1, x2, color='black')
            # check that h1 and h2 are not higher than y. if they are, don't add that vline
            if h1 < y:
                ax.vlines(x1, y, h1, color='black')  # vertical line left
            if h2 < y:
                ax.vlines(x2, y, h2, color='black')  # vertical line right
    return ax


def boxplot_paired_no_brLNr(
    data,
    x_col,
    y_col,
    x_group_col,
    stats_df=None,
    ax=None,
    y_lim=None,
    x_group_order=None,
    markers=None,
    condition_colours=None,
    omit_ns=False,
    box_kwargs=None,
    scatter_kwargs=None,
    mouse_id_female=None,
    mouse_id_col='Mouse_ID',
):
    """
    Make a boxplot + scatterplot with paired data, in the style as developed in plotting_style_exploration_2.py
    If a statistical test is wanted, provide a stats_df of the form (paired condition 1, paired condition 2, p-value,
    with length of x-axis categories).

    Parameters
    ----------
    data : pd.DataFrame
        Dataframe with the data.
    x_col : str
        Column name of the x-axis.
    y_col : str
        Column name of the y-axis.
    x_group_col : str
        Column used for paired grouping within each x-axis category.
    stats_df : pd.DataFrame or None
        Dataframe with statistical test results. If None, no significance bars are plotted.
    ax : matplotlib.axes.Axes or None
        Axis to plot on. If None, uses current axis.
    y_lim : float or None
        Limit for the y-axis, used for proportioning of significance bars.
        If None, ceil of the highest value in the y-axis is used.
    x_group_order : list-like or None
        Order in which to plot the x_group_col categories.
        If None, the order in which they appear in the data is used.
    markers : list or None
        Custom markers for the scatterplot.
        If None, defaults to ['o', 's', '^', 'v', 'D', 'P', 'X'] and cycles as needed.
    condition_colours : list or None
        Colours for the conditions. If None, uses seaborn tab20 palette.
    omit_ns : bool
        If True, omit significance bars for non-significant p-values.
    box_kwargs : dict or None
        Extra kwargs for sns.boxplot.
    scatter_kwargs : dict or None
        Extra kwargs for ax.scatter.
    mouse_id_female : list-like or None
        Optional list of mouse IDs that should be plotted in black instead of the condition colour.
        If None, original behaviour is preserved.
    mouse_id_col : str
        Column name containing mouse IDs. Only used if mouse_id_female is not None.

    Returns
    -------
    matplotlib.axes.Axes
        The plotted axis.
    """
    np.random.seed(0)

    if scatter_kwargs is None:
        scatter_kwargs = {}
    if box_kwargs is None:
        box_kwargs = {}
    if ax is None:
        ax = plt.gca()
    if y_lim is None:
        y_lim = math.ceil(data[y_col].max())

    if x_group_order is None:
        x_group_order = list(data[x_group_col].unique())
    else:
        x_group_order = list(x_group_order)
        data = data.copy()
        data[x_group_col] = pd.Categorical(
            data[x_group_col],
            categories=x_group_order,
            ordered=True
        )

    if markers is None:
        markers = ['o', 's', '^', 'v', 'D', 'P', 'X']
        markers = markers * (len(data[x_group_col].unique()) // len(markers) + 1)

    if condition_colours is None:
        condition_colours = sns.color_palette('tab20', len(data[x_col].unique()))

    if mouse_id_female is not None:
        mouse_id_female = set(mouse_id_female)

    # boxplot
    median_width = 0.75
    percentile_width = 0.5

    sns.boxplot(
        x=x_col,
        y=y_col,
        data=data,
        ax=ax,
        hue=x_group_col,
        fliersize=0,
        whis=0,
        showcaps=False,
        boxprops={'facecolor': 'none', 'linewidth': 0},
        medianprops={'color': 'black', 'linewidth': 2},
        width=median_width,
        **box_kwargs
    )

    x_positions = np.arange(len(data[x_col].unique()))

    def calculate_offset(median_width, n, i):
        return (median_width / n) * (i - (n - 1) / 2)

    offset1 = calculate_offset(median_width, 2, 0)
    offset2 = calculate_offset(median_width, 2, 1)
    offsets = [offset1, offset2]

    categories = data[x_col].unique()
    y_coords_list = []

    for i_cond, category in enumerate(categories):
        for i_ln, ln in enumerate(x_group_order):
            subset_df = data[(data[x_col] == category) & (data[x_group_col] == ln)].copy()
            subset = subset_df[y_col]

            if len(subset) == 0:
                continue

            q1 = np.percentile(subset, 25)
            q3 = np.percentile(subset, 75)

            i_cond_ln = i_cond + offsets[i_ln]
            xmin = i_cond_ln - (percentile_width / 6)
            xmax = i_cond_ln + (percentile_width / 6)

            ax.hlines(y=[q1, q3], xmin=xmin, xmax=xmax, color='black')
            ax.vlines(x=i_cond_ln, ymin=q1, ymax=q3, color='black')

            x_coords = np.repeat(x_positions[i_cond] + offsets[i_ln], len(subset_df))
            x_coords_jittered = x_coords + np.random.uniform(-0.05, 0.05, len(subset_df))
            y_coords = subset_df[y_col].to_numpy()
            y_coords_list.append(y_coords)

            # backward-compatible default: all points use condition colour
            if mouse_id_female is None:
                ax.scatter(
                    x_coords_jittered,
                    y_coords,
                    marker=markers[i_ln],
                    color=condition_colours[i_cond],
                    zorder=0,
                    edgecolor='black',
                    **scatter_kwargs
                )
            else:
                point_colours = [
                    'black' if mouse_id in mouse_id_female else condition_colours[i_cond]
                    for mouse_id in subset_df[mouse_id_col]
                ]

                ax.scatter(
                    x_coords_jittered,
                    y_coords,
                    marker=markers[i_ln],
                    c=point_colours,
                    zorder=0,
                    edgecolor='black',
                    **scatter_kwargs
                )

    sns.despine()

    legend = ax.get_legend()
    if legend is not None:
        legend.remove()

    if stats_df is not None:
        ln_to_compare = [(stats_df.iloc[i, 0], stats_df.iloc[i, 1]) for i in range(len(stats_df))]
        i_ln1 = x_group_order.index(ln_to_compare[0][0])
        ln1 = ln_to_compare[0][0]
        i_ln2 = x_group_order.index(ln_to_compare[0][1])
        ln2 = ln_to_compare[0][1]
        p_vals_list = stats_df.iloc[:, 2].tolist()
        y_highest_point = data[y_col].max()

        for i_cond, condition in enumerate(data[x_col].unique()):
            data_condition = data[data[x_col] == condition]
            p_val = p_vals_list[i_cond]
            stars, text_fontsize, text_height_modifier = p_val_formatting(p_val)

            if omit_ns and stars.lower() == "ns":
                continue

            data_cond_ln1 = data_condition[data_condition[x_group_col] == ln1]
            data_cond_ln2 = data_condition[data_condition[x_group_col] == ln2]
            h1 = data_cond_ln1[y_col].max()
            h2 = data_cond_ln2[y_col].max()
            y_horizontal_bar_low = y_highest_point + (y_lim * 0.07)
            h1 = h1 + (y_lim * 0.05)
            h2 = h2 + (y_lim * 0.05)
            y_max_text = y_horizontal_bar_low + (y_highest_point * text_height_modifier)
            y = y_horizontal_bar_low
            x1 = x_positions[i_cond] + offsets[i_ln1]
            x2 = x_positions[i_cond] + offsets[i_ln2]
            x_mid = float((x1 + x2) / 2)

            ax.text(x_mid, y_max_text, stars, fontsize=text_fontsize, ha='center', va='bottom')
            ax.hlines(y, x1, x2, color='black')

            if h1 < y:
                ax.vlines(x1, y, h1, color='black')
            if h2 < y:
                ax.vlines(x2, y, h2, color='black')

    return ax

def single_paired_plot(
    data,
    x_col,
    y_col,
    id_col,
    condition=None,
    organ_col=None,
    x_pair=None,
    stats_df=None,
    p_val = 'p-value', # allow specifying the column name for the p-value in stats_df
    markers=None,
    ax=None,
    scatter_kwargs=None,
    y_lim=None,
    remove_unpaired=True,
    mouse_id_female=None,
):
    """
    Make a paired plot of one treatment condition between two lymph nodes (pairs).
    """

    np.random.seed(0)

    if scatter_kwargs is None:
        scatter_kwargs = {}
    if ax is None:
        ax = plt.gca()
    if y_lim is None:
        y_lim = math.ceil(data[y_col].max())
    if condition is None:
        raise ValueError(
            'Please provide a condition to plot. If you want to plot all conditions, use the '
            'boxplot_paired function.'
        )
    if x_pair is None:
        raise ValueError('Please provide a pair of lymph nodes to compare.')
    if markers is None:
        markers = ['o', 's']

    if mouse_id_female is not None:
        mouse_id_female = set(mouse_id_female)

    # subset the data
    data_1cond = data[data[x_col] == condition]
    data_1cond_2ln = data_1cond[data_1cond[organ_col].isin(x_pair)].copy()

    if remove_unpaired:
        data_1cond_2ln = data_1cond_2ln.groupby(id_col).filter(lambda x: len(x) == 2)

    # enforce order
    data_sorted = data_1cond_2ln.sort_values(
        organ_col,
        key=lambda x: x.map({x_pair[0]: 0, x_pair[1]: 1})
    ).copy()

    organ_to_x = {x_pair[0]: 0, x_pair[1]: 1}
    organ_to_marker = {x_pair[0]: markers[0], x_pair[1]: markers[1]}

    default_scatter_kwargs = dict(edgecolor='black', zorder=10)
    default_scatter_kwargs.update(scatter_kwargs)

    # plot points
    for organ in x_pair:
        subset = data_sorted[data_sorted[organ_col] == organ].copy()
        if len(subset) == 0:
            continue

        x_coords = np.full(len(subset), organ_to_x[organ], dtype=float)
        x_coords_jittered = x_coords + np.random.uniform(-0.03, 0.03, len(subset))
        y_coords = subset[y_col].to_numpy()

        if mouse_id_female is None:
            ax.scatter(
                x_coords_jittered,
                y_coords,
                marker=organ_to_marker[organ],
                **default_scatter_kwargs
            )
        else:
            base_colour = default_scatter_kwargs.get('color', None)
            if base_colour is None:
                base_colour = default_scatter_kwargs.get('c', None)

            point_colours = [
                'black' if mid in mouse_id_female else base_colour
                for mid in subset[id_col]
            ]
            # check how often the mouse id is in mouse_id_female
            print(f"Mouse IDs in {organ} that are in mouse_id_female: {sum(subset[id_col].isin(mouse_id_female))} out of {len(subset)}")
            # which ones?
            for mid in subset[id_col]:
                if mid in mouse_id_female:
                    print(f"Mouse ID {mid} in {organ} is in mouse_id_female")


            kwargs = default_scatter_kwargs.copy()
            kwargs.pop('color', None)
            kwargs.pop('c', None)

            ax.scatter(
                x_coords_jittered,
                y_coords,
                marker=organ_to_marker[organ],
                c=point_colours,
                **kwargs
            )

    sns.despine()
    ax.set_xlim(-0.5, 1.5)

    # paired lines
    for line in data_sorted[id_col].unique():
        data_line = data_sorted[data_sorted[id_col] == line].copy()
        if len(data_line) < 2:
            continue

        data_line = data_line.sort_values(
            organ_col,
            key=lambda x: x.map(organ_to_x)
        )

        x_coords = [organ_to_x[val] for val in data_line[organ_col].tolist()]
        y_coords = data_line[y_col].tolist()
        ax.plot(x_coords, y_coords, color='black')

    legend = ax.get_legend()
    if legend is not None:
        legend.remove()

    ax.set_xticks([0, 1])
    ax.set_xticklabels(x_pair)

    # significance bars
    if stats_df is not None:
        stats_df_cond = stats_df[(stats_df['Cell type'] == y_col) & (stats_df[x_col] == condition)]
        stats_df_cond = stats_df_cond[
            (stats_df_cond['LN 1'] == x_pair[0]) & (stats_df_cond['LN 2'] == x_pair[1])
        ]

        if len(stats_df_cond) > 0:
            p_val = stats_df_cond[p_val].values[0]
            stars, text_fontsize, text_height_modifier = p_val_formatting(p_val)

            y_highest_point = data[y_col].max()
            h1 = data_sorted[data_sorted[organ_col] == x_pair[0]][y_col].max()
            h2 = data_sorted[data_sorted[organ_col] == x_pair[1]][y_col].max()

            y_horizontal_bar_low = y_highest_point + (y_lim * 0.07)
            h1 = h1 + (y_lim * 0.05)
            h2 = h2 + (y_lim * 0.05)
            y_max_text = y_horizontal_bar_low + (y_highest_point * text_height_modifier)

            ax.text(0.5, y_max_text, stars, fontsize=text_fontsize, ha='center', va='bottom')
            ax.hlines(y_horizontal_bar_low, 0, 1, color='black')

            if h1 < y_horizontal_bar_low:
                ax.vlines(0, y_horizontal_bar_low, h1, color='black')
            if h2 < y_horizontal_bar_low:
                ax.vlines(1, y_horizontal_bar_low, h2, color='black')

    ax.set_ylabel(y_col)
    return ax
