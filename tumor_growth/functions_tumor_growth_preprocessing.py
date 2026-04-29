#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
# '''=================================================
# @Project: ACT_melanoma_neutrophil (data: act_neutrophil_data_repository)
# @Author : Gemma van der Voort
# @Desc   :
# @Desc updated: Tumor growth preprocessing helpers.
# '''=================================================
import pandas as pd
from datetime import datetime


def calculate_area_from_L_W(df):
    """
    input is a df taken as one of the tables in Maike's growth curve Excel sheets. The function cleans up this table
    and calculates the area from the L and W columns, then returns a df with the area per day. All metadata except
    the ID and cage information is kept (to prevent double index columns).
    :param df: table as extracted from the Excel sheet. usually one of 1-5 tables in a sheet, concat them later.
    :return: area df.
    """
    # drop the first column
    df = df.iloc[:, 1:]
    df.columns = df.iloc[0]
    df = df.iloc[1:]  # remove the first row after setting it as column names
    # starting from where row 2 is not nan for the first time,
    L_W_row = df.iloc[2, :]
    first_non_nan_position = L_W_row.notna().values.argmax()
    length_indexes = list(range(first_non_nan_position, len(df.columns), 2))
    width_indexes = list(range(first_non_nan_position + 1, len(df.columns), 2))
    # get experiment day, which are all the non-nan entries in the 2nd row of df_87_1_part_1_clean
    exp_day = df.iloc[1].dropna()
    exp_day = exp_day.astype(str).tolist()  # get as list of strings to use as column names
    # assert exp_day and length_indexes and width_indexes are the same length
    assert len(exp_day) == len(length_indexes) == len(width_indexes)
    shared_columns = ['ID', 'Strain', 'DOB', 'Age (wks)', 'Ear notch', 'Sex', 'Cell line injected']
    area_df = pd.DataFrame(columns=shared_columns + exp_day)
    area_df[shared_columns] = df[shared_columns]
    for exp_day, length_idx, width_idx in zip(exp_day, length_indexes, width_indexes):
        length = df.iloc[:, length_idx]
        width = df.iloc[:, width_idx]
        # check that the 3rd row contains L for length and W for width
        assert length.iloc[2] == 'L'
        assert width.iloc[2] == 'W'
        area = length.iloc[3:, ] * width.iloc[3:, ]
        area_df[exp_day] = area
    area_df = area_df.iloc[3:]  # drop the first 3 rows, as they now only contain nans.
    return area_df


def extract_coloured_events(ws, events_list, colour_list, colour_type_list):
    """
    extracts metadata events encoded in Excel using colors. Returns metadata events. Additionally, a list of found
    colors is returned that is intended to help troubleshooting when events are missing due to slight color variations.
    :param ws:
    :param events_list:
    :param colour_list:
    :param colour_type_list:
    :return:
    """
    assert len(events_list) == len(colour_list) == len(colour_type_list), \
        "events_list, colour_list, and colour_type_list must be the same length"

    def is_date(value):
        return isinstance(value, (datetime,))

    def get_color(cell, color_type):
        if color_type == 'font':
            color = cell.font.color
        elif color_type == 'fill':
            color = cell.fill.start_color
        else:
            raise ValueError(f"Invalid color_type: {color_type}")
        if color and color.type == 'rgb':
            return color.rgb
        return None

    def find_date_above(ws, cell):
        col_letter = cell.column_letter
        row = cell.row
        for r in range(row - 1, 0, -1):
            val = ws[f"{col_letter}{r}"].value
            if is_date(val):
                return val, ws[f"{col_letter}{r}"]
        raise ValueError(f"No date found above cell {cell.coordinate}")

    def find_experiment_day(ws, date_cell):
        col = date_cell.column
        row = date_cell.row
        below_cell = ws.cell(row=row + 1, column=col)
        if below_cell.value is None:
            raise ValueError(f"No experiment day below date cell {date_cell.coordinate}")
        return below_cell.value, below_cell

    def find_mouse_id(ws, cell):
        row = cell.row
        col = cell.column
        for c in range(col - 1, 0, -1):
            value = ws.cell(row=row, column=c).value
            if isinstance(value, str) and (value.upper().startswith('BAL') or value.upper().startswith('HOE')):
                return value
        raise ValueError(f"No mouse ID found to the left of cell {cell.coordinate}")

    # Normalize colour_list to be list-of-lists
    normalized_colour_list = [
        cl if isinstance(cl, list) else [cl]
        for cl in colour_list
    ]

    records = []
    found_colours = set()


    processed_rows = set()

    for row in ws.iter_rows():
        row_idx = row[0].row
        if row_idx in processed_rows:
            continue

        for cell in row:
            # Track all colors found (font and fill)
            font_rgb = get_color(cell, 'font')
            fill_rgb = get_color(cell, 'fill')
            if font_rgb: found_colours.add(f'font:{font_rgb}')
            if fill_rgb: found_colours.add(f'fill:{fill_rgb}')

            # Check each event
            for event, allowed_colors, color_type in zip(events_list, normalized_colour_list, colour_type_list):
                cell_color = get_color(cell, color_type)
                if cell_color in allowed_colors:
                    try:
                        date_val, date_cell = find_date_above(ws, cell)
                        experiment_day, day_cell = find_experiment_day(ws, date_cell)
                        mouse_id = find_mouse_id(ws, cell)

                        records.append({
                            'date': date_val,
                            'experiment date': experiment_day,
                            'experiment date cell': day_cell.coordinate,
                            'mouse_ID': mouse_id,
                            'value': cell.value,
                            'value cell': cell.coordinate,
                            'event': event,
                            'sheet': ws.title
                        })

                        processed_rows.add(row_idx)
                        break  # stop processing more events for this row

                    except ValueError as e:
                        print(f"[{ws.title}] Error: {e}")
                        processed_rows.add(row_idx)
    return pd.DataFrame(records), found_colours
