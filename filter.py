import pandas as pd
import numpy as np


def get_colnames():

    path = "./output/"
    unfiltered = pd.read_csv(f"{path}/consensus_unfiltered.csv")

    active_substrings =  [
        "WT",
        "OE",
        "KO",
    ]

    blank_substrings = [
        "blank1"

    ]

    media_substrings = [
        "YES"
    ]

    conditions = {string: [] for string in substrings}
    cols = unfiltered.columns
    
    for col in cols:
        for string in active_substrings:
            if string in col:
                conditions[string].append(col)
    
    return conditions, unfiltered, path

def filter_cols(conditions, unfiltered, min_ratio):
    data_cols = unfiltered.iloc[:, 3:]
    
    ratio_matrix = pd.DataFrame(index=unfiltered.index)

    for condition, columns in conditions.items():
        valid_cols = [c for c in columns if c in data_cols.columns]
        replicate = data_cols[valid_cols]
        valid_features = (replicate > 0).sum(axis=1)

        condition_ratio = valid_features / len(valid_cols)
        ratio_matrix[condition] = condition_ratio

    keep_rows = (ratio_matrix >= min_ratio).any(axis=1)

    filtered_data = unfiltered[keep_rows].copy()

    return filtered_data

def remove_blank():



def main():
    min_ratio = 0.6

    conditions, unfiltered, path = get_colnames()
    filtered_data = filter_cols(conditions,unfiltered,min_ratio)
    output_filepath = f"{path}consensus_filtered.csv"

    filtered_data.to_csv(output_filepath, index=False)
    print(f"saved to {output_filepath}")

main()
