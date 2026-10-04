import pandas as pd
from pathlib import Path
from collections import defaultdict


def filter_cols(conditions, unfiltered, min_ratio):
    data_cols = unfiltered.iloc[:, 4:]
    ratio_matrix = pd.DataFrame(index=unfiltered.index)

    for (group, media), columns in conditions.items():
        columns = list(dict.fromkeys(columns))
        present = (data_cols[columns] > 0).sum(axis=1)
        ratio_matrix[f"{group}_{media}"] = present / len(columns)

    keep_rows = (ratio_matrix >= min_ratio).any(axis=1)
    return unfiltered[keep_rows].copy(), ratio_matrix


def main(config, arg_dic, classified_files,conditions):
    path = arg_dic["paths"]["out"]
    unfiltered = pd.read_csv(f"{path}/consensus_unfiltered.csv")
    min_ratio = float(config["filter"]["present_ratio"])


    filtered, ratio_matrix = filter_cols(conditions, unfiltered, min_ratio)

    filtered.to_csv(f"{path}/consensus_filtered.csv", index=False)
    print(f"{len(filtered)}/{len(unfiltered)} features survived the filter")

    return conditions