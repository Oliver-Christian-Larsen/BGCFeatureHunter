import pandas as pd
from collections import defaultdict
from pathlib import Path


def get_conditions(classified_files, columns):
    available = set(columns)
    conditions = defaultdict(list)

    for info in classified_files:
        group = info["group"]
        if group is None or group == "blank":
            continue

        col = info["filename"].replace(".mzML", "")
        if col not in available:
            print(f"WARNING: no column '{col}' in consensus table, skipping")
            continue

        conditions[(group, info["media"])].append(col)

    return dict(conditions)

def filter_cols(conditions, unfiltered, min_ratio):
    data_cols = unfiltered.iloc[:, 4:]
    ratio_matrix = pd.DataFrame(index=unfiltered.index)

    for (group, media), columns in conditions.items():
        columns = list(dict.fromkeys(columns))
        present = (data_cols[columns] > 0).sum(axis=1)
        ratio_matrix[f"{group}_{media}"] = present / len(columns)

    keep_rows = (ratio_matrix >= min_ratio).any(axis=1)
    return unfiltered[keep_rows].copy(), ratio_matrix


def main(config, arg_dic, classified_files):
    path = arg_dic["paths"]["out"]
    unfiltered = pd.read_csv(f"{path}/consensus_unfiltered.csv")
    min_ratio = float(config["filter"]["present_ratio"])

    conditions = get_conditions(classified_files, unfiltered.columns[4:])

    filtered, ratio_matrix = filter_cols(conditions, unfiltered, min_ratio)

    filtered.to_csv(f"{path}/consensus_filtered.csv", index=False)
    print(f"{len(filtered)}/{len(unfiltered)} features survived the filter")