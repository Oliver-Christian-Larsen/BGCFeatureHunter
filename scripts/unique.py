import pandas as pd
import numpy as np

def unique_features():

    data = pd.read_csv("output/consensus_filtered.csv")

    data_cols = data.iloc[:, 3:]
    positive_groups = [
    "OE_Y",
    "OE_P",

    "WT_Y",
    "WT_P",
]
    negative_groups = [
    "KO_Y",
    "KO_P"
]

    pos_cols = [col for col in data_cols if any(sub in col for sub in positive_groups)]
    neg_cols = [col for col in data_cols if any(sub in col for sub in negative_groups)]

    if not pos_cols or not neg_cols:
        raise ValueError("Group substrings did not match any columns in the dataframe.")

    overlap = set(pos_cols).intersection(set(neg_cols))
    if overlap:
        raise ValueError(f"Ambiguous definitions: Columns {overlap} match both groups.")

    is_in_positive = (data[pos_cols] > 0.0).any(axis=1)

    is_absent_in_negative = (data[neg_cols] == 0.0).all(axis=1)

    unique_features = is_in_positive & is_absent_in_negative

    metadata_cols = data.columns[:3].tolist()
    cols_to_keep = metadata_cols + pos_cols + neg_cols

    filtered_df = data.loc[unique_features, cols_to_keep].copy()

    summary_df = filtered_df[metadata_cols].copy()

    summary_df['mean_intensity_positive'] = filtered_df[pos_cols].mean(axis=1)

    summary_df = summary_df.sort_values(by='mean_intensity_positive', ascending=False)


    return summary_df

def known_features(unique_data):
    known = [
    ]

    for i in known:
        mass_tolerance = i * 5 * 1e-6
        lower_mass = i - mass_tolerance
        upper_mass = i + mass_tolerance

        if unique_data['mz'].between(lower_mass, upper_mass).any():
            print(f"Target mass {i} is found!")
        else:
            print(f"Target mass {i} is not found!")

def main():
    unique_data = unique_features()

    unique_data.to_csv('./output/unique_features.csv', index=False)

    #known_features(unique_data)
    
if __name__ == "__main__":
    main()
