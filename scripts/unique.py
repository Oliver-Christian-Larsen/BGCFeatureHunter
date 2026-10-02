import pandas as pd
import numpy as np


def get_cols(arg_dic,classified_files):
    pos_cols = []
    neg_cols = []
    if arg_dic.get("OE") is not None:
        oe_cols = []
    for c in classified_files:
        filename = c["filename"].split(".")[0]
        print(filename)
        if c["group"] == "WT":
            pos_cols.append(filename)
        if c["group"] == "KO":
            neg_cols.append(filename)
        if c["group"] == "OE" and arg_dic.get("OE") is not None:
            oe_cols.append(filename)
    print(pos_cols,neg_cols)

    if arg_dic.get("OE") is not None:
        print("wrong path")
        return pos_cols, neg_cols, oe_cols
    else:
        print("right path")
        return pos_cols, neg_cols

def unique_features(pos_cols, neg_cols, data,config,oe_cols):

    if config["unique"]["allow_KO_detection"] == True:
        wt = data[pos_cols].mean(axis=1)
        ko = data[neg_cols].mean(axis=1)
        log2_wt_over_ko = np.log2(wt + 1) - np.log2(ko + 1)
        is_LFC = log2_wt_over_ko > config["unique"]["WT_LFC_to_KO"]

        unique_features = is_LFC

    else:
        is_absent_in_negative = (data[neg_cols] == 0.0).all(axis=1)
        is_in_positive = (data[pos_cols] > 0.0).any(axis=1)
        unique_features = is_in_positive & is_absent_in_negative

    metadata_cols = data.columns[:4].tolist()
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

def main(config,arg_dic,classified_files):
    if arg_dic.get("OE") is not None:
        pos_cols, neg_cols, oe_cols = get_cols(arg_dic,classified_files)
    else:
        pos_cols, neg_cols = get_cols(arg_dic,classified_files)

    
    data = pd.read_csv(f"{arg_dic["paths"]["out"]}/consensus_filtered.csv")
    if arg_dic.get("OE") is not None:
        unique_data = unique_features(pos_cols,neg_cols,data,config,oe_cols)
    else:
        unique_data = unique_features(pos_cols,neg_cols,data,config,oe_cols=None)

    unique_data.to_csv(f'{arg_dic["paths"]["out"]}/unique_features.csv', index=False)

    #known_features(unique_data)

if __name__ == "__main__":
    main()
