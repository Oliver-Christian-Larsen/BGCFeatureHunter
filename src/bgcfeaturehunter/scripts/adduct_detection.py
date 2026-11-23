import pandas as pd
import numpy as np


def detect_adducts(config,unique_data,arg_dic):    
    mz = unique_data["mz"].to_numpy()
    rt = unique_data["rt"].to_numpy()
    cid = unique_data["consensus_id"].to_numpy()

    adducts = config["adducts"]
    rt_diff = config["adduct_settings"]["adduct_rt_diff"]
    ppm_diff = config["adduct_settings"]["adduct_ppm_diff"]

    mz_diff = mz[None, :] - mz[:, None]
    rt_ok = np.abs(rt[None, :] - rt[:, None]) <= rt_diff
    mz_tol = mz[:, None] * ppm_diff * 1e-6


    adduct = np.full(len(mz), None, dtype=object)
    parent_id = np.full(len(mz), None, dtype=object)


    for name, delta in config["adducts"].items():
        match = (np.abs(mz_diff - delta) <= mz_tol) & rt_ok
        i_idx, j_idx = np.nonzero(match)
        adduct[j_idx] = name
        parent_id[j_idx] = cid[i_idx]

    adducts_df = unique_data.copy()
    adducts_df["adduct"] = adduct
    adducts_df["adduct_parent_id"] = parent_id
    print(adducts_df)

    
    return adducts_df


def main(config,arg_dic):
    unique_data = pd.read_csv(f'{arg_dic["paths"]["out"]}/unique_features.csv')
    adducts_df = detect_adducts(config,unique_data,arg_dic)
    adducts_df.to_csv(f"{arg_dic["paths"]["out"]}/unique_features.csv",index=False)
