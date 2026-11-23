import pandas as pd
import numpy as np
import matplotlib.pyplot as plt



def plot(unique_df,filtered_df,arg_dic):
    unique_data = unique_df["mz"].to_list()
    unique_data_int = unique_df["mean_intensity_positive"].to_numpy()
    unique_data_int_log2 = np.log2(unique_data_int + 1)
    filtered_data = filtered_df["mz"].to_list()

    dec_unique, num_unique = np.modf(unique_data)
     
    dec_filt, num_filt = np.modf(filtered_data)
    fig, ax = plt.subplots(figsize=(8, 5))

    
    sc = ax.scatter(
        num_unique,
        dec_unique,
        c=unique_data_int_log2,
        alpha=0.6,
        s=7,
        cmap="viridis",
    )

    ax.scatter(
        num_filt,
        dec_filt,
        c="grey",
        alpha=0.05,
        s=5,
    )
    ax.set_title("Mass Defect plot of unique and filtered features")
    ax.set_ylabel("Decimal part")
    ax.set_xlabel("Integer part")
    cbar = plt.colorbar(sc)
    cbar.set_label("Log2(intensity + 1) of unique features")
    plt.savefig(f"{arg_dic["paths"]["out"]}/mass_defect.png",dpi=300)


def main(arg_dic):
    unique_df = pd.read_csv(f"{arg_dic["paths"]["out"]}/unique_features.csv")
    filtered_df = pd.read_csv(f"{arg_dic["paths"]["out"]}/consensus_filtered.csv")
    plot(unique_df,filtered_df,arg_dic)

if __name__ == "__main__":
    main()