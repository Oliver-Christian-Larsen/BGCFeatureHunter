import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA


def load_data(arg_dic):
    filtered_df = pd.read_csv(f"{arg_dic["paths"]["out"]}/consensus_filtered.csv")

    filtered_df.drop(["rt","mz","quality","consensus_id"],axis=1, inplace = True)

    return filtered_df
    





def pca(pe_df):
    pe_pca_df = pe_df.drop(["id"], axis=1)
    pca = PCA(n_components=20)
    pe_PCA = pca.fit_transform(pe_pca_df)
    prpca = pca.explained_variance_ratio_
    cumulative = np.cumsum(pca.explained_variance_ratio_)

    plt.plot(cumulative)
    plt.title("Kummulativ forklaring af embeddings i PC components")
    plt.savefig("cumulative.png")
    plt.close()

    plt.plot(prpca)
    plt.savefig("prpca.png")
    plt.close()

    plt.scatter(pe_PCA[:,0], pe_PCA[:,1], s=1, alpha = 0.1)
    plt.savefig("PCA_plot.png")
    plt.close()





def main(arg_dic,classified_files):
    for f in classified_files:
        print(f)
    f_df = load_data(arg_dic)
    print(f_df.head())