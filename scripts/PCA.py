import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler


def load_data(arg_dic):
    filtered_df = pd.read_csv(
        f'{arg_dic["paths"]["out"]}/consensus_filtered.csv'
    )

    filtered_df.drop(
        ["rt", "mz", "quality", "consensus_id"],
        axis=1,
        inplace=True
    )

    return filtered_df


def pca(filtered_df, conditions, arg_dic):

    samples = filtered_df.columns

    X = StandardScaler().fit_transform(filtered_df.T)

    model = PCA()
    scores = model.fit_transform(X)

    explained = model.explained_variance_ratio_
    cumulative = np.cumsum(explained)

    sample_conditions = {}

    for (strain, medium), condition_samples in conditions.items():
        for sample in condition_samples:
            sample_conditions[sample] = (strain, medium)

    strains = sorted(set(strain for strain, medium in conditions))
    media = sorted(set(medium for strain, medium in conditions))

    colors = plt.cm.tab10(np.linspace(0, 1, len(strains)))
    strain_colors = dict(zip(strains, colors))

    markers = ["o", "^", "s", "D", "v", "P", "X", "<", ">"]

    if len(media) > len(markers):
        print(f"WARNING: The length of media surpasses the number of available shapes \n"
                "The PCA-plot will have duplicate shape entries for different media types \n"
                f"Only {len(media)} shapes available")

    media_markers = {
        medium: markers[i % len(markers)]
        for i, medium in enumerate(media)
    }

    fig, ax = plt.subplots(figsize=(6, 4))

    ax.plot(
        np.arange(1, len(cumulative) + 1),
        cumulative * 100,
        color="black"
    )

    ax.set_xlabel("Number of principal components")
    ax.set_ylabel("Cumulative explained variance (%)")

    fig.tight_layout()
    fig.savefig(
        f'{arg_dic["paths"]["out"]}/PCA_cumulative.png', dpi=300,
        bbox_inches="tight"
    )
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6, 5))

    for i, sample in enumerate(samples):

        if sample not in sample_conditions:
            continue

        strain, medium = sample_conditions[sample]

        ax.scatter(
            scores[i, 0],
            scores[i, 1],
            color=strain_colors[strain],
            marker=media_markers[medium],
            s=70,
            edgecolor="black",
            linewidth=0.5
        )

    ax.set_xlabel(f"PC1 ({explained[0] * 100:.1f}%)")
    ax.set_ylabel(f"PC2 ({explained[1] * 100:.1f}%)")

    strain_handles = []

    for strain in strains:
        label = arg_dic["display_names"].get(strain, strain)

        strain_handles.append(
            plt.Line2D(
                [0], [0],
                marker="o",
                linestyle="",
                markerfacecolor=strain_colors[strain],
                markeredgecolor="black",
                markersize=8,
                label=label
            )
        )

    media_handles = []

    for medium in media:
        media_handles.append(
            plt.Line2D(
                [0], [0],
                marker=media_markers[medium],
                linestyle="",
                color="black",
                markerfacecolor="white",
                markersize=8,
                label=medium
            )
        )

    strain_legend = ax.legend(
        handles=strain_handles,
        title="Strain",
        frameon=False,
        loc="upper left",
        bbox_to_anchor=(1.02, 1)
    )

    ax.add_artist(strain_legend)

    ax.legend(
        handles=media_handles,
        title="Medium",
        frameon=False,
        loc="lower left",
        bbox_to_anchor=(1.02, 0)
    )

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    fig.tight_layout()


    fig.savefig(
        f'{arg_dic["paths"]["out"]}/PCA_plot.png',
        dpi=300,
        bbox_inches="tight"
    )

    plt.close(fig)


def main(arg_dic, conditions):
    filtered_df = load_data(arg_dic)
    pca(filtered_df, conditions, arg_dic)

if __name__ == "__main__":
    main()