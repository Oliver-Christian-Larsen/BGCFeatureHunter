import copy
import itertools
import math
import random
from collections import Counter

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from scripts import filter_features
from scripts import unique


def load_data(arg_dic, conditions):
    path = f'{arg_dic["paths"]["out"]}/consensus_unfiltered.csv'
    data = pd.read_csv(path, dtype={"consensus_id": str})

    for cols in conditions.values():
        for c in cols:
            if c not in data.columns:
                raise KeyError(f"Column '{c}' not found in {path}")
    return data


def get_sweep_points(config):
    sweep = config["fdr"].get("sweep") or {}
    if not sweep:
        return [("baseline", config)]

    keys = list(sweep.keys())
    points = []
    for values in itertools.product(*[sweep[k] for k in keys]):
        new_config = copy.deepcopy(config)
        label_parts = []
        for key, value in zip(keys, values):
            section, name = key.split(".")
            new_config[section][name] = value
            label_parts.append(f"{name}={value}")
        points.append((", ".join(label_parts), new_config))
    return points


def get_permutable_media(conditions):
    media = {}
    for (group, medium), cols in conditions.items():
        if group in ("WT", "KO"):
            media.setdefault(medium, {"WT": [], "KO": []})[group] = cols

    permutable = {}
    for medium, groups in media.items():
        if len(groups["WT"]) == 0 or len(groups["KO"]) == 0:
            print(f"WARNING: medium {medium} does not have both WT and KO, it is not permuted")
            continue
        permutable[medium] = (groups["WT"] + groups["KO"], tuple(groups["WT"]))
    return permutable


def make_labelings(permutable, fdr_cfg):
    media = list(permutable.keys())
    options = {}
    for m in media:
        all_samples, real_wt = permutable[m]
        options[m] = list(itertools.combinations(all_samples, len(real_wt)))

    total = math.prod(len(options[m]) for m in media)
    n_perm = fdr_cfg["n_permutations"]

    if total <= n_perm:
        candidates = itertools.product(*[options[m] for m in media])
    else:
        rng = random.Random(fdr_cfg["seed"])
        candidates = random_labelings(rng, [options[m] for m in media], total)

    n_samples = sum(len(permutable[m][0]) for m in media)
    labelings = []
    for combo in candidates:
        labeling = dict(zip(media, combo))

        distance = 0
        for m in media:
            distance += len(set(labeling[m]) ^ set(permutable[m][1]))

        if distance < fdr_cfg["min_distance"]:
            continue
        if distance == n_samples and not fdr_cfg["keep_complement"]:
            continue

        labelings.append((labeling, distance))
        if len(labelings) >= n_perm:
            break

    print(f"Permutation space: {total} labelings, using {len(labelings)}")
    return labelings


def random_labelings(rng, options_per_medium, total):
    seen = set()
    while len(seen) < total:
        pick = tuple(rng.randrange(len(o)) for o in options_per_medium)
        if pick in seen:
            continue
        seen.add(pick)
        yield tuple(o[i] for o, i in zip(options_per_medium, pick))


def labeling_to_conditions(conditions, permutable, labeling):
    new_conditions = dict(conditions)
    for medium, wt_samples in labeling.items():
        all_samples = permutable[medium][0]
        new_conditions[("WT", medium)] = [s for s in all_samples if s in wt_samples]
        new_conditions[("KO", medium)] = [s for s in all_samples if s not in wt_samples]
    return new_conditions


def run_labeling(data, config, conditions):
    min_ratio = float(config["filter"]["present_ratio"])
    filtered, ratio_matrix = filter_features.filter_cols(conditions, data, min_ratio)

    pos_cols = []
    neg_cols = []
    for (group, medium), cols in conditions.items():
        if group == "WT":
            pos_cols += cols
        if group == "KO":
            neg_cols += cols

    hits = unique.unique_features(pos_cols, neg_cols, filtered, config, oe_cols=None)
    return hits["consensus_id"].astype(str).tolist()


def summarize(label, real_ids, null_counts, feature_counter):
    null = np.array(null_counts)
    real = len(real_ids)

    null_median = float(np.median(null))
    est_fdr = null_median / real if real > 0 else np.nan
    real_seen_in_null = sum(1 for i in real_ids if feature_counter[i] > 0)

    return {
        "sweep": label,
        "real_n_unique": real,
        "n_permutations": len(null),
        "null_mean": float(np.mean(null)),
        "null_median": null_median,
        "null_q95": float(np.quantile(null, 0.95)),
        "null_max": int(null.max()),
        "est_FDR": est_fdr,
        "real_hits_seen_in_null": real_seen_in_null,
    }


def feature_table(data, label, real_ids, feature_counter, n_perm):
    ids = set(feature_counter.keys()) | set(real_ids)
    meta = data.iloc[:, :4].copy()
    meta["consensus_id"] = meta["consensus_id"].astype(str)
    meta = meta[meta["consensus_id"].isin(ids)].copy()

    meta.insert(0, "sweep", label)
    meta["n_null_hits"] = meta["consensus_id"].map(feature_counter).fillna(0).astype(int)
    meta["null_freq"] = meta["n_null_hits"] / n_perm
    meta["is_real_hit"] = meta["consensus_id"].isin(real_ids)
    return meta.sort_values(["is_real_hit", "null_freq"], ascending=False)


def plot_boxplot(null_by_sweep, summary, arg_dic, path):
    labels = list(null_by_sweep.keys())
    positions = list(range(1, len(labels) + 1))
    real = summary["real_n_unique"].tolist()
    est_fdr = summary["est_FDR"].tolist()

    fig, (ax, ax2) = plt.subplots(2, 1, figsize=(max(6, 1.3 * len(labels) + 3.5), 7.5),
                                  sharex=True, gridspec_kw={"height_ratios": [3, 1.3]})

    box = ax.boxplot(list(null_by_sweep.values()), positions=positions, widths=0.55,
                     patch_artist=True, medianprops={"color": "black"})
    for patch in box["boxes"]:
        patch.set_facecolor("#9ecae1")

    real_points = ax.scatter(positions, real, marker="D", s=60, color="#d62728", zorder=5)
    for x, y in zip(positions, real):
        ax.annotate(str(y), (x, y), textcoords="offset points", xytext=(10, 4), color="#d62728")

    wt = arg_dic["display_names"]["WT"]
    ko = arg_dic["display_names"]["KO"]
    n_perm = summary["n_permutations"].iloc[0]
    ax.set_title(f"{wt} vs {ko}: unique features under within-medium label permutation\n(n={n_perm} permutations)", fontsize=10)
    ax.set_ylabel("# unique features")
    ax.legend([box["boxes"][0], real_points], ["Permuted labels (null)", "Real labels"], loc="upper right")

    fdr_plot = [0 if np.isnan(f) else min(f, 1) for f in est_fdr]
    ax2.bar(positions, fdr_plot, width=0.55, color="#bdbdbd")
    for x, f, h in zip(positions, est_fdr, fdr_plot):
        ax2.text(x, h + 0.02, "n/a" if np.isnan(f) else f"{f:.2f}", ha="center")
    ax2.set_ylim(0, 1.15)
    ax2.set_ylabel("est. FDR\n(median null / real)")
    ax2.set_xticks(positions)
    ax2.set_xticklabels(labels, rotation=30 if len(labels) > 3 else 0, ha="right" if len(labels) > 3 else "center")

    fig.tight_layout()
    fig.savefig(path, dpi=300)
    plt.close(fig)


def main(config, arg_dic, conditions):
    fdr_cfg = config["fdr"]
    out = arg_dic["paths"]["out"]

    data = load_data(arg_dic, conditions)
    permutable = get_permutable_media(conditions)
    if not permutable:
        raise ValueError("No medium contains both WT and KO samples")

    labelings = make_labelings(permutable, fdr_cfg)
    sweep_points = get_sweep_points(config)

    summary_rows = []
    feature_tables = []
    count_rows = []
    null_by_sweep = {}

    for label, sweep_config in sweep_points:
        print(f"Running: {label}")
        real_ids = run_labeling(data, sweep_config, conditions)

        null_counts = []
        feature_counter = Counter()
        for labeling, distance in labelings:
            perm_conditions = labeling_to_conditions(conditions, permutable, labeling)
            ids = run_labeling(data, sweep_config, perm_conditions)
            null_counts.append(len(ids))
            feature_counter.update(set(ids))

        null_by_sweep[label] = null_counts
        summary_rows.append(summarize(label, real_ids, null_counts, feature_counter))
        feature_tables.append(feature_table(data, label, set(real_ids), feature_counter, len(labelings)))

        count_rows.append({"sweep": label, "perm_id": "real", "n_unique": len(real_ids)})
        for i, n in enumerate(null_counts, start=1):
            count_rows.append({"sweep": label, "perm_id": i, "n_unique": n})

    summary = pd.DataFrame(summary_rows)

    permutations = []
    for i, (labeling, distance) in enumerate(labelings, start=1):
        row = {"perm_id": i, "distance_from_real": distance}
        for medium, wt_samples in labeling.items():
            row[f"WT_samples[{medium}]"] = ";".join(wt_samples)
        permutations.append(row)

    summary.to_csv(f"{out}/fdr_summary.csv", index=False)
    pd.DataFrame(count_rows).to_csv(f"{out}/fdr_counts.csv", index=False)
    pd.concat(feature_tables).to_csv(f"{out}/fdr_feature_frequency.csv", index=False)
    pd.DataFrame(permutations).to_csv(f"{out}/fdr_permutations.csv", index=False)
    plot_boxplot(null_by_sweep, summary, arg_dic, f"{out}/fdr_boxplot.png")

    print(summary.to_string(index=False, float_format=lambda x: f"{x:.3g}"))
    print(f"FDR results saved to {out}")


if __name__ == "__main__":
    main()