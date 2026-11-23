import os
import glob
from collections import defaultdict



def get_mzML_paths(arg_dic):
    mzml_folder = arg_dic["paths"]["mzml"]
    search_pattern = os.path.join(mzml_folder, "*.mzML")

    mzml_files = glob.glob(search_pattern)
    mzml_files.sort()

    return mzml_files

def find_matches(filename, substrings):

    matches = []
    filename_lower = filename.lower()

    for label, substring in substrings.items():
        if substring is None:
            continue

        if substring.lower() in filename_lower:
            matches.append(label)

    return matches

def classify_one_file(file_path, group_substrings, media_substrings):
    filename = os.path.basename(file_path)

    group_matches = find_matches(filename, group_substrings)

    if len(group_matches) > 1:
        raise ValueError(
            f"Ambiguous file: {filename} matches several groups: {group_matches}"
        )

    if len(group_matches) == 1:
        group = group_matches[0]
    else:
        group = None

    media_dict = {}
    for media_name in media_substrings:
        media_dict[media_name] = media_name

    media_matches = find_matches(filename, media_dict)

    if len(media_matches) > 1:
        raise ValueError(
            f"Ambiguous file: {filename} matches several media: {media_matches}"
        )

    if len(media_matches) == 1:
        media = media_matches[0]
    else:
        media = None
    return {
        "filepath": file_path,
        "filename": filename,
        "group": group,
        "media": media,
    }

def print_summary(classified_files):
    print("\n--- File classification ---")
    print(f"A total of {len(classified_files)} has been identified")
    for info in classified_files:
        print(f"{info['filename']:50} group={info['group']}  media={info['media']}")

    unclassified = []
    for info in classified_files:
        if info["group"] is None:
            unclassified.append(info["filename"])

    if unclassified:
        print("\nWARNING: these files matched no group. They will be preprocessed, but not used in subsequent analysis:")
        for name in unclassified:
            print(f"{name}")


def get_conditions(classified_files, mzml_files):
    mzml_files_no_ext = []
    for f in mzml_files:
        fname = os.path.basename(f)
        mzml_files_no_ext.append(os.path.splitext(fname)[0])

    conditions = defaultdict(list)

    for info in classified_files:
        group = info["group"]
        if group is None or group == "blank":
            continue

        col = os.path.splitext(info["filename"])[0]
        if col not in mzml_files_no_ext:
            print(f"WARNING: no column '{col}' in consensus table, skipping")
            continue

        conditions[(group, info["media"])].append(col)

    return dict(conditions)

def check_groups(classified_files, arg_dic):
    counts = {}
    for info in classified_files:
        counts[info["group"]] = counts.get(info["group"], 0) + 1

    required = {"WT": arg_dic["strain"]["WT"], "KO": arg_dic["strain"]["KO"]}
    if arg_dic["strain"]["OE"] is not None:
        required["OE"] = arg_dic["strain"]["OE"]

    problems = [
        f"{label}: no filenames contain '{substring}' (set with --{label.lower()})"
        for label, substring in required.items()
        if counts.get(label, 0) == 0
    ]
    if problems:
        raise SystemExit(
            "\nERROR: some groups matched no .mzML files:\n"
            + "\n".join(problems)
            + "\nCheck the substrings against the filenames listed above "
            "(matching is case-insensitive)."
        )



def main(arg_dic):
    mzml_files = get_mzML_paths(arg_dic)
    if len(mzml_files) == 0:
        raise FileNotFoundError(f"No .mzML files found in {arg_dic['paths']['mzml']}")

    group_substrings = {
        "WT": arg_dic["strain"]["WT"],
        "KO": arg_dic["strain"]["KO"],
        "OE": arg_dic["strain"]["OE"],
        "blank": arg_dic["blank"],
    }
    media_substrings = arg_dic["media"]

    classified_files = []
    for file_path in mzml_files:
        file_info = classify_one_file(file_path, group_substrings, media_substrings)
        classified_files.append(file_info)

    conditions = get_conditions(classified_files, mzml_files)

    print_summary(classified_files)
    check_groups(classified_files, arg_dic)
    return classified_files,mzml_files,conditions



if __name__ == "__main__":
    main()

