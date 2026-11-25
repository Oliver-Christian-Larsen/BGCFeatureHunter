import argparse
import os
from importlib.resources import files
import yaml

def load_config():
    config = files("bgcfeaturehunter").joinpath("config.yaml").read_text()

    return yaml.safe_load(config)



def main():

    config = load_config()

    parser = argparse.ArgumentParser(description='Allow user to specify specific conditions',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--media', type=str, nargs='+', help='One or more names of the media used. Use substrings, that are present in the mzML filename, eg. PDA.')
    
    parser.add_argument('--blank', type=str, help='The substring used to identify the samples from the blank if used. It will be subtracted, so both injection blanks and media backgrounds are accepted (currently does nothing)')
    
    parser.add_argument('--wt', type=str, default="WT", help='The substring used to identify the samples from the WT, or background')
    parser.add_argument('--ko', type=str, default="KO", help='The substring used to identify the samples from the deletion mutant')
    parser.add_argument('--oe', type=str, help='The substring used to identify the samples from an overexpression mutant')
    
    parser.add_argument('--wt-display', type=str, default="WT", help='The showname in the plot for WT')
    parser.add_argument('--ko-display', type=str, default="KO", help='The showname in the plot for KO')
    parser.add_argument('--oe-display', type=str, default="OE", help='The showname in the plot for OE')
    
    parser.add_argument('--mzml', type=str, default= "./data", help='The filepath to the calibrated centroided .mzML files')
    parser.add_argument('--out', type=str,  default= "./output", help='The output filepath')
    
    parser.add_argument('--no_pdf', action='store_true', help='Flag to disable the production of the waterfall pdf')
    parser.add_argument('--adducts', action='store_true', help='Detect adducts. Only runs on the features in the unique features list')
    parser.add_argument('--fdr', action='store_true', help='Run the FDR permutation test after the pipeline (settings in config.yaml)')
    
    parser.add_argument('--elemental', action='store_true', help='If you want molecular formula calculation. By default, only runs if a feature has 1 > detected ions')
    parser.add_argument('--elemental_single', action='store_true', help='If you want molecular formula calculation for single ions. Only runs if --elemental is used')
    
    parser.add_argument('--elemental_feature_file', type=str, default=None, help='The feature file you want to calculate molecular formulas for. Defaults to <out>/unique_features.csv')
    
    args = parser.parse_args()
    
    
    if args.elemental_feature_file is not None and not args.elemental:
        parser.error("--elemental_feature_file requires --elemental")
    
    if args.elemental_feature_file is None:
        args.elemental_feature_file = os.path.join(args.out, "unique_features.csv")
    
    arg_dic = {
        "strain" : {
            "KO" : args.ko, 
            "WT" : args.wt,
            "OE" : args.oe,
        },
        "display_names": {
            "KO": args.ko_display,
            "WT": args.wt_display,
            "OE": args.oe_display,
        },
        "media" : args.media or [],
        "blank": args.blank,
        "paths" : {
            "mzml" : args.mzml,
            "out" : args.out,
            "elemental_featurefile" : args.elemental_feature_file,
        },
        "elemental_single" : args.elemental_single
    }
    
    
    print("Loading relevant libraies")
    from bgcfeaturehunter.scripts import classify_files
    from bgcfeaturehunter.scripts import preprocessing as pp
    from bgcfeaturehunter.scripts import filter_features
    from bgcfeaturehunter.scripts import unique
    from bgcfeaturehunter.scripts import waterfall
    from bgcfeaturehunter.scripts import PCA
    from bgcfeaturehunter.scripts import mass_deficit
    from bgcfeaturehunter.scripts import FDR
    from bgcfeaturehunter.scripts import adduct_detection

    
    print("Loaded libraries, running pipeline!")
    
    
    classified_files,mzml_files,conditions = classify_files.main(arg_dic)
    
    ref_file = pp.main(config,arg_dic,mzml_files)

    filter_features.main(config,arg_dic,classified_files,conditions)
    
    PCA.main(arg_dic,conditions)
    
    unique.main(config,arg_dic,classified_files)
    
    if args.fdr:
        FDR.main(config, arg_dic, conditions)
    
    mass_deficit.main(arg_dic)
    if args.adducts:
        adduct_detection.main(config,arg_dic)
    
    if args.elemental:
        from bgcfeaturehunter.elemental import calculate_MF
        from bgcfeaturehunter.elemental import get_relative_intensities
        ratio_data = get_relative_intensities.main(arg_dic)
        calculate_MF.main(ratio_data,arg_dic)
    
    if args.no_pdf == True:
        print("Job done! No visualization rendered")
    else:
        unique_csv = os.path.join(args.out, "unique_features.csv")
        with open(unique_csv) as fh:
            n_unique = sum(1 for _ in fh) - 1
        if n_unique == 0:
            print("No unique features found, so no waterfall PDF was made. "
                  "Consider lowering filter.present_ratio or setting "
                  "unique.allow_KO_detection: true in the config."
                  "Alternatively, lower the noise threshold.")
        else:
            pipe = waterfall.WaterfallPipeline(arg_dic, classified_files)
            pipe.run(arg_dic, config)
    
if __name__ == "__main__":
    main()