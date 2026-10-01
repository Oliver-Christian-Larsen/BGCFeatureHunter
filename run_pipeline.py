from scripts import classify_files
from scripts import preproccessing as pp
from scripts import filter
from scripts import unique
from scripts import waterfall

import yaml
import argparse

parser = argparse.ArgumentParser(description='Allow user to specify specific conditions',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--media', type=str, nargs='+', help='One or more names of the media used. Use substrings, that are present in the mzML filename, eg. PDA.')

parser.add_argument('--blank', type=str, help='The substring used to identify the samples from the blank if used. It will be subtracted, so both injection blanks and media backgrounds are accepted')

parser.add_argument('--wt', type=str, default="WT", help='The substring used to identify the samples from the WT, or background')
parser.add_argument('--ko', type=str, default="KO", help='The substring used to identify the samples from the deletion mutant')
parser.add_argument('--oe', type=str, help='The substring used to identify the samples from an overexpression mutant')

parser.add_argument('--mzml', type=str, default= "./data", help='The filepath to the calibrated centroided .mzML files')
parser.add_argument('--out', type=str,  default= "./output", help='The output filepath')

args = parser.parse_args()

arg_dic = {
    "strain" : {
        "KO" : args.ko,
        "WT" : args.wt,
        "OE" : args.oe,
    },
    "media" : args.media or [],
    "blank": args.blank,
    "paths" : {
        "mzml" : args.mzml,
        "out" : args.out,
    }
}

classified_files,mzml_files = classify_files.main(arg_dic)

with open('config.yaml', 'r') as file:
    config = yaml.safe_load(file)

pp.main(config,arg_dic,mzml_files)
filter.main(config)

unique.main(config)

pipe = waterfall.WaterfallPipeline(arg_dic)
pipe.run()