# BGCFeatureHunter

A simple PyOpenMS-based pipeline for detection and visualization of candidate features originating from a genetically engineered biosynthetic gene cluster (gene deletions, overexpression, heterologous expression).

The input is calibrated, centroided .mzML files. msconvert form proteowizard is a usefull tool for converting your files.
Briefly, the workflow is as following:

1. Preprocessing with PyOpenMS, using the standard settings in config.yaml
2. Filtering of the feature table
3. Detection of unique features (default: present in WT, absent in KO)
4. A waterfall plot for each candidate feature, combined in one PDF

## Installation

Requires Python 3.12 or 3.13.


pip install git+https://github.com/Oliver-Christian-Larsen/BGCFeatureHunter
bgcfeaturehunter --help


To edit the code or settings, install from a clone instead:


git clone https://github.com/Oliver-Christian-Larsen/BGCFeatureHunter
cd BGCFeatureHunter
pip install -e .


The settings file is then 'src/bgcfeaturehunter/config.yaml'.

## Usage

Put your '.mzML' files in one folder. The group and medium options are case-insensitive **substrings of the filenames**, so make sure they identify your files unambiguously.

bgcfeaturehunter --mzml ./data --out ./output --wt WT --ko Delta

Useful extras: '--media PDA YES', '--adducts', '--fdr', '--elemental', '--no_pdf'. Run 'bgcfeaturehunter --help' for all options.

Results are written to the output folder, including 'unique_features.csv' and 'waterfall.pdf'.

Feel free to reach out.