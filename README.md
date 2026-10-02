A simple PyOpenMS-based pipeline for detection and visualization of candidate features, originating from a genetically engineered biosynthetic gene cluster. 

The input to the pipeline is calibrated centroid .mzML files. Brifly, the workflow is as follows:
Preprocessing with PyOpenMS, using a set of standard configurations, found in config.yaml. After Preprocessing, the feature table is filtered, and scanned for features marked as unique (default is present in WT, absent in KO). For all detected features, a waterfall plot is generated, in a combined .pdf file.

Ensure your filenames of the input .mzml files are distinguishable, and run. You can add flags, corresponding to your specific substrings (run python3 run_pipeline.py -h).

To run, the program, simply install dependencies, and run python3 run_pipeline.py [enter specific settings]


Feel free to reach out
