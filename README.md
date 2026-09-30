A simple PyOpenMS-based pipeline for detection and visualization of candidate features, originating from a genetically engineered biosynthetic gene cluster.

The input to the pipeline is calibrated centroid .mzML files. Brifly, the workflow is as follows:

1) Enter centroided and calibrated .mzML files into a folder, with replicates of the same condition having distinct, identifiable names (default is './cal_centroid_mzML').
2) run preproccessing.py, which outputs a featureXML file for each .mzML file, and an unfiltered .csv file
3) run filter.py on the unfiltered .csv file
4) run unique.py on the filtered .csv file
5) run waterfall.py on the unique features .csv file. This yields a pdf file, with potential hits.
