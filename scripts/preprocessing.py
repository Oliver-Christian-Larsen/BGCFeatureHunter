import pyopenms as oms
import glob
import os
import pandas as pd
import gc

def detect_mass_traces(exp,config):
    mass_traces = []
    mtd = oms.MassTraceDetection()
    mtd.setLogType(oms.LogType.NONE)
    mtd_params = mtd.getDefaults()
    mtd_params.setValue("mass_error_ppm", float(config["preprocessing"]["mass_trace_ppm"]))  
    mtd_params.setValue("noise_threshold_int", float(config["preprocessing"]["noise_threshold"])) 
    mtd.setParameters(mtd_params)
    mtd.run(exp, mass_traces, 0)
    return mass_traces

def detect_elution_peaks(mass_traces):
    mass_traces_split, mass_traces_final = [], []
    epd = oms.ElutionPeakDetection()
    epd_params = epd.getDefaults()
    epd.setLogType(oms.LogType.NONE)
    epd.getParameters().getValue("width_filtering") == "auto"
    epd.setParameters(epd_params)
    epd.detectPeaks(mass_traces, mass_traces_split)

    return mass_traces_final


def detect_features_from_traces(mass_traces_final,config):
    fmap = oms.FeatureMap()
    feat_chrom = []
    ffm = oms.FeatureFindingMetabo()
    ffm.setLogType(oms.LogType.NONE)
    ffm_params = ffm.getDefaults()
    ffm_params.setValue("isotope_filtering_model", "none")
    ffm_params.setValue("remove_single_traces", str(config["preprocessing"]["remove_single_traces"])) 
    ffm_params.setValue("report_convex_hulls", "true")
    ffm.setParameters(ffm_params)
    ffm.run(mass_traces_final, fmap, feat_chrom)
    fmap.setUniqueIds()

    return fmap


def save_features(feature_map, filename, path, output_path):
    base_name = os.path.splitext(filename)[0]
    output_dir = os.path.join(output_path, './featureXML')
    os.makedirs(output_dir, exist_ok=True)

    output_file = os.path.join(output_dir, f"{base_name}.featureXML")
    oms.FeatureXMLFile().store(output_file, feature_map)
    return output_file

def align_feature_maps(feature_maps, filenames):
    print("Aligning feature maps")
    ref_index = max(range(len(feature_maps)), key=lambda i: feature_maps[i].size())
    ref_map = feature_maps[ref_index]
    ref_filename = filenames[ref_index]
    print(f"Using '{ref_filename}' as reference (features: {ref_map.size()})")

    aligner = oms.MapAlignmentAlgorithmPoseClustering()
    aligner.setLogType(oms.LogType.NONE)
    aligner.setReference(ref_map)
    params = aligner.getParameters()
    aligner.setParameters(params)

    for i, fmap in enumerate(feature_maps):
        if i == ref_index: continue
        trafo = oms.TransformationDescription()
        aligner.align(fmap, trafo)
        transformer = oms.MapAlignmentTransformer()
        transformer.transformRetentionTimes(fmap, trafo, True) 
    return feature_maps, ref_filename, ref_index

    
def link_features(feature_maps, filenames,config):
    for i, (fmap, filename) in enumerate(zip(feature_maps, filenames)):
        fmap.setPrimaryMSRunPath([filename.encode()])
        for feature in fmap:
            feature.setMetaValue("map_index", i)
    
    consensus_map = oms.ConsensusMap()
    column_headers = consensus_map.getColumnHeaders()
    for i, filename in enumerate(filenames):
        file_desc = oms.ColumnHeader()
        file_desc.filename = filename
        file_desc.label = filename.replace('.mzML', '')
        file_desc.size = feature_maps[i].size()
        column_headers[i] = file_desc
    consensus_map.setColumnHeaders(column_headers)  

    feature_grouper = oms.FeatureGroupingAlgorithmQT()
    params = feature_grouper.getDefaults()
    params.setValue("distance_MZ:unit", "ppm")  
    params.setValue("distance_MZ:max_difference", float(config["preprocessing"]["max_ppm_diff_combine_features"]))
    params.setValue("distance_RT:max_difference", float(config["preprocessing"]["max_rt_diff_combine_features"]))
    params.setValue("ignore_charge", "true")
    feature_grouper.setParameters(params)
    
    print("Running feature grouping algorithm")
    feature_grouper.group(feature_maps, consensus_map)
    consensus_map.setUniqueIds()
    return consensus_map

def save_consensus_csv_unfiltered(consensus_map, output_path):
    consensus_xml_file = os.path.join(output_path, 'consensus_features.consensusXML')
    oms.ConsensusXMLFile().store(consensus_xml_file, consensus_map)

    intensities = consensus_map.get_intensity_df()
    meta_data = consensus_map.get_metadata_df()[["rt", "mz", "quality"]]
    consensus_ids = [c_feat.getUniqueId() for c_feat in consensus_map]
    meta_data.insert(3,'consensus_id',consensus_ids)

    cm_df = pd.concat([meta_data, intensities], axis=1).reset_index(drop=True)
    cm_df.columns = [col.replace('.mzML', '') if '.mzML' in col else col for col in cm_df.columns]

    unfiltered_csv_path = os.path.join(output_path, 'consensus_unfiltered.csv')
    cm_df.to_csv(unfiltered_csv_path, index=False)
    print(f"Saved Unfiltered CSV to: {unfiltered_csv_path}")
    return unfiltered_csv_path

def main(config,arg_dic,mzML_files):        
    path = arg_dic["paths"]["mzml"]

    if not mzML_files:
        print("No files found!")
        return

    filenames = [os.path.basename(f) for f in mzML_files]
    all_features = []

    print(f"Processing {len(mzML_files)} files")

    output_path = arg_dic["paths"]["out"]
    os.makedirs(output_path, exist_ok=True)

    for f_path in mzML_files:
        fname = os.path.basename(f_path)
        print(f"Processing {fname}")
        
        exp = oms.MSExperiment()
        oms.MzMLFile().load(f_path, exp)
        exp.sortSpectra(True)

        mass_traces = detect_mass_traces(exp,config)
        mass_traces_final = detect_elution_peaks(mass_traces)
        features = detect_features_from_traces(mass_traces_final,config)
        
        all_features.append(features)
        save_features(features, fname, path, output_path)

        del exp
        del mass_traces
        del mass_traces_final
        gc.collect()


    all_features, ref_file, ref_index = align_feature_maps(all_features, filenames)
  
    consensus_map = link_features(all_features, filenames,config)

    del all_features
    gc.collect()
    
    unfiltered_csv_path = save_consensus_csv_unfiltered(consensus_map, output_path)
    print("Preprocessing Complete")

if __name__ == "__main__":
    main()
