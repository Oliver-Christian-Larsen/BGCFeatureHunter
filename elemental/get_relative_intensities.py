import pandas as pd
import numpy as np
import glob
import os
import pyopenms as oms


def get_files(arg_dic):

    f_xml_files = glob.glob(os.path.join(arg_dic["paths"]["out"],"featureXML/","*.featureXML"))

    cf_xml_file = os.path.join(arg_dic["paths"]["out"],"consensus_features.consensusXML") 

    feature_csv = pd.read_csv(arg_dic["paths"]["elemental_featurefile"], dtype={'consensus_id': str})

    return f_xml_files, cf_xml_file, feature_csv


def make_lookup(cf_xml_file):

    mapping_dic = {}

    cm = oms.ConsensusMap()
    oms.ConsensusXMLFile().load(cf_xml_file,cm)

    headers = cm.getColumnHeaders()

    for map_id, header in headers.items():
        file_name_no_ext = header.filename
        filenameXML = file_name_no_ext.replace(".mzML",".featureXML")
        
        mapping_dic[map_id] = filenameXML

    return mapping_dic, cm


def feature_iteration(consensus_map, feature_csv):
    consensus_ids = feature_csv['consensus_id'].dropna().apply(int)
    
    consensus_dict = {feature.getUniqueId(): feature for feature in consensus_map}
    extracted_data = []

    for cid in consensus_ids:
        if cid in consensus_dict:
            target_feature = consensus_dict[cid]
            
            map_feature_mapping = {}
            
            for handle in target_feature.getFeatureList():
                map_index = handle.getMapIndex()
                original_feature_id = handle.getUniqueId()
                map_feature_mapping[map_index] = original_feature_id

            extracted_data.append({
            'consensus_id': cid,
            'map_mapping': map_feature_mapping
        })

    return extracted_data

def build_feature_index(mapping_dic, f_xml_files):

    indexed_maps = {}

    for map_id, expected_filename in mapping_dic.items():
        target_path = next((path for path in f_xml_files if expected_filename in path), None)

        fmap = oms.FeatureMap()

        oms.FeatureXMLFile().load(target_path, fmap)

        feature_dict = {feature.getUniqueId(): feature for feature in fmap}

        indexed_maps[map_id] = feature_dict

    return indexed_maps

def safe_weighted_avg(ratios, weights_list):
    valid_idx = ~np.isnan(ratios)
    if not np.any(valid_idx):
        return np.nan
    
    valid_ratios = np.array(ratios)[valid_idx]
    valid_weights = np.array(weights_list)[valid_idx]
    
    if np.sum(valid_weights) == 0:
        return np.nan
        
    return np.average(valid_ratios, weights=valid_weights)


def compute_masstrace_ratios(feature_translation_dict, master_index):
    final_data = []

    for entry in feature_translation_dict:
        cid = str(entry['consensus_id'])
        map_mapping = entry['map_mapping']
        
        ratios_m0 = []
        ratios_m1 = []
        ratios_m2 = []
        ratios_m3 = []
        mz_m0 = []
        mz_m1 = []
        mz_m2 = []
        mz_m3 = []
        weights = []
        charge = []

        for map_idx, original_feat_id in map_mapping.items():
    
            if map_idx in master_index and original_feat_id in master_index[map_idx]:
                target_feature = master_index[map_idx][original_feat_id]
                
                if target_feature.metaValueExists("num_of_masstraces"):
                    num_traces = int(target_feature.getMetaValue("num_of_masstraces"))
                    
                    if num_traces > 0 and target_feature.metaValueExists("masstrace_intensity"):
                        weights.append(target_feature.getIntensity())
                        charge.append(target_feature.getCharge())

                        intensities = list(target_feature.getMetaValue("masstrace_intensity"))
                        mz_list = list(target_feature.getMetaValue("masstrace_centroid_mz"))
                        
                        if len(intensities) > 0 and intensities[0] > 0:
                            ratios_m0.append(1)
                            mz_m0.append(mz_list[0]) 
                        else:
                            ratios_m0.append(np.nan)
                            mz_m0.append(np.nan)

                        if len(intensities) > 1 and intensities[1] > 0 and intensities[0] > 0:    
                            ratios_m1.append(intensities[1] / intensities[0])
                            mz_m1.append(mz_list[1])
                        else:
                            ratios_m1.append(np.nan)
                            mz_m1.append(np.nan)

                        if len(intensities) > 2 and intensities[2] > 0 and intensities[0] > 0:
                            ratios_m2.append(intensities[2] / intensities[0])
                            mz_m2.append(mz_list[2])
                        else:
                            ratios_m2.append(np.nan)
                            mz_m2.append(np.nan)
                            
                        if len(intensities) > 3 and intensities[3] > 0 and intensities[0] > 0:
                            ratios_m3.append(intensities[3] / intensities[0])
                            mz_m3.append(mz_list[3])
                        else:
                            ratios_m3.append(np.nan)
                            mz_m3.append(np.nan)


        final_data.append({
            'consensus_id': cid,
            'charge' : round(safe_weighted_avg(charge, weights), 0),
            'mz_m0': safe_weighted_avg(mz_m0, weights),
            'ratio_m0': safe_weighted_avg(ratios_m0, weights),

            'mz_m1': safe_weighted_avg(mz_m1, weights),
            'ratio_m1': safe_weighted_avg(ratios_m1, weights),

            'mz_m2': safe_weighted_avg(mz_m2, weights),
            'ratio_m2': safe_weighted_avg(ratios_m2, weights),

            'mz_m3': safe_weighted_avg(mz_m3, weights),
            'ratio_m3': safe_weighted_avg(ratios_m3, weights),
        })

    return final_data


def main(arg_dic):
    f_xml_files, cf_xml_file, feature_csv = get_files(arg_dic)

    mapping_dic, consensus_map = make_lookup(cf_xml_file)

    feature_translation_dict = feature_iteration(consensus_map, feature_csv)

    master_index = build_feature_index(mapping_dic, f_xml_files)

    ratio_data = compute_masstrace_ratios(feature_translation_dict, master_index)

    return ratio_data
 

if __name__ == "__main__":
    main(arg_dic)

