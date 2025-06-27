import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

DATDIR = '/projects/p32818/metagenomic_data/data'

"""
get_filepath
input: prefix (sample), file_type (abundance .bed or annotation .tsv)
output: a filepath string
"""
def get_filepath(prefix, file_type):
    base_path = '/projects/p32818/metagenomic_data/data/raw_data'
    if file_type == 'annotation':
        filename = f"{prefix}.coassembly_ORFid_1stClusterDB_2ndClusterDB_KO_annotations_250316.tsv"
    elif file_type == 'abundance':
        filename = f"{prefix}_all_samples_ORF_count_regions_rm0_ORF_ID_changed.bed"
    else:
        raise ValueError("file_type must be either 'annotation' or 'abundance'")
    
    return f'{base_path}/{filename}'

"""
get_clustered_orfs
input: filepath, KO number
output: a list of cluster IDs which correspond to that KO in a given sample
"""
def find_clustered_orfs(file_path, ko_number):
    df = pd.read_csv(file_path, sep='\t', header=None)
    matching_orfs = df[df[3] == ko_number][2].tolist()
    return matching_orfs

"""
find_orfs
input: filepath, KO number
output: a list of unclostered ORFs which correspond to that KO in a given sample
"""
def find_orfs(file_path, ko_number):
    df = pd.read_csv(file_path, sep='\t', header=None)
    matching_orfs = df[df[3] == ko_number][0].tolist()
    return matching_orfs

"""
find_orfs_from_cluster
input: a cluster ID
output: all corresponding ORFs 
"""
def find_orfs_from_cluster(cluster):
    prefixes = ['T0', 'Soil3', 'Soil5', 'Soil6', 'Soil9', 'Soil11', 'Soil12', 'Soil14', 'Soil15', 'Soil16', 'Soil17']
    orfs = []
    for prefix in prefixes: 
        FPATH = get_filepath(prefix, 'annotation')
        df = pd.read_csv(FPATH, sep='\t', header=None)
        temp = df[df[2] == cluster][0].tolist()
        orfs.extend(temp)
    return orfs

"""
find_cluster_from_orf
input: a cluster ID
output: all corresponding ORFs 
"""
def find_cluster_from_orf(orf):
    prefix = orf.split('.')[0]
    FPATH = get_filepath(prefix, 'annotation')
    df = pd.read_csv(FPATH, sep='\t', header=None)
    cluster = df[df[0] == orf][2].tolist()
    return cluster[0]

"""
perturbed_pHs
given a soil id in the list Soil3, 5, 6, 9, 11, 12, 14, 15, 16, 17
output: an array containing the perturbed pHs from that soil
"""
def perturbed_pHs(soil):
    pHs = []
    metadata = pd.read_csv(f'{DATDIR}/metadata.tsv', sep='\t')
    metadata = metadata.set_index('sample')
    
    soil_ = soil + '_'
    for sample in metadata.index:
        if 'None' in sample and 'T9' in sample and soil_ in sample and 'Nitrate' not in sample:
            pHs.append(metadata.loc[sample, 'pH'])
    
    return sorted(pHs)

"""
samples_from_soils
given a soil id in the list Soil3, 5, 6, 9, 11, 12, 14, 15, 16, 17
output: an array containing the sample IDs for that soil after perturbation
"""
def samples_from_soils(soil):
    sample_IDs = []
    metadata = pd.read_csv(f'{DATDIR}/metadata.tsv', sep='\t')
    metadata = metadata.set_index('sample')
    
    soil_ = soil + '_'
    for sample in metadata.index:
        if 'None' in sample and 'T9' in sample and soil_ in sample and 'Nitrate' not in sample:
            sample_IDs.append(sample)
    
    return sample_IDs

"""
native_pH
input: any soil id 
output: the native pH of that soil
"""
def native_pH(soil):
    metadata = pd.read_csv(f'{DATDIR}/metadata.tsv', sep='\t')
    metadata = metadata.set_index('sample')
    for sample in metadata.index:
        soil_= soil + '_'
        if 'T0' in sample and soil_ in sample:
            pH = metadata.loc[sample, 'pH']
            
    return pH

'''
enriched_native_pH
inputs: cluster ID, filtered data, cluseter list
output: soil which is most enriched at T0
'''

def enriched_native_pH(CID, data, cluster_IDs):
    soils = ['Soil3', 'Soil5', 'Soil6', 'Soil9', 'Soil11', 'Soil12', 'Soil14', 'Soil15', 'Soil16', 'Soil17']
    CIDX = cluster_IDs.index(CID)
    soil_idx = np.argmax(data[CIDX])
    return soils[soil_idx]
