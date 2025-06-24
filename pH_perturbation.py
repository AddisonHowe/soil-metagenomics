import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys

DATDIR = '/projects/p32818/metagenomic_data/data'

"""
get_filepath
input: prefix (sample), file_type (abundance .bed or annotation .tsv)
output: a filepath string
"""
def get_filepath(prefix, file_type):
    base_path = '/projects/p32818/metagenomic_data/data'
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
    
    for sample in metadata.index:
        if 'None' in sample and 'T9' in sample and soil in sample and 'Nitrate' not in sample:
            pHs.append(metadata.loc[sample, 'pH'])
    
    return pHs

"""
samples_from_soils
given a soil id in the list Soil3, 5, 6, 9, 11, 12, 14, 15, 16, 17
output: an array containing the sample IDs for that soil after perturbation
"""
def samples_from_soils(soil):
    sample_IDs = []
    metadata = pd.read_csv(f'{DATDIR}/metadata.tsv', sep='\t')
    metadata = metadata.set_index('sample')
    
    for sample in metadata.index:
        if 'None' in sample and 'T9' in sample and soil in sample and 'Nitrate' not in sample:
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
        if 'T0' in sample and soil in sample:
            pH = metadata.loc[sample, 'pH']
            
    return pH


import sys

def main():
    

    soil = sys.argv[1]
    print(f"Received input: {soil}")
    
    
    ORFs = find_orfs(get_filepath(soil, 'annotation'), 'K00370')

    #Unique_IDs is the array, defined earlier, containing the list of cluster IDs

    sample_list = samples_from_soils(soil)

    metadata = pd.read_csv(f'{DATDIR}/metadata.tsv', sep = '\t')
    metadata = metadata.set_index('sample')
    cluster_IDs = pd.read_csv('out/cluster_ids.tsv', sep = '\t')
    cluster_IDs = cluster_IDs.values
    cluster_IDs = cluster_IDs.tolist()

    chunk_size = 100000

    data = np.zeros((len(cluster_IDs), 11)) #data for plot stored here, each row is a dinstinct cluster

    for chunk in pd.read_csv(get_filepath(soil, 'abundance'), sep='\s+', header=None,  chunksize = chunk_size):
        filtered_chunk = chunk[chunk.iloc[:, 1].isin(ORFs)]
        for i in range(len(filtered_chunk)):
            sample_id = filtered_chunk.iloc[i, 0]
            if sample_id in sample_list:
                orf = filtered_chunk.iloc[i, 1]
                rel_abundance = filtered_chunk.iloc[i, 2]
                spikein = metadata.loc[sample_id, 'spikein_sum']
                cluster = find_cluster_from_orf(orf)
                if cluster in cluster_IDs:
                    row_idx = cluster_IDs.index(find_cluster_from_orf(orf))
                    col_idx = sample_list.index(sample_id)
                    data[row_idx, col_idx] += rel_abundance/spikein
                    print(rel_abundance/spikein)
    
    print(data)
    np.savetxt(f"out/{soil}data.tsv", data, delimiter = '\t', fmt = '%0.6f')
        
        
if __name__ == "__main__":
    main()
