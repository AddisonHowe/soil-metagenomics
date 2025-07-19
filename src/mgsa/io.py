"""Input/Output functions

"""
import pandas as pd
import numpy as np


def get_filepath(prefix, file_type, KO = None, DATDIR = '../data'):
    """
    get_filepath
    input: prefix ('soil' or 'T0'), file_type (abundance .bed or annotation .tsv)
    output: a filepath string
    """
    if file_type == 'annotation' and KO is None:
        filename = f"raw_data/{prefix}.coassembly_ORFid_1stClusterDB_2ndClusterDB_KO_annotations_250316.tsv"
    elif file_type == 'abundance' and KO is None:
        filename = f"raw_data/{prefix}_all_samples_ORF_count_regions_rm0_ORF_ID_changed.bed"
    elif file_type == 'annotation':
        filename = f"subset_{KO}/{prefix}.coassembly_annotations_{KO}.tsv"
    elif file_type == f'abundance':
        filename = f"subset_{KO}/{prefix}_all_samples_{KO}.bed"
    else:
        raise ValueError("file_type must be either 'annotation' or 'abundance'")
    
    return f'{DATDIR}/{filename}'


def pH_sample(sample, DATDIR = '../data'):
    """
    pH_sample
    input: sample ID
    output: the perturbed pH
    """

    metadata = pd.read_csv(f'{DATDIR}/metadata.tsv', sep='\t')
    metadata = metadata.set_index('sample')
    return metadata.loc[sample, 'pH']


def samples_from_soils(soil, drug = 'None', DATDIR = '../data'):
    """
    samples
    input: soil name (eg 'Soil3') and drug (default: None, 'CHL' other option)
    output: a list of sample IDs corresponding to that soil, in order of increasing pH
    """
    
    sample_IDs = []
    metadata = pd.read_csv(f'{DATDIR}/metadata.tsv', sep='\t')
    metadata = metadata.set_index('sample')
    soil_ = soil + '_'
    pHs = []
    for sample in metadata.index:
        if drug in sample and 'T9' in sample and soil_ in sample and 'Nitrate' not in sample:
            sample_IDs.append(sample)
            pHs.append(pH_sample(sample, DATDIR))
            
    #make sure the samples are organized by increaing perturbed pH
    pHs = np.array(pHs)
    sample_IDs = np.array(sample_IDs)
    
    indices = np.argsort(pHs)
    
    sample_IDs = sample_IDs[indices]
    sample_IDs = sample_IDs.tolist()
    
    
    return sample_IDs


def find_orfs(file_path, ko_number):
    """
    find_orfs
    input: filepath, KO number
    output: a list of unclostered ORFs which correspond to that KO in a given sample
    """
    df = pd.read_csv(file_path, sep='\t', header=None)
    matching_orfs = df[df[3] == ko_number][0].tolist()
    return matching_orfs

def find_orfs_from_cluster(cluster, KO, DATDIR):
    """
    find_orfs_from_cluster
    input: a cluster ID
    output: all corresponding ORFs 
    """
    prefixes = ['T0', 'Soil3', 'Soil5', 'Soil6', 'Soil9', 'Soil11', 'Soil12', 'Soil14', 'Soil15', 'Soil16', 'Soil17']
    orfs = []
    for prefix in prefixes: 
        FPATH = get_filepath(prefix, 'annotation', KO, DATDIR)
        df = pd.read_csv(FPATH, sep='\t', header=None)
        temp = df[df[2] == cluster][0].tolist()
        orfs.extend(temp)
    return orfs

def find_cluster_from_orf(orf, KO = 'K00370', DATDIR = '../data'):
    """
    find_cluster_from_orf
    input: a cluster ID
    output: all corresponding ORFs 
    """
    prefix = orf.split('.')[0]
    FPATH = get_filepath(prefix, 'annotation', KO, DATDIR)
    df = pd.read_csv(FPATH, sep='\t', header=None)
    cluster = df[df[0] == orf][2].tolist()
    return cluster[0]

def perturbed_pHs(soil, DATDIR = '../data'):
    """
    perturbed_pHs
    given a soil id in the list Soil3, 5, 6, 9, 11, 12, 14, 15, 16, 17
    output: an array containing the perturbed pHs from that soil
    """
    pHs = []
    metadata = pd.read_csv(f'{DATDIR}/metadata.tsv', sep='\t')
    metadata = metadata.set_index('sample')
    
    soil_ = soil + '_'
    for sample in metadata.index:
        if 'None' in sample and 'T9' in sample and soil_ in sample and 'Nitrate' not in sample:
            pHs.append(metadata.loc[sample, 'pH'])
    
    return sorted(pHs)