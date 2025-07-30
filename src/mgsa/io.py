"""Input/Output functions

"""
import pandas as pd
import numpy as np
from Bio.PDB import PDBParser
from Bio.PDB.Structure import Structure


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


def load_pdb_structure(fpath: str, id: str) -> Structure:
    """Load a PDB structure from a pdb file.
    
    Args:
        fpath (str): path to pdb file.
        id (str): the id for the returned structure.
    
    Returns:
        (Structure) Protein structure.
    """
    parser = PDBParser()
    structure = parser.get_structure(id, fpath)
    return structure


def get_data(orf, DATDIR = '../out', KO = 'K00370', drug = 'None', map = '09', ):
    """
    Args:
        orf (string): an orf
        DATDIR = path to out directory
        KO = KO number (eg 'K00370')
        drug = 'None' or 'CHL'
        map = '09' or '08'
    Returns:
        (1) a 20-length array with T0 data
        (2) a 10-row, 11-column array with T9 data 
    """
    soils = ['Soil3', 'Soil5', 'Soil6', 'Soil9', 'Soil11', 'Soil12', 'Soil14', 'Soil15', 'Soil16', 'Soil17'] 
    if KO == 'K02567':
        soils = ['Soil3', 'Soil5', 'Soil9', 'Soil11', 'Soil12', 'Soil14', 'Soil15', 'Soil16', 'Soil17'] 
        
    id_list = pd.read_csv(f'{DATDIR}/orf_ids/cluster_ids_{map}_{KO}.tsv', sep = '\t', header = None)
    id_list = id_list.values
    id_list = [id[0] for id in id_list]

    idx = id_list.index(orf)

    T0data_all = pd.read_csv(f"{DATDIR}/{KO}abundances/T0data_{map}_{drug}_{KO}.tsv", sep='\t', header=None)
    T0data_all = T0data_all.values
    T0data = T0data_all[idx]

    T9data = []
    for soil in soils: 
        data_all = pd.read_csv(f"{DATDIR}/{KO}abundances/{soil}data_{map}_{drug}_{KO}.tsv", sep='\t', header=None)
        data_all = data_all.values
        T9data.append(data_all[idx])
        
    return T0data, T9data

def get_function(sample, DATDIR = '../data'):
    """
    Args:
        sample (string): a string name of a sample
    Outputs:
        an array of floats, 6 rows by 10 columns
        The first three rows are for replicates 1 through 3, NO3 concentraiton in mM
        The last three rows are for replicates 1 through 3, NO2 concentration in mM
    """
    
    data= pd.read_csv(f'{DATDIR}/function.tsv', keep_default_na=False, sep='\t')
    output = np.zeros((6, 10))

    parts = sample.split('_')

    soil = parts[0]

    if 'N' in parts[3]:
        unit = int(parts[4]) 
    else:
        unit = int(parts[3])

    drug = parts[-2]
    print('drug', drug)
    print('unit', unit)
    print('soil', soil)


    for i in [1,2,3]:
        temp1 = data[(data['Chloramphenicol'] == drug) & (data['Unit'] == unit) & (data['Soil'] == soil) & (data['Replicate'] == i)]
        temp2 = temp1['NO3_mM']
        temp2 = temp2.tolist()
        output[i-1] = temp2
        temp2 = temp1['NO2_mM']
        temp2 = temp2.tolist()
        output[i+2] = temp2
        
    return output

