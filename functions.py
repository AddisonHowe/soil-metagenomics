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
perturbed_PH_sample
given a sample ID
output: the perturbed pH
"""

def perturbed_pH_sample(sample):
    metadata = pd.read_csv(f'{DATDIR}/metadata.tsv', sep='\t')
    metadata = metadata.set_index('sample')
    return metadata.loc[sample, 'pH']

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

'''
plot_native_perturbed
inputs: an array with 10 rows (native pHs) and 11 columns (perturbed pHs)
produces plot
'''

def plot(data, title):
    soils = ['Soil3', 'Soil5', 'Soil6', 'Soil9', 'Soil11', 'Soil12', 'Soil14', 'Soil15', 'Soil16', 'Soil17']
    native = np.zeros(len(soils))
    for i, soil in enumerate(soils):
        native[i] = native_pH(soil)
 
    plt.imshow(data, aspect='auto', cmap='Blues', origin='lower')
    plt.colorbar()

    plt.xlabel('Perturbed pH (approximate)')
    plt.ylabel('Native pH')

    x = np.linspace(3.8, 8.4, 11)
    plt.xticks(ticks=np.linspace(0, 10, 11), labels=[f"{val:.1f}" for val in x], rotation=45)
    plt.yticks(ticks=np.linspace(0, 9, 10), labels=[f"{val:.1f}" for val in native], rotation=45)


    plt.title(label=title)


    plt.tight_layout()
    plt.show()
    
'''
information_1D
input: 1 array with with 10 rows (native pHs) and 11 columns (perturbed pHs)
output: 2 1D arrays one for mutual information across native pH, one mutual information across perturbed pH
'''
def information_1D(array):
    
    from sklearn.feature_selection import mutual_info_regression

    soils = ['Soil3', 'Soil5', 'Soil6', 'Soil9', 'Soil11', 'Soil12', 'Soil14', 'Soil15', 'Soil16', 'Soil17']

    native = np.zeros(len(soils))
    for i, soil in enumerate(soils):
        native[i] = native_pH(soil)
        
    perturbed = np.linspace(3.8, 8.4, 11)
    
    #first, MI between native pH and array
    
    x = native.reshape(-1, 1)
    MI_native_and_array = np.zeros(11)
    Y = array
    for i in range(len(MI_native_and_array)):
        
        
        mi = mutual_info_regression(x, Y[:, i], random_state=0) 
        
        MI_native_and_array[i] = mi 
    
    x = np.array(perturbed).reshape(-1,1)
    MI_perturbed_and_array = np.zeros(10)
    for i in range(len(MI_perturbed_and_array)):
        

        mi = mutual_info_regression(x, Y[i], random_state=0) 
        
        MI_perturbed_and_array[i] = mi 
        
    return MI_native_and_array, MI_perturbed_and_array

'''
information_1D_2
input: 2 arrays with with 10 rows (native pHs) and 11 columns (perturbed pHs)
output: 2 1D arrays one for mutual information across native pH, one mutual information across perturbed pH
'''
def information_1D_2(array1, array2):
    
    from sklearn.feature_selection import mutual_info_regression

    soils = ['Soil3', 'Soil5', 'Soil6', 'Soil9', 'Soil11', 'Soil12', 'Soil14', 'Soil15', 'Soil16', 'Soil17']

    native = np.zeros(len(soils))
    for i, soil in enumerate(soils):
        native[i] = native_pH(soil)
        
    perturbed = np.linspace(3.8, 8.4, 11)
    
    MI_vs_perturbed = np.zeros(11)
    MI_vs_native = np.zeros(10)
    
    for i in range(len(MI_vs_native)):
        
        X = array1[i].reshape(-1, 1)
        Y = array2[i]
        
        mi = mutual_info_regression(X, Y, random_state=0) 
        
        MI_vs_native[i] = mi
        
    for i in range(len(MI_vs_perturbed)):
        
        X = array1[:,i].reshape(-1, 1)
        Y = array2[:,i]
        
        mi = mutual_info_regression(X, Y, random_state=0)
        
        MI_vs_perturbed[i] = mi
        
    return MI_vs_perturbed, MI_vs_native

        

'''
information_2D
input: 2 arrays with with 10 rows (native pHs) and 11 columns (perturbed pHs)
output: a new array with the approximate mutual information at each point. 
Using a binning approach
'''

def information_2D(array1, array2):
    from sklearn.feature_selection import mutual_info_regression
    data = np.zeros((10, 11))
    for i in range(10):
        for j in range(11):
            
            #if in corner, average mi from only 4 points
            if (i == 0 and j == 0):
                X = [array1[0][0], array1[1][0], array1[1][1], array1[0][1]]
                Y = [array2[0][0], array2[1][0], array2[1][1], array2[0][1]]
                
            elif (i == 9 and j == 0):
                X = [array1[9][0], array1[8][0], array1[8][1], array1[9][1]]
                Y = [array2[9][0], array2[8][0], array2[8][1], array2[9][1]]
                
            elif (i == 9 and j == 10):
                X = [array1[9][9], array1[8][9], array1[8][10], array1[9][10]]
                Y = [array2[9][9], array2[8][9], array2[8][10], array2[9][10]]
                
            elif (i == 0 and j == 10):
                X = [array1[0][9], array1[1][9], array1[1][10], array1[0][10]]
                Y = [array2[0][9], array2[1][9], array2[1][10], array2[0][10]]
                
            #if on a side, average 6 points
            elif (i == 0):
                X = [array1[0][j], array1[0][j+1], array1[0][j-1], array1[1][j], array1[1][j+1], array1[1][j-1]]
                Y = [array2[0][j], array2[0][j+1], array2[0][j-1], array2[1][j], array2[1][j+1], array2[1][j-1]]
                
            elif (i == 9):
                X = [array1[9][j], array1[9][j+1], array1[9][j-1], array1[8][j], array1[8][j+1], array1[8][j-1]]
                Y = [array2[9][j], array2[9][j+1], array2[9][j-1], array2[8][j], array2[8][j+1], array2[8][j-1]]
                
            elif (j == 0):
                X = [array1[i][0], array1[i-1][0], array1[i+1][0], array1[i][1], array1[i-1][1], array1[i+1][1]]
                Y = [array2[i][0], array2[i-1][0], array2[i+1][0], array2[i][1], array2[i-1][1], array2[i+1][1]]
                
            elif (j == 10):
                X = [array1[i][10], array1[i-1][10], array1[i+1][10], array1[i][9], array1[i-1][9], array1[i+1][9]]
                Y = [array2[i][10], array2[i-1][10], array2[i+1][10], array2[i][9], array2[i-1][9], array2[i+1][9]]
                
            #if in the middle, average 9 points
            else:
                X = [array1[i][j], array1[i + 1][j], array1[i - 1][j], array1[i][j+1], array1[i][j-1], array1[i+1][j+1], array1[i+1][j-1], array1[i-1][j+1], array1[i-1][j-1]]
                Y = [array2[i][j], array2[i + 1][j], array2[i - 1][j], array2[i][j+1], array2[i][j-1], array2[i+1][j+1], array2[i+1][j-1], array2[i-1][j+1], array2[i-1][j-1]]
                
                
            X = np.array(X).reshape(-1, 1)
            Y = np.array(Y)
            
            mi = mutual_info_regression(X, Y)
            data[i][j] = mi
    return data 