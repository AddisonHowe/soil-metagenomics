import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from sklearn.feature_selection import mutual_info_regression
#from npeet_plus import mi, cmi
from scipy.spatial import KDTree


DATDIR = 'data'

"""
get_filepath
input: prefix (sample), file_type (abundance .bed or annotation .tsv)
output: a filepath string
"""
def get_filepath(prefix, file_type):
    base_path = 'data'
    if file_type == 'annotation':
        filename = f"raw_data/{prefix}.coassembly_ORFid_1stClusterDB_2ndClusterDB_KO_annotations_250316.tsv"
    elif file_type == 'abundance':
        filename = f"raw_data/{prefix}_all_samples_ORF_count_regions_rm0_ORF_ID_changed.bed"
    elif file_type == 'annotation_K00370':
        filename = f"subset_K00370/{prefix}.coassembly_annotations_K00370.tsv"
    elif file_type == 'abundance_K00370':
        filename = f"subset_K00370/{prefix}_all_samples_K00370.bed"
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
        FPATH = get_filepath(prefix, 'annotation_K00370')
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
    FPATH = get_filepath(prefix, 'annotation_K00370')
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
    metadata = pd.read_csv(f'../data/metadata.tsv', sep='\t')
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
    metadata = pd.read_csv(f'data/metadata.tsv', sep='\t')
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

def plot(data, title, cmap = 'Blues', vmin = None, vmax = None):
    soils = ['Soil3', 'Soil5', 'Soil6', 'Soil9', 'Soil11', 'Soil12', 'Soil14', 'Soil15', 'Soil16', 'Soil17']
    native = np.zeros(len(soils))
    for i, soil in enumerate(soils):
        native[i] = native_pH(soil)
 
    plt.imshow(data, aspect='auto', cmap=cmap, origin='lower', vmin = vmin, vmax = vmax)
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

        
import numpy as np
from sklearn.feature_selection import mutual_info_regression

def bootstrap_mi(X, Y, k=3, n_bootstrap=1000, confidence=0.95):
    """
    Compute bootstrapped mutual information using KSG estimator.
    
    Parameters:
    - X, Y: 1D arrays of shape (n_samples,) with your data.
    - k: k for kNN-based MI estimation (default=3).
    - n_bootstrap: Number of bootstrap resamples (default=1000).
    - confidence: Confidence interval width (default=0.95 for 95% CI).
    
    Returns:
    - mi_mean: Mean MI across bootstraps.
    - ci_low, ci_high: Confidence interval bounds.
    """
    X = np.array(X).reshape(-1, 1)  
    Y = np.array(Y).reshape(-1, 1)
    n_samples = X.shape[0]
    mi_boot = np.zeros(n_bootstrap)
    
    for i in range(n_bootstrap):
        # Resample with replacement
        idx = np.random.choice(n_samples, n_samples, replace=True)
        X_resampled = X[idx]
        Y_resampled = Y[idx]
        
        # Compute MI for this resample
        mi_boot[i] = mi(X_resampled, Y_resampled, k=k)
    
    # Compute summary statistics
    mi_mean = np.mean(mi_boot)
    ci_low = np.percentile(mi_boot, 100 * (1 - confidence) / 2)
    ci_high = np.percentile(mi_boot, 100 * (1 + confidence) / 2)
    
    return mi_mean, ci_low, ci_high


def information_2D(array1, array2, neighborhood_size=20, k=3):
    """
    Compute 2D mutual information map using KDTree for neighborhoods and NPEET for MI.
    
    Parameters:
    - array1, array2: 2D arrays (10x11) with your data
    - neighborhood_size: Target number of points to use for each MI estimate
    - k: k for kNN-based MI estimation
    
    Returns:
    - data: MI estimate at each point
    - data_low: Lower CI bound at each point
    - data_high: Upper CI bound at each point
    """
    height, width = array1.shape
    data = np.zeros((height, width))
    data_low = np.zeros((height, width))
    data_high = np.zeros((height, width))
    
    # Create coordinates and flatten arrays
    coords = np.array([[i, j] for i in range(height) for j in range(width)])
    values1 = array1.reshape(-1, 1)  
    values2 = array2.reshape(-1, 1)
    
    # Build KDTree for fast neighborhood queries
    tree = KDTree(coords)
    
    for i in range(height):
        for j in range(width):
            # Find nearest neighbors (including self)
            _, neighbor_indices = tree.query(
                [[i, j]], 
                k=neighborhood_size
            )
            
            # Extract local samples
            X_local = values1[neighbor_indices[0]]
            Y_local = values2[neighbor_indices[0]]
            
            # Compute MI with bootstrap CIs
            info, ci_low, ci_high = bootstrap_mi(X_local, Y_local, k=k)
            
            data[i][j] = info
            data_low[i][j] = ci_low
            data_high[i][j] = ci_high
            
    return data, data_low, data_high

    """
    Compute pairwise MI and CMI for 2D grids X, Y, Z using local neighborhoods.
    
    Args:
        X, Y, Z: 2D numpy arrays of shape (height, width) with continuous values.
        neighborhood_size: Number of nearest points to include in local estimation.
        k_neighbors: k for KSG estimation (n_neighbors).
        
    Returns:
        Dict containing MI and CMI grids (each a 2D array of same shape as X).
    """

def pairwise_cmi_analysis(X, Y, Z, neighborhood_size=50, k_neighbors=4):
    assert X.shape == Y.shape == Z.shape, "All input matrices must have the same shape"
    height, width = X.shape
    total_points = height * width
    
    # Flatten grids and get coordinates for spatial neighborhoods
    coords = np.array([[i, j] for i in range(height) for j in range(width)])
    values_X = X.reshape(-1)
    values_Y = Y.reshape(-1)
    values_Z = Z.reshape(-1)
    
    # Build KDTree for fast neighborhood queries
    tree = KDTree(coords)
    
    # Initialize output grids
    mi_xy_grid = np.zeros((height, width))
    mi_xz_grid = np.zeros_like(mi_xy_grid)
    mi_yz_grid = np.zeros_like(mi_xy_grid)
    cmi_xy_z_grid = np.zeros_like(mi_xy_grid)
    cmi_xz_y_grid = np.zeros_like(mi_xy_grid)
    cmi_yz_x_grid = np.zeros_like(mi_xy_grid)
    
    for idx in range(total_points):
        i, j = coords[idx]
        
        # Find nearest neighbors (including self)
        distances, neighbor_indices = tree.query(
            [i, j], 
            k=neighborhood_size
        )
        
        # Extract local samples
        x_local = values_X[neighbor_indices].reshape(-1, 1)
        y_local = values_Y[neighbor_indices]
        z_local = values_Z[neighbor_indices]
        
        # Compute MI and CMI (using NPEET)
        mi_xy_grid[i, j] = mi(x_local, y_local, k=k_neighbors)
        mi_xz_grid[i, j] = mi(x_local, z_local, k=k_neighbors)
        mi_yz_grid[i, j] = mi(y_local, z_local, k=k_neighbors)
        cmi_xy_z_grid[i, j] = mi(x_local, y_local, z = z_local, k=k_neighbors)
        cmi_xz_y_grid[i, j] = mi(x_local, z_local, z = y_local, k=k_neighbors)
        cmi_yz_x_grid[i, j] = mi(y_local, z_local, z = x_local, k=k_neighbors)
    
    return {
        'mi_xy': mi_xy_grid,
        'mi_xz': mi_xz_grid,
        'mi_yz': mi_yz_grid,
        'cmi_xy_z': cmi_xy_z_grid,
        'cmi_xz_y': cmi_xz_y_grid,
        'cmi_yz_x': cmi_yz_x_grid
    }

'''
get_info
input: two 2D arrays
output: (1) the mutual information between the two variables
        (2) 95% confidence error bar width from bootstrapping
        (3) a p value 
'''
def get_info(X, Y, permute = 1000, k = 6):
    X_flat = X.reshape(-1, 1)
    Y_flat = Y.reshape(-1, 1)
    info = mi(X_flat, Y_flat)
    mi_mean, ci_low, ci_high = bootstrap_mi(X_flat, Y_flat, k = k)
    error = (ci_high - ci_low)/2
    p = 0
    for i in range(permute):
        Y_shuffle = np.random.permutation(Y_flat)
        test_info = mi(X_flat, Y_shuffle)
        if test_info > info:
            p += 1
    p = (p + 1)/(permute + 1)
    return info, error, p 

'''
shuffle_Y_within_X_strata
used to perform shuffling within X bins within function get_info_control
'''
def shuffle_Y_within_X_strata(Y, X, n_bins=10):
    """Shuffle 2D Y within bins of 2D X."""
    # Flatten X and Y for binning
    X_flat = X.ravel()
    Y_flat = Y.ravel()
    
    # Bin X into quantiles (use 2D binning if spatial correlation matters)
    bins = np.quantile(X_flat, np.linspace(0, 1, n_bins + 1))
    X_binned = np.digitize(X_flat, bins[:-1])
    
    # Shuffle Y within each X-bin
    Y_shuffled_flat = Y_flat.copy()
    for bin_id in np.unique(X_binned):
        mask = (X_binned == bin_id)
        Y_shuffled_flat[mask] = np.random.permutation(Y_shuffled_flat[mask])
    
    # Reshape back to 2D
    return Y_shuffled_flat.reshape(Y.shape)

'''
get_info_control
inputs: X, Y, Z 
output: I(X ; Y | Z), error, p 
'''
def get_info_control(X, Y, Z, permute = 1000, k = 6):
    X_flat = X.reshape(-1, 1)
    Y_flat = Y.reshape(-1, 1)
    Z_flat = Z.reshape(-1, 1)
    info = cmi(X_flat, Y_flat, Z_flat)
    mi_mean, ci_low, ci_high = bootstrap_mi(X_flat, Y_flat, k = k)
    error = (ci_high - ci_low)/2
    p = 0
    for i in range(permute):
        Y_shuffle = shuffle_Y_within_X_strata(Y, Z)
        test_info = cmi(X_flat, Y_shuffle.reshape(-1, 1), z = Z_flat)
        if test_info > info:
            p += 1
    p = (p + 1)/(permute + 1)
    return info, error, p 