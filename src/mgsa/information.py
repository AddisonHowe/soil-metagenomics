import numpy as np
from sklearn.feature_selection import mutual_info_regression
from npeet_plus import mi, cmi, entropy
from scipy.spatial import KDTree
import pandas as pd
from sklearn.preprocessing import StandardScaler

def native_pH(soil, DATDIR = '../data'):
    metadata = pd.read_csv(f'{DATDIR}/metadata.tsv', sep='\t')
    metadata = metadata.set_index('sample')
    for sample in metadata.index:
        soil_= soil + '_'
        if 'T0' in sample and soil_ in sample:
            pH = metadata.loc[sample, 'pH']
            
    return pH

def normalize_data(X, Y=None):
    """
    Normalize data for kNN-based information estimation.
    If Y is provided, fit scaler on X and transform both X and Y.
    """
    X = np.array(X).reshape(-1, 1)
    
    if Y is not None:
        Y = np.array(Y).reshape(-1, 1)
        # For MI, normalize X and Y independently
        X_norm = StandardScaler().fit_transform(X)
        Y_norm = StandardScaler().fit_transform(Y)
        return X_norm, Y_norm
    else:
        # For entropy, just normalize X
        return StandardScaler().fit_transform(X)

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
    # Normalize the data first
    X_norm, Y_norm = normalize_data(X, Y)
    n_samples = X_norm.shape[0]
    mi_boot = np.zeros(n_bootstrap)
    
    for i in range(n_bootstrap):
        # Resample with replacement
        idx = np.random.choice(n_samples, n_samples, replace=True)
        X_resampled = X_norm[idx]
        Y_resampled = Y_norm[idx]
        
        # Compute MI for this resample
        mi_boot[i] = mi(X_resampled, Y_resampled, k=k)
    
    # Compute summary statistics
    mi_mean = np.mean(mi_boot)
    ci_low = np.percentile(mi_boot, 100 * (1 - confidence) / 2)
    ci_high = np.percentile(mi_boot, 100 * (1 + confidence) / 2)
    
    return mi_mean, ci_low, ci_high

def information_2D(array1, array2, neighborhood_size=20, k=3, n_bootstrap=1000):
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
    entropy_map = np.zeros((height, width))
    
    # Create coordinates and flatten arrays
    coords = np.array([[i, j] for i in range(height) for j in range(width)])
    
    # Normalize the input arrays
    values1_norm = normalize_data(array1)
    values2_norm = normalize_data(array2)
    
    # Build KDTree for fast neighborhood queries
    tree = KDTree(coords)
    
    for i in range(height):
        for j in range(width):
            # Find nearest neighbors (including self)
            _, neighbor_indices = tree.query(
                [[i, j]], 
                k=neighborhood_size
            )
            
            # Extract local samples (already normalized)
            X_local = values1_norm[neighbor_indices[0]]
            Y_local = values2_norm[neighbor_indices[0]]
            
            # Compute MI with bootstrap CIs
            info, ci_low, ci_high = bootstrap_mi(X_local, Y_local, k=k, n_bootstrap=n_bootstrap)
            
            e1 = entropy(X_local, k=k)
            e2 = entropy(Y_local, k=k)
            entropy_map[i, j] = max(e1, e2)
            data[i][j] = info
            data_low[i][j] = ci_low
            data_high[i][j] = ci_high
    
    with np.errstate(divide='ignore', invalid='ignore'):
        efficiency_map = data / (entropy_map + 1e-5)
        efficiency_map = np.clip(efficiency_map, 0, 1)    
        efficiency_map[np.isnan(efficiency_map)] = 0.0
            
    return data, data_low, data_high, efficiency_map, entropy_map

def pairwise_cmi_analysis(X, Y, Z, neighborhood_size=50, k_neighbors=4):
    """
    Compute pairwise MI and CMI for 2D grids X, Y, Z using local neighborhoods.
    
    Args:
        X, Y, Z: 2D numpy arrays of shape (height, width) with continuous values.
        neighborhood_size: Number of nearest points to include in local estimation.
        k_neighbors: k for KSG estimation (n_neighbors).
        
    Returns:
        Dictionary containing MI and CMI grids (each a 2D array of same shape as X).
    """
    assert X.shape == Y.shape == Z.shape, "All input matrices must have the same shape"
    height, width = X.shape
    total_points = height * width
    
    # Normalize all input arrays
    X_norm = normalize_data(X).reshape(height, width)
    Y_norm = normalize_data(Y).reshape(height, width)
    Z_norm = normalize_data(Z).reshape(height, width)
    
    # Flatten grids and get coordinates for spatial neighborhoods
    coords = np.array([[i, j] for i in range(height) for j in range(width)])
    values_X = X_norm.reshape(-1)
    values_Y = Y_norm.reshape(-1)
    values_Z = Z_norm.reshape(-1)
    
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
        
        # Extract local samples (already normalized)
        x_local = values_X[neighbor_indices].reshape(-1, 1)
        y_local = values_Y[neighbor_indices]
        z_local = values_Z[neighbor_indices]
        
        # Compute MI and CMI (using NPEET)
        mi_xy_grid[i, j] = mi(x_local, y_local, k=k_neighbors)
        mi_xz_grid[i, j] = mi(x_local, z_local, k=k_neighbors)
        mi_yz_grid[i, j] = mi(y_local, z_local, k=k_neighbors)
        cmi_xy_z_grid[i, j] = mi(x_local, y_local, z=z_local, k=k_neighbors)
        cmi_xz_y_grid[i, j] = mi(x_local, z_local, z=y_local, k=k_neighbors)
        cmi_yz_x_grid[i, j] = mi(y_local, z_local, z=x_local, k=k_neighbors)
    
    return {
        'mi_xy': mi_xy_grid,
        'mi_xz': mi_xz_grid,
        'mi_yz': mi_yz_grid,
        'cmi_xy_z': cmi_xy_z_grid,
        'cmi_xz_y': cmi_xz_y_grid,
        'cmi_yz_x': cmi_yz_x_grid
    }

def shuffle_Y_within_X_strata(Y, X, n_bins=10):
    '''
    shuffle_Y_within_X_strata
    used to perform shuffling within X bins within function get_info_control
    '''
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

def get_info(X, Y, permute=1000, k=4):
    '''
    get_info
    input: two 2D arrays
    output: (1) the mutual information between the two variables
            (2) 95% confidence error bar width from bootstrapping
            (3) a p value 
    '''
    # Normalize both variables
    X_norm, Y_norm = normalize_data(X, Y)
    
    info = mi(X_norm, Y_norm, k=k)
    mi_mean, ci_low, ci_high = bootstrap_mi(X_norm, Y_norm, k=k)
    error = (ci_high - ci_low) / 2
    p = 0
    
    for i in range(permute):
        Y_shuffle = np.random.permutation(Y_norm)
        test_info = mi(X_norm, Y_shuffle, k=k)
        if test_info > info:
            p += 1
            
    p = (p + 1) / (permute + 1)
    return info, error, p

def get_entropy(X, k=4):
    '''
    get_entropy
    input: 2D array
    output: entropy estimate in nats
    '''
    X_norm = normalize_data(X)
    e = entropy(X_norm, k=k)
    return e

def get_info_control(X, Y, Z, permute=1000, k=4):
    '''
    get_info_control
    inputs: X, Y, Z 
    output: I(X ; Y | Z), error, p 
    '''
    # Normalize all three variables
    X_norm = normalize_data(X)
    Y_norm = normalize_data(Y)
    Z_norm = normalize_data(Z)
    
    info = cmi(X_norm, Y_norm, Z_norm, k=k)
    mi_mean, ci_low, ci_high = bootstrap_mi(X_norm, Y_norm, k=k)
    error = (ci_high - ci_low) / 2
    p = 0
    
    for i in range(permute):
        Y_shuffle = shuffle_Y_within_X_strata(Y_norm.reshape(Y.shape), Z_norm.reshape(Z.shape))
        test_info = cmi(X_norm, Y_shuffle.reshape(-1, 1), z=Z_norm, k=k)
        if test_info > info:
            p += 1
            
    p = (p + 1) / (permute + 1)
    return info, error, p

def normalized_info_2D(array1, array2, neighborhood_size=20, k=3):
    """
    Compute normalized 2D mutual information map I(x;y)/H(y) using KDTree for neighborhoods.
    
    Parameters:
    - array1, array2: 2D arrays (10x11) with your data
    - neighborhood_size: Target number of points to use for each MI estimate
    - k: k for kNN-based MI estimation
    
    Returns:
    - normalized_mi: Normalized MI estimate I(x;y)/H(y) at each point
    - raw_mi: Raw mutual information I(x;y) at each point  
    - entropy_y: Entropy H(y) at each point
    """
    height, width = array1.shape
    normalized_mi = np.zeros((height, width))
    raw_mi = np.zeros((height, width))
    entropy_y = np.zeros((height, width))
    
    # Create coordinates and flatten arrays
    coords = np.array([[i, j] for i in range(height) for j in range(width)])
    
    # Normalize the input arrays
    values1_norm = normalize_data(array1)
    values2_norm = normalize_data(array2)
    
    # Build KDTree for fast neighborhood queries
    tree = KDTree(coords)
    
    for i in range(height):
        for j in range(width):
            # Find nearest neighbors (including self)
            _, neighbor_indices = tree.query(
                [[i, j]], 
                k=neighborhood_size
            )
            
            # Extract local samples (already normalized)
            X_local = values1_norm[neighbor_indices[0]]
            Y_local = values2_norm[neighbor_indices[0]]
            
            # Compute raw mutual information
            info = mi(X_local, Y_local, k=k)
            raw_mi[i, j] = info
            
            # Compute entropy of second variable (H(y))
            h_y = entropy(Y_local, k=k)
            entropy_y[i, j] = h_y
            
            # Normalize: I(x;y)/H(y)
            if h_y > 1e-10:  # Avoid division by zero
                normalized_mi[i, j] = info / h_y
            else:
                normalized_mi[i, j] = 0.0
    
    # Clip to [0, 1] since I(x;y) ≤ H(y)
    normalized_mi = np.clip(normalized_mi, 0, 1)
    
    return normalized_mi, raw_mi, entropy_y