import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys

DATDIR = '/projects/p32818/metagenomic_data/data'

from functions import find_cluster_from_orf, get_filepath, find_orfs

def main():
    soils = ['Soil3', 'Soil5', 'Soil6', 'Soil9', 'Soil11', 'Soil12', 'Soil14', 'Soil15', 'Soil16', 'Soil17']

    cluster_IDs = pd.read_csv('out/cluster_ids.tsv', sep='\t', header=None)
    cluster_IDs = cluster_IDs.values
    cluster_IDs = [item[0] for item in cluster_IDs]

    cluster_tracker = np.zeros((len(cluster_IDs), len(soils)))


    for soil in soils:
        yes = 0
        no = 0
        col_idx = soils.index(soil)
        ORFS = find_orfs(get_filepath(soil, 'annotation'), 'K00370')
        for orf in ORFS:
            cluster = find_cluster_from_orf(orf)
            if cluster in cluster_IDs:
                row_idx = cluster_IDs.index(cluster)
                yes += 1
                print('yes')
                cluster_tracker[row_idx, col_idx] = 1
            else:
                no += 1
                print('no')
        print(soil, 'yes = ', yes, 'no = ', no)
        
    print(cluster_tracker)
    
    np.savetxt("out/track_clusters.tsv", cluster_tracker, delimiter = '\t', fmt = '%0.6f')
    
    
if __name__ == "__main__":
    main()