
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys
from collections import defaultdict

DATDIR = '/projects/p32818/metagenomic_data/data'



import numpy as np
import pandas as pd
from functions import find_orfs, get_filepath

soils = ['Soil3', 'Soil5', 'Soil6', 'Soil9', 'Soil11', 'Soil12', 'Soil14', 'Soil15', 'Soil16', 'Soil17']


#K00370 narG
#K02567 napP

def main():
    
    prefixes = ['T0', 'Soil3', 'Soil5', 'Soil6', 'Soil9', 'Soil11', 'Soil12', 'Soil14', 'Soil15', 'Soil16', 'Soil17']
    file_path_list = []
    for prefix in prefixes:
        file_path_list.append(get_filepath(prefix, 'annotation', 'K02567'))

    #Want to produce "map" with 1336 rows corresponding to K00370 ORFs
    
    map = []
    
    #The 2nd Column is the round 1 cluster ORFs
    
    #The 3rd Column is the round 2 clustered ORFs, with the new clustering 
    
    #produce a list of round 1 ORFs
    def find_1st_orf(orf, file_path):
        df = pd.read_csv(file_path, sep = '\t', header = None)
        return df[df[0] == orf][1].iloc[0]

    orf_1_list = []
    for i in range(len(prefixes)):
        ORF_list = find_orfs(file_path_list[i], 'K02567')
        for orf in ORF_list:
            orf_1 = find_1st_orf(orf, file_path_list[i])
            map.append([orf, orf_1, 'blank'])
            if orf_1 not in orf_1_list:
                orf_1_list.append(orf_1)
                
    print('orf_1_list: ', orf_1_list)
    print('len: ', len(orf_1_list))
            
    print('Completed 1st Round')

    def build_targeted_lookup(file_path, target_orfs):
        """Memory-efficient lookup builder"""
        lookup_dict = {}
        target_orfs = set(target_orfs)
        
        # Use low-memory chunking
        for chunk in pd.read_csv(
            file_path,
            sep='\t',
            header=None,
            usecols=[0, 1],
            dtype=str,
            chunksize=10_000_000
        ):
            # Vectorized filtering
            mask = chunk[1].isin(target_orfs)
            filtered = chunk[mask]
            lookup_dict.update(zip(filtered[1], filtered[0]))
        
        return lookup_dict

    # Usage:
    file_path = '../data/raw_data/all.coassembly_proteins_1st_ClusterDB_repseq_2ndClusterDB.tsv'
    lookup_dict = build_targeted_lookup(file_path, orf_1_list)  # Do this ONCE
    
    print('Dictionary Built')

    def find_2nd_orf(id, lookup_dict):
        """O(1) lookup from preloaded dictionary."""
        return lookup_dict.get(id, None)
    
    for entry in map:
        entry[2] = find_2nd_orf(entry[1], lookup_dict)
        print(entry[2])
    
    np.savetxt("..out/cluster08map_nap.tsv", map, delimiter = '\t', fmt = '%s')
        
if __name__ == "__main__":
    main()
