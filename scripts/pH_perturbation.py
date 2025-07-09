
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys

DATDIR = '/projects/p32818/metagenomic_data/data'

from functions import find_cluster_from_orf, get_filepath, find_orfs, samples_from_soils 

def main():
    
    soils = ['Soil3', 'Soil5', 'Soil6', 'Soil9', 'Soil11', 'Soil12', 'Soil14', 'Soil15', 'Soil16', 'Soil17']
    #soils = ['Soil5', 'Soil6', 'Soil9', 'Soil12', 'Soil15', 'Soil16', 'Soil17']
    
    cluster_IDs = pd.read_csv('../out/cluster_ids_nap.tsv', sep='\t', header=None)
    cluster_IDs = cluster_IDs.values
    cluster_IDs = [item[0] for item in cluster_IDs]
    
    metadata = pd.read_csv(f'{DATDIR}/metadata.tsv', sep = '\t')
    metadata = metadata.set_index('sample')
    
    for soil in soils: 
        print(soil)

        #narG: K00370
        #napA: K02567
        ORFs = find_orfs(get_filepath(soil, 'annotation'), 'K02567')

        #Unique_IDs is the array, defined earlier, containing the list of cluster IDs

        sample_list = samples_from_soils(soil)


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
                        row_idx = cluster_IDs.index(cluster)
                        col_idx = sample_list.index(sample_id)
                        data[row_idx, col_idx] += rel_abundance/spikein
            
        
        print(data)
        np.savetxt(f"../out/{soil}data_nap.tsv", data, delimiter = '\t', fmt = '%0.6f')
        
        
if __name__ == "__main__":
    main()
