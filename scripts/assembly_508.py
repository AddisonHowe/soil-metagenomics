
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys

DATDIR = '../data'

from functions import find_cluster_from_orf, get_filepath, find_orfs, samples_from_soils 

def main():
    
    soils = ['Soil3', 'Soil5', 'Soil6', 'Soil9', 'Soil11', 'Soil12', 'Soil14', 'Soil15', 'Soil16', 'Soil17']
    
    cluster_IDs = pd.read_csv('../out/complete_ids.tsv', sep='\t', header=None)
    cluster_IDs = cluster_IDs.values
    cluster_IDs = [item[0] for item in cluster_IDs]
    
    metadata = pd.read_csv(f'{DATDIR}/metadata.tsv', sep = '\t')
    metadata = metadata.set_index('sample')
    
    def samples(soil):
        sample_IDs = []
        metadata = pd.read_csv('../data/metadata.tsv', sep='\t')
        metadata = metadata.set_index('sample')
        soil_ = soil + '_'
        for sample in metadata.index:
            if 'None' in sample and 'T9' in sample and soil_ in sample and 'Nitrate' not in sample:
                sample_IDs.append(sample)
        
        return sample_IDs
    
    def find_cluster(orf):
        prefix = orf.split('.')[0]
        FPATH = f"../data/subset_K00370/{soil}.coassembly_annotations_K00370.tsv"
        df = pd.read_csv(FPATH, sep='\t', header=None)
        cluster = df[df[0] == orf][2].tolist()
        return cluster[0]
    
    for soil in soils: 
        print(soil)

        #narG: K00370
        #napA: K02567
        ORFs = find_orfs(f"../data/subset_K00370/{soil}.coassembly_annotations_K00370.tsv", 'K00370')

        #Unique_IDs is the array, defined earlier, containing the list of cluster IDs

        sample_list = samples(soil)


        data = np.zeros((len(cluster_IDs), 11)) #data for plot stored here, each row is a dinstinct cluster

        abundance_data = pd.read_csv(f"../data/subset_K00370/{soil}_all_samples_K00370.bed", sep='\s+', header=None)
        filtered_data = abundance_data[abundance_data.iloc[:, 1].isin(ORFs)]
        for i in range(len(filtered_data)):
            sample_id = filtered_data.iloc[i, 0]
            if sample_id in sample_list:
                orf = filtered_data.iloc[i, 1]
                rel_abundance = filtered_data.iloc[i, 2]
                spikein = metadata.loc[sample_id, 'spikein_sum']
                cluster = find_cluster(orf)
                if cluster in cluster_IDs:
                    row_idx = cluster_IDs.index(cluster)
                    col_idx = sample_list.index(sample_id)
                    data[row_idx, col_idx] += rel_abundance/spikein
            
        
        print(data)
        np.savetxt(f"../out/{soil}data_508_nar.tsv", data, delimiter = '\t', fmt = '%0.6f')
        
        
if __name__ == "__main__":
    main()
