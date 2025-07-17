'''
Assembly using a new clustering, from a cluster map 
The cluster map can be assembled by build_new_map, using an output from MMSeq
Input: a KO number (eg usage: python assembly08.py --ko 'K00370')
Creates files for T0 and T9 data containing the abundances of each cluster in each soil assembly
'''


import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys
import argparse

DATDIR = '../data'

from functions import find_cluster_from_orf, get_filepath, find_orfs, samples_from_soils, perturbed_pH_sample

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--ko', required=True, help='KO identifier (e.g. K00370)')
    args = parser.parse_args()
    
    #soils = ['Soil3', 'Soil5', 'Soil6', 'Soil9', 'Soil11', 'Soil12', 'Soil14', 'Soil15', 'Soil16', 'Soil17']
    soils = ['Soil3', 'Soil5', 'Soil9', 'Soil11', 'Soil12', 'Soil14', 'Soil15', 'Soil16', 'Soil17'] #no soil6 for nap?
    
    cluster_IDs = pd.read_csv('../out/cluster_ids_08_nap.tsv', sep='\t', header=None)
    cluster_IDs = cluster_IDs.values
    cluster_IDs = [item[0] for item in cluster_IDs]
    
    metadata = pd.read_csv(f'{DATDIR}/metadata.tsv', sep = '\t')
    metadata = metadata.set_index('sample')
    
    samples = pd.read_csv(f'{DATDIR}/T0_sampleIDs.tsv', header=None)[0]
    sample_list = list(samples)
    
    def find_cluster(orf):
        prefix = orf.split('.')[0]
        df = pd.read_csv('../out/cluster08map_nap.tsv', sep ='\t', header=None)
        cluster = df[df[0] == orf][2].tolist()
        return cluster[0]
    
    

    #narG: K00370
    #napA: K02567
    ORFs = find_orfs(f"../data/subset_{args.ko}/T0.coassembly_annotations_{args.ko}.tsv", args.ko)


    data = np.zeros((len(cluster_IDs), 20)) #data for plot stored here, each row is a dinstinct cluster

    abundance_data = pd.read_csv(f"../data/subset_{args.ko}/T0_all_samples_{args.ko}.bed", sep='\s+', header=None)
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
    np.savetxt(f"../out/T0data_08_{args.ko}.tsv", data, delimiter = '\t', fmt = '%0.6f')
    
    def samples(soil):
        sample_IDs = []
        metadata = pd.read_csv('../data/metadata.tsv', sep='\t')
        metadata = metadata.set_index('sample')
        soil_ = soil + '_'
        pHs = []
        for sample in metadata.index:
            if 'None' in sample and 'T9' in sample and soil_ in sample and 'Nitrate' not in sample:
                sample_IDs.append(sample)
                pHs.append(perturbed_pH_sample(sample))
                
        #make sure the samples are organized by increaing perturbed pH
        pHs = np.array(pHs)
        sample_IDs = np.array(sample_IDs)
        
        indices = np.argsort(pHs)
        
        sample_IDs = sample_IDs[indices]
        sample_IDs = sample_IDs.tolist()
        
        
        return sample_IDs
    
    
    for soil in soils: 
        print(soil)

        #narG: K00370
        #napA: K02567
        ORFs = find_orfs(f"../data/subset_{args.ko}/{soil}.coassembly_annotations_{args.ko}.tsv", args.ko)

        #Unique_IDs is the array, defined earlier, containing the list of cluster IDs

        sample_list = samples(soil)


        data = np.zeros((len(cluster_IDs), 11)) #data for plot stored here, each row is a dinstinct cluster

        abundance_data = pd.read_csv(f"../data/subset_{args.ko}/{soil}_all_samples_{args.ko}.bed", sep='\s+', header=None)
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
        np.savetxt(f"../out/{soil}data_08_{args.ko}.tsv", data, delimiter = '\t', fmt = '%0.6f')
        
        
if __name__ == "__main__":
    main()
