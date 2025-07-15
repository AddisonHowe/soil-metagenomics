from Bio import SeqIO
import pandas as pd

cluster_IDs = pd.read_csv('../out/cluster_ids.tsv', sep='\t', header=None)
cluster_IDs = set(cluster_IDs[0].values)  # Faster lookups

# Create a disk-based index (memory-efficient)
seq_index = SeqIO.index("../data/raw_data/all.coassembly_proteins_1st_ClusterDB_repseq.fasta", "fasta")

orf_sequence_pairs = []
for orf in cluster_IDs:
    if orf in seq_index:
        sequence = str(seq_index[orf].seq)
        orf_sequence_pairs.append({'orf': orf, 'aa_sequence': sequence})
        print(orf, 'found', sequence)

seq_index.close()  # Important to free resources

# Save to TSV
df = pd.DataFrame(orf_sequence_pairs)
df.to_csv('../out/orf_sequences2.tsv', sep='\t', index=False)