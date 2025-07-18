from Bio import SeqIO
import pandas as pd

# Load cluster IDs
cluster_IDs = pd.read_csv('../out/orf_ids/cluster_ids_09_nap.tsv', sep='\t', header=None)
cluster_IDs = set(cluster_IDs[0].values)  # Using a set for faster lookups

orf_sequence_pairs = []

# Process the FASTA file record by record
with open("../data/raw_data/all.coassembly_proteins_1st_ClusterDB_repseq.fasta", "r") as fasta_file:
    for record in SeqIO.parse(fasta_file, "fasta"):
        if record.id in cluster_IDs:
            orf_sequence_pairs.append({'orf': record.id, 'aa_sequence': str(record.seq)})
            print(record.id, 'sequence found', record.seq)

# Save to TSV
df = pd.DataFrame(orf_sequence_pairs)
df.to_csv('../out/aaseqs/orf_sequences_09_nap.tsv', sep='\t', index=False)
