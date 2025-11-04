"""Preprocess ORFs

Generalized script version of the `nb_aaseq_explore_<KO>.ipynb` notebooks.
See files in scripts/orf_preprocessing for example usage.

"""

import argparse
import os, sys
import numpy as np
import pandas as pd
import scipy.cluster.hierarchy as sch
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib import gridspec

from collections import Counter
import tqdm as tqdm

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.SeqRecord import SeqRecord
from Bio.Align import PairwiseAligner


AA_LIST = sorted([
    'A',  # Alanine
    'R',  # Arginine
    'N',  # Asparagine
    'D',  # Aspartic acid
    'C',  # Cysteine
    'E',  # Glutamic acid
    'Q',  # Glutamine
    'G',  # Glycine
    'H',  # Histidine
    'I',  # Isoleucine
    'L',  # Leucine
    'K',  # Lysine
    'M',  # Methionine
    'F',  # Phenylalanine
    'P',  # Proline
    'S',  # Serine
    'T',  # Threonine
    'W',  # Tryptophan
    'Y',  # Tyrosine
    'V',  # Valine
])
AA_TO_IDX = {aa: i for i, aa in enumerate(AA_LIST)}


def parse_args(args):
    parser = argparse.ArgumentParser()
    parser.add_argument("-ko", type=str, required=True)
    parser.add_argument("-i", "--infpath", type=str, required=True,
                        help="Path to ORF sequence tsv file, containing " \
                        "columns `orf` and `aa_sequence`")
    parser.add_argument("-o", "--outdir", type=str, required=True)
    parser.add_argument("-t", "--length_threshold", type=int, required=True)
    parser.add_argument("-r", "--reference_fpath", type=str, default=None,
                        help="Path to fasta file containing any reference " \
                        "sequences to include in output file " \
                        "`long_sequences.fasta`")
    parser.add_argument("--alphafold_batchsize", type=int, default=6)

    return parser.parse_args(args)


def main(args):
    # Process command line args
    KO_NUMBER = args.ko
    ORF_FPATH = args.infpath
    OUTDIR = args.outdir  # f"../out/aaseqs/{KO_NUMBER}"
    REF_SEQ_FPATH = args.reference_fpath

    LONG_SEQ_THRESH = args.length_threshold
    ALPHA_FOLD_BATCH_SIZE = args.alphafold_batchsize
    
    # Housekeeping
    os.makedirs(OUTDIR, exist_ok=True)

    IMGDIR = os.path.join(OUTDIR, "images")
    os.makedirs(IMGDIR, exist_ok=True)
    
    # Load the aa sequence data
    DF_SEQ = pd.read_csv(ORF_FPATH, delimiter="\t").set_index("orf")

    # Compute sequence lengths, and whether the sequence is complete.
    DF_SEQ["length"] = DF_SEQ["aa_sequence"].map(len)
    DF_SEQ["start_codon"] = DF_SEQ["aa_sequence"].map(lambda x: x[0] == "M")
    DF_SEQ["stop_codon"] = DF_SEQ["aa_sequence"].map(lambda x: x[-1] == "*")
    DF_SEQ["complete"] = DF_SEQ["start_codon"] & DF_SEQ["stop_codon"]

    long_seqs = DF_SEQ[DF_SEQ["length"] > LONG_SEQ_THRESH].copy()
    long_seqs["aa_sequence"] = long_seqs["aa_sequence"].str.replace(
        "*", "", regex=False
    )

    # Write the long sequences to an output file
    long_seqs_outfpath = os.path.join(OUTDIR, "long_sequences.fasta")
    with open(long_seqs_outfpath, "w") as fasta_file:
        for orf_id, row in long_seqs.iterrows():
            sequence = row["aa_sequence"]
            fasta_file.write(f">{orf_id}\n{sequence}\n")
        
        # Include reference sequence(s), if specified
        if REF_SEQ_FPATH:
            for record in SeqIO.parse(REF_SEQ_FPATH, "fasta"):
                refseqid = record.id
                refseq = record.seq
                fasta_file.write(f">{refseqid}\n{refseq}\n")

    print(
        f"Found {np.sum(DF_SEQ["complete"])}/{len(DF_SEQ)} sequences beginning "
        f"with start codon (M) and ending with stop codon (*)"
    )

    COMPLETE_DF = DF_SEQ.loc[DF_SEQ["complete"]]

    ###########################################################################
    ##  Plot lengths of all sequences, then the subset of the complete ones
    
    fig, ax = plt.subplots(1, 1)
    ax.hist(DF_SEQ["length"], bins=20)
    ax.set_xlabel(f"length")
    ax.set_ylabel(f"count")
    ax.set_title(f"All aa sequence lengths")
    plt.savefig(
        f"{IMGDIR}/all_seqs_length_dist.png", 
        transparent=False, bbox_inches=None,
    )
    plt.close()


    fig, ax = plt.subplots(1, 1)
    # fig, axes = plt.subplots(1, 2, figsize=(8,3))
    # ax = axes[0]
    ax.hist(COMPLETE_DF["length"], bins=20)
    ax.set_xlabel(f"length")
    ax.set_ylabel(f"count")
    ax.set_title(f"Complete aa sequence lengths")
    # x0, y0 = 1160, 0
    # x1, y1 = 1275, 50
    # rect = patches.Rectangle(
        # (x0, y0),           # (x, y) = lower-left corner
        # x1 - x0,            # width
        # y1 - y0,            # height
        # linewidth=0.5,
        # edgecolor='black',
        # facecolor='none'
    # )
    # ax.add_patch(rect)
    # ax = axes[1]
    # ax.hist(COMPLETE_DF.loc[COMPLETE_DF["length"] > 1100, "length"], bins=20)
    # ax.set_xlabel(f"length")
    # ax.set_ylabel(f"count")
    # fig.suptitle(f"Complete aa sequence lengths")
    plt.savefig(
        f"{IMGDIR}/complete_seqs_length_dist.png", 
        transparent=False, bbox_inches=None,
    )
    plt.close()


    ###########################################################################
    ##  Additional information using biopython
    
    # Initialize the instability and p[AA] columns
    DF_SEQ['instability'] = None
    for aa in AA_LIST:
        DF_SEQ['p' + aa] = 0.

    # Compute instability index and percent composition for each amino acid
    for idx, row in DF_SEQ.iterrows():
        seq = row['aa_sequence']
        seq = Seq(seq[0:-1] if seq[-1] == '*' else seq) # Don't include final '*'
        res = ProteinAnalysis(seq)
        aa_comp = res.amino_acids_percent
        for aa in aa_comp:
            p = aa_comp.get(aa, 0.)
            DF_SEQ.loc[idx, 'p' + aa] = p
        instability = res.instability_index()
        DF_SEQ.loc[idx, 'instability'] = instability

    COMPLETE_DF = DF_SEQ[DF_SEQ.complete]


    ##############################################################################
    ##  Plotting
    
    # Plot the composition statistics and the instabilities for all sequences
    aa_comp_matrix = DF_SEQ[["p" + aa for aa in AA_LIST]].values
    instabilities = DF_SEQ["instability"].values
    seq_lengths = DF_SEQ["length"].values
    sorted_idxs = np.argsort(seq_lengths)

    # Plot composition
    fig, ax = plt.subplots(1, 1, figsize=(5, 6))
    sc = ax.imshow(aa_comp_matrix[sorted_idxs], aspect="auto", cmap="twilight")
    ax.set_xticks(range(len(AA_LIST)), labels=AA_LIST)
    ax.set_xlabel(f"aa")
    ax.set_ylabel(f"sequence (sorted by length)")
    ax.set_title(f"aa% (all sequences)")
    fig.colorbar(sc)
    plt.savefig(
        f"{IMGDIR}/all_seqs_composition.png", 
        transparent=False, bbox_inches=None,
    )
    plt.close()

    # Pot instability
    fig, ax = plt.subplots(1, 1)
    ax.plot(seq_lengths[sorted_idxs], instabilities[sorted_idxs], '.')
    ax.set_xlabel(f"length")
    ax.set_ylabel(f"instability index")
    ax.set_title(f"Instability by sequence length (all sequences)")
    ax.axhline(40, color='r', linestyle="--", label="stability threshold")
    ax.legend()
    plt.savefig(
        f"{IMGDIR}/all_seqs_instability.png", 
        transparent=False, bbox_inches=None,
    )
    plt.close()

    # Plot the composition statistics and the instabilities for complete seqs
    aa_comp_matrix = COMPLETE_DF[["p" + aa for aa in AA_LIST]].values
    instabilities = COMPLETE_DF["instability"].values
    seq_lengths = COMPLETE_DF["length"].values
    sorted_idxs = np.argsort(seq_lengths)

    # Plot composition
    fig, ax = plt.subplots(1, 1, figsize=(5, 6))
    sc = ax.imshow(aa_comp_matrix[sorted_idxs], aspect="auto", cmap="twilight")
    ax.set_xticks(range(len(AA_LIST)), labels=AA_LIST)
    ax.set_xlabel(f"aa")
    ax.set_ylabel(f"sequence (sorted by length)")
    ax.set_title(f"aa% (complete sequences)")
    fig.colorbar(sc)
    plt.savefig(
        f"{IMGDIR}/complete_seqs_composition.png", 
        transparent=False, bbox_inches=None,
    )
    plt.close()

    # Pot instability
    fig, ax = plt.subplots(1, 1)
    ax.plot(seq_lengths[sorted_idxs], instabilities[sorted_idxs], '.')
    ax.set_xlabel(f"length")
    ax.set_ylabel(f"instability index")
    ax.set_title(f"Instability by sequence length (complete sequences)")
    ax.axhline(40, color='r', linestyle="--", label="stability threshold")
    ax.legend()
    plt.savefig(
        f"{IMGDIR}/complete_seqs_instability.png", 
        transparent=False, bbox_inches=None,
    )
    plt.close()

    ###########################################################################
    ##  Sequence alignment of long, complete sequences

    df = COMPLETE_DF.loc[COMPLETE_DF["length"] > LONG_SEQ_THRESH].copy()

    seqs = df['aa_sequence'].values
    
    # Select scoring scheme
    # see https://blast.ncbi.nlm.nih.gov/html/sub_matrix.html
    aligner = PairwiseAligner(scoring="blastp")
    print(aligner)

    alignments = [[None for _ in range(len(seqs))] for _ in range(len(seqs))]
    similarity_matrix = np.nan * np.ones([len(seqs), len(seqs)])
    for i in tqdm.trange(len(seqs)):
        seq1 = seqs[i]
        for j in range(i, len(seqs)):
            seq2 = seqs[j]
            alignment = aligner.align(seq1, seq2)
            alignments[i][j] = alignment
            score = alignment.score
            similarity_matrix[i, j] = score
            similarity_matrix[j, i] = score

    # Compute the distance matrix from the similarity matrix
    distance_matrix = np.max(similarity_matrix) - similarity_matrix

    # Plot pairwise alignment scores
    fig, ax = plt.subplots(1, 1, figsize=(6, 6))
    sc = ax.imshow(similarity_matrix, aspect="equal", cmap="viridis", origin="lower")
    ax.set_xlabel(f"sequence $i$")
    ax.set_ylabel(f"sequence $j$")
    ax.set_title(f"Pairwise alignment score")
    fig.colorbar(sc)
    plt.savefig(
        f"{IMGDIR}/long_complete_pairwise_alignments.png", 
        transparent=False, bbox_inches="tight",
    )
    plt.close()
    
    ##############################################################################
    ##  Hierarchical clustering with scipy
    ##  https://docs.scipy.org/doc/scipy/reference/cluster.hierarchy.html
    
    # Perform hierarchical clustering
    linkage = sch.linkage(
        np.triu(distance_matrix), 
        method='single',
    )
    SIM_THRESH = 0.7 * np.max(linkage[:,2])

    # Determine an ordering through a dendrogram
    dendro = sch.dendrogram(
        linkage,
        no_plot=True,
    )
    order = dendro['leaves']
    print("Dendrogram ordering:", order)

    # Cluster sequences based on linkages
    cluster_labels = sch.fcluster(
        linkage, 
        t=SIM_THRESH, 
        criterion='distance'
    )

    print("Cluster assignments:")
    print("--------------------")
    for i in range(5):
        print(f"Sequence {i} → Cluster {cluster_labels[i]}")
    print("...")

    # Reorder the similarity matrix
    ordered_similarity_matrix = similarity_matrix[np.ix_(order, order)]
    ordered_cluster_labels = cluster_labels[order]

    # Store computed cluster in dataframe
    df.loc[:,"cluster"] = cluster_labels
    df.head()

    for idx in order[0:5]:
        print(
            f"Seq {idx} (cluster {cluster_labels[idx]}):", 
            df.iloc[idx]['aa_sequence']
        )


    ###########################################################################
    ##  Plot sequence alignment with dendrogram
    fig = plt.figure(figsize=(10, 6))

    gs = gridspec.GridSpec(1, 2, width_ratios=[1, 4], wspace=0.05)
    ax0 = plt.subplot(gs[0])
    ax1 = plt.subplot(gs[1])

    # Plot dendrogram
    dendro = sch.dendrogram(
        linkage, 
        orientation="left",
        no_labels=True,
        color_threshold=SIM_THRESH,
        above_threshold_color='k',
        ax=ax0, 
    )
    assert np.allclose(order, dendro["leaves"]), "Order changed!"
    ax0.axvline(SIM_THRESH, color='k', linestyle=':')
    ax0.axis('off')

    # Plot similarity matrix
    # to_plot = ordered_similarity_matrix.copy()
    # np.fill_diagonal(to_plot, np.nan)
    to_plot = ordered_similarity_matrix
    sc = ax1.imshow(
        to_plot, 
        aspect="equal", 
        cmap="viridis", 
        origin="lower"
    )

    fig.colorbar(sc)
    ax1.set_title(f"Pairwise sequence alignment and clustering")
    ax1.set_xticks([])
    ax1.set_yticks(
        range(len(seqs)), 
        labels=[f"({c}) {i}" for c, i in zip(ordered_cluster_labels, order)], 
        fontsize=6
    )
    plt.savefig(
        f"{IMGDIR}/long_complete_dendrogram.png", 
        transparent=False, bbox_inches="tight",
    )
    plt.close()

    ###########################################################################
    ##  Save output
    
    np.savetxt(f"{OUTDIR}/long_complete_orf_ids.txt", 
            df.index.values, "%s")
    np.savetxt(f"{OUTDIR}/long_complete_aa_seqs.txt", 
            [s[0:-1] for s in df.aa_sequence.values], "%s")
    np.save(f"{OUTDIR}/long_complete_pwalignment_scores.npy", similarity_matrix)
    df.to_csv(f"{OUTDIR}/long_complete_dataframe.csv")
    np.savetxt(f"{OUTDIR}/length_threshold.txt", [LONG_SEQ_THRESH], "%d")

    # Save fasta files for sequences
    records = []
    for idx, row in df.iterrows():
        seq = row["aa_sequence"][0:-1]
        id = str(idx)
        record = SeqRecord(Seq(seq), id=id, name=id, description="")
        records.append(record)
        
    SeqIO.write(records, f"{OUTDIR}/long_complete_aa_seqs.fasta", "fasta")


    # Save single-sequence fasta files for alphafold2
    records = []
    for idx, row in df.iterrows():
        seq = row["aa_sequence"][0:-1]
        id = str(idx)
        record = SeqRecord(Seq(seq), id=id, name=id, description="")
        records.append(record)

    subdir = f"{OUTDIR}/alphafold_batches_{KO_NUMBER}_complete"
    os.makedirs(subdir, exist_ok=True)
    counter = 1
    for i, record in enumerate(records):
        if i % ALPHA_FOLD_BATCH_SIZE == 0:
            outdir = f"{subdir}/alphafold_{KO_NUMBER}_complete_batch{counter}"
            counter += 1
        os.makedirs(outdir, exist_ok=True)
        SeqIO.write(record, f"{outdir}/{record.id}.fasta", "fasta")


    ###########################################################################
    ##  Map incomplete sequences to complete ones
    ##  Take all short or incomplete sequences and map them to each of the long
    ##  complete sequences.

    target_df = COMPLETE_DF.loc[COMPLETE_DF["length"] > LONG_SEQ_THRESH]
    query_df = DF_SEQ[
        (~DF_SEQ["complete"]) | (DF_SEQ["length"] < LONG_SEQ_THRESH)
    ]

    target_seqs = target_df['aa_sequence'].to_numpy()
    query_seqs = query_df["aa_sequence"].to_numpy()

    target_orfs = target_df.index.values
    query_orfs = query_df.index.values

    ntargets = len(target_seqs)
    nqueries = len(query_seqs)

    print(f"{ntargets} target sequences (long and complete)")
    print(f"{nqueries} query sequences (short or incomplete)")

    assert len(DF_SEQ) == ntargets + nqueries, "Lengths don't add to total!"


    aligner = PairwiseAligner(scoring="blastp", mode="local")
    print(aligner)

    scores = np.zeros([nqueries, ntargets])
    for target_idx in tqdm.trange(ntargets, desc="target sequence"):
        target_seq = target_seqs[target_idx]
        for query_idx in range(nqueries):
            query_seq = query_seqs[query_idx]
            score = aligner.score(query_seq, target_seq)
            scores[query_idx, target_idx] = score

    np.save(f"{OUTDIR}/incomplete_vs_complete_alignment_scores.npy", scores)


    # Normalize scores by query length
    scores_normalized = scores / query_df.length.values[:,None]

    best_match_idxs = np.argmax(scores, axis=1)
    best_match_orfs = target_orfs[best_match_idxs]
    best_match_scores = np.zeros(nqueries)
    with open(f"{OUTDIR}/incomplete_to_complete_seq_mapping.tsv", "w") as f:
        f.write("query (incomplete)\ttarget (complete)\tscore\n", )
        for i in range(nqueries):
            q_orf = query_orfs[i]
            t_orf = best_match_orfs[i]
            score = scores_normalized[i,best_match_idxs[i]]
            best_match_scores[i] = score
            f.write(f"{q_orf}\t{t_orf}\t{score}\n")

    
    # Save top alignments
    with open(f"{OUTDIR}/incomplete_to_complete_top_alignments.txt", "w") as f:
        for i in range(nqueries):
            q_orf = query_orfs[i]
            t_orf = best_match_orfs[i]
            query_seq = query_seqs[i]
            target_seq = target_seqs[best_match_idxs[i]]
            score = scores_normalized[i,best_match_idxs[i]]
            # Realign to save alignment
            alignment = aligner.align(query_seq, target_seq)[0]
            f.write(f"query={q_orf}\ntarget={t_orf}\nscore={score}\n\n")
            f.write(str(alignment) + "\n\n")

    
    ###########################################################################
    ##  Plotting

    fig, ax = plt.subplots(1, 1)
    ax.hist(best_match_scores)
    ax.set_xlabel(f"normalized score (score / query length)")
    ax.set_ylabel(f"count")
    ax.set_title(f"Top match score distribution")
    plt.savefig(
        f"{IMGDIR}/short_seqs_top_match_distribution.png", 
        transparent=False, bbox_inches=None,
    )
    plt.close()


    fig, ax = plt.subplots(1, 1)
    counts = Counter(best_match_idxs)
    x = np.arange(ntargets)
    y = [counts[i] for i in x]
    plt.bar(x, y)
    ax.set_xlabel(f"target query idx")
    ax.set_ylabel(f"number of query matches")
    ax.set_title(f"title")
    plt.savefig(
        f"{IMGDIR}/short_seqs_num_matches.png", 
        transparent=False, bbox_inches=None,
    )
    plt.close()

    print("Done!")


if __name__ == "__main__":
    args = parse_args(sys.argv[1:])
    main(args)
