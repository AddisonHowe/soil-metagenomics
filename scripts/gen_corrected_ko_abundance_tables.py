"""Generate corrected KO abundance tables, handling spike-ins


Produces the following csv files in data/KO_tables_corrected:
    chl_neg_samples_corrected.csv
    chl_neg_samples_corrected_normalized.csv
    chl_pos_samples_corrected.csv
    chl_pos_samples_corrected_normalized.csv
    no_nitrate_corrected.csv
    no_nitrate_corrected_normalized.csv
    t0_samples_corrected.csv
    t0_samples_corrected_normalized.csv

These files show the sum of accepted and rejected reads, with the spike-in 
strains not included, i.e.
    corrected = accepted + rejected - spikein1 - spikein2
The normalized version of each file then divides this number by the total 
spike-in sum. Note that we have not adjusted the spike-in sums in any way.

The input files for this script are the accepted and rejected tables produced
in the KO-wrangling pipeline, and which are found in data/tables.

T0 abundance levels reflect the values that are found in the Soil batches, not
the T0 batch. That is, the values here should match those found in the KO tables
located in data/KO_table/Soil<N>.KO_absolute_average_depth.tsv.

"""

import os
import pandas as pd

tables_dir = "data/tables"
spikeins_dir = "data/KO_table"
metadata_fpath = "data/metadata.tsv"
OUTDIR = "data/KO_tables_corrected"


os.makedirs(OUTDIR, exist_ok=True)

outdir_raw = OUTDIR
outdir_norm = OUTDIR

os.makedirs(outdir_raw, exist_ok=True)
os.makedirs(outdir_norm, exist_ok=True)


TAXID_SPIKE1 = 1
TAXID_SPIKE2 = 2849180

SAMPLE_SET_KEYS = [
    "chl_neg_samples", "chl_pos_samples", "no_nitrate", "t0_samples",
]


# Load the sample metadata file containing spikein sums
df_metadata = pd.read_csv(
    metadata_fpath, sep="\t", index_col=0,
)
soil_spikein_sums = {}
for idx, row in df_metadata.iterrows():
    soil_spikein_sums[idx] = float(row["spikein_sum"])


for i, sample_key in enumerate(SAMPLE_SET_KEYS):

    acc_tab_fpath = f"{tables_dir}/{sample_key}_accepted.csv"
    rej_tab_fpath = f"{tables_dir}/{sample_key}_rejected.csv"
    spike1_fpath = f"{tables_dir}/by_taxa/{sample_key}_{TAXID_SPIKE1}_rejected.csv"
    spike2_fpath = f"{tables_dir}/by_taxa/{sample_key}_{TAXID_SPIKE2}_rejected.csv"


    df_acc = pd.read_csv(
        acc_tab_fpath, sep=",", index_col=0,
    ).fillna(0)

    df_rej = pd.read_csv(
        rej_tab_fpath, sep=",", index_col=0,
    ).fillna(0)

    df_spike1 = pd.read_csv(
        spike1_fpath, sep=",", index_col=0,
    ).fillna(0)

    df_spike2 = pd.read_csv(
        spike2_fpath, sep=",", index_col=0,
    ).fillna(0)

    df_corrected = df_acc + df_rej - df_spike1 - df_spike2

    df_corrected.to_csv(
        f"{outdir_raw}/{sample_key}_corrected.csv", 
        float_format="%.6f"
    )

    common = df_corrected.index.intersection(df_metadata.index)
    df_corrected = df_corrected.loc[common]
    df_meta = df_metadata.loc[common]

    df_corrected_normalized = df_corrected.div(df_meta["spikein_sum"], axis=0)
    df_corrected_normalized.to_csv(
        f"{outdir_norm}/{sample_key}_corrected_normalized.csv", 
        float_format="%.6f"
    )
