"""

"""

import argparse
import os
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import csv
import tqdm as tqdm
import time


_SOIL_INTS = [11, 12, 14, 15, 16, 17, 3, 5, 6, 9]

parser = argparse.ArgumentParser()
parser.add_argument("-s", "--soil", type=str, required=True, 
                    choices=[f"Soil{i}" for i in _SOIL_INTS])
parser.add_argument("--seqdir", type=str, required=True)
parser.add_argument("-o", "--outdir", type=str, required=True)
parser.add_argument("--pbar", action="store_true")
args = parser.parse_args()



SOIL_KEY = args.soil
DISABLE_PBAR = not args.pbar

SEQDIR = args.seqdir  # "/Users/addisonhowe/mnts/raw_data/sequencing"
DATDIR = f"{SEQDIR}/coverage/arrays"  # location of .npz soil coverage arrays
sam_header_fname = "header.sam"

OUTDIR = args.outdir # "../out/coverage"
os.makedirs(OUTDIR, exist_ok=True)




# Need to fix this name since in the basecov file there is a comma.
FIX_NAMES = {
    "NC_010943.1_Stenotrophomonas_maltophilia_K279a__strain_K279a":
        "NC_010943.1_Stenotrophomonas_maltophilia_K279a,_strain_K279a",
}




RAWNAME_TO_SAMPLE = {}
SAMPLE_TO_RAWNAME = {}
with open(f"{SEQDIR}/raw_data_name_mapping.tsv", "r") as f:
    csvreader = csv.reader(f, delimiter="\t")
    for row in csvreader:  # process each row
        RAWNAME_TO_SAMPLE[row[0]] = row[1]
        SAMPLE_TO_RAWNAME[row[1]] = row[0]




REGION_MAPPING = {}  # Map name to gene to interval

for d in os.listdir("../igv/out/search_results"):
    search_file = f"../igv/out/search_results/{d}"
    with open(search_file, "r") as f:
        for line in f:
            # Skip comment lines (optional)
            if line.startswith("#") or not line.strip():
                continue
            parts = line.strip().split("\t")
            if len(parts) >= 5:
                seqname = parts[0]
                start = int(parts[3])
                end = int(parts[4])
                if seqname not in REGION_MAPPING:
                    REGION_MAPPING[seqname] = {}
                if d not in REGION_MAPPING[seqname]:
                    REGION_MAPPING[seqname][d] = []
                REGION_MAPPING[seqname][d].append((start, end))

# for k in REGION_MAPPING:
#     print(k, REGION_MAPPING[k])




ALL_SOIL_SAMPLES = [f for f in os.listdir(DATDIR) if f.endswith(".npz")]

def order_samples(sample_list, searchkey):
    if searchkey:
        subset = [s for s in sample_list if searchkey + "_" in s]
    else:
        subset = sample_list

    t0_set = [s for s in subset if "T0" in s]
    t9_set = [s for s in subset if "T9" in s]
    no_nitrate_t0_set = [s for s in t0_set if "No_Nitrate" in s]
    no_nitrate_t9_set = [s for s in t9_set if "No_Nitrate" in s]
    t0_set = [s for s in t0_set if "No_Nitrate" not in s]
    t9_set = [s for s in t9_set if "No_Nitrate" not in s]
    chl_set = [s for s in t9_set if "CHL" in s]
    none_set = [s for s in t9_set if "None" in s]
    
    def sortfunc(s):
        sp = s.split("_")
        if sp[5][0] == "N":
            v = float(sp[6])
        else:
            v = float(sp[5])
        return v

    chl_set = sorted(chl_set, key=sortfunc)
    none_set = sorted(none_set, key=sortfunc)
    no_nitrate_t0_set = sorted(no_nitrate_t0_set, key=sortfunc)
    no_nitrate_t9_set = sorted(no_nitrate_t9_set, key=sortfunc)
    ord = t0_set + chl_set + none_set + no_nitrate_t0_set + no_nitrate_t9_set
    assert set(ord) == set(subset) and len(ord) == len(subset), \
        f"Messed up! Expected {subset}. Got {ord}"
    return ord


SAMPLES_BY_SOIL = {}
for i in [11, 12, 14, 15, 16, 17, 3, 5, 6, 9]:
    k = f"Soil{i}"
    SAMPLES_BY_SOIL[k.lower()] =  order_samples(ALL_SOIL_SAMPLES, k)




sample_list = SAMPLES_BY_SOIL[SOIL_KEY.lower()]

# for s in sample_list:
#     print(s)




def find_spiked_regions(x, dx=1):
    """Identify spiked regions."""
    if len(x) == 0:
        return []
    regions = []
    r0 = x[0]
    r1 = x[0] + 1
    idx_prev = x[0]
    for i in range(1, len(x)):
        idx = x[i]
        # Check if the current value is an extension of the current region
        if idx - idx_prev <= dx:
            # Extend the region
            r1 = idx + 1
        else:
            # Reached end of regions. Store and reset.
            regions.append([r0, r1])
            r0 = idx
            r1 = idx + 1
        idx_prev = idx
    # Append final region
    regions.append([r0, r1])
    return np.array(regions)

def get_spike_height(x, regions):
    spike_heights = np.zeros(len(regions))
    for i, region in enumerate(regions):
        r0, r1 = region
        if np.isnan(r1):
            r1 = r0
        heights = x[int(r0):int(r1 + 1)]
        spike_heights[i] = np.max(heights)
    return spike_heights




##############################################################################
##  Tests for find_spiked_regions

x = np.concatenate([
    np.zeros(10),
    np.zeros(10),
    np.arange(6),
    np.flip(np.arange(5)),
    np.zeros(10),
    np.zeros(10),
    np.zeros(10),
    np.arange(4)
], dtype=float)

kappa = 5
dx = 1
# print(len(x))
med = np.median(x)
# print(f"kappa={kappa}, med={med}")
spike_locations = np.argwhere(x > kappa * med).flatten()
# print("Spike locations:", spike_locations)
regs = find_spiked_regions(spike_locations, dx=dx)
# print("spike regions:\n", regs)
heights = get_spike_height(x, regs)
# print("spike heights:\n", heights)




dx = 200  # Padding for regions used in `find_spiked_regions`
bin_target = 500
logscale=lambda x: np.where(x == 0, np.nan, np.log2(x, where=(x != 0)))

spike_marker_size = None

for sampidx, sample in tqdm.tqdm(
        enumerate(sample_list),
        total=len(sample_list), 
        desc="Sample",
        leave=False,
        disable=DISABLE_PBAR,
):
    arr = np.load(f"{DATDIR}/{sample}", allow_pickle=True)
    scaffold_keys = list(arr.keys())
    for scaffold_key in tqdm.tqdm(
            scaffold_keys, 
            total=len(scaffold_keys), 
            desc="scaffold", 
            leave=False,
            disable=DISABLE_PBAR,
    ):
        x = arr[scaffold_key]
        outdir = f"{OUTDIR}/by_scaffold/{scaffold_key}/arrays"
        imgdir = f"{OUTDIR}/by_scaffold/{scaffold_key}/images"
        os.makedirs(outdir, exist_ok=True)
        os.makedirs(imgdir, exist_ok=True)

        #####################################################################
        ##  Compute and save spike information, etc.

        bin_frac = bin_target / len(x)
        xpos = x[x > 0]
        med = np.median(xpos)  # store the median value
        meanlog = np.mean(np.log2(xpos))  # store the mean log2 of positives
        stdlog = np.std(np.log2(xpos))  # store the std log2 of positives
        p25 = np.percentile(xpos, 25)  # store the 25th percentile
        p75 = np.percentile(xpos, 75)  # store the 75th percentile
        iqr = p75 - p25  # compute the IQR
        # threshold = med + 1.5 * iqr
        threshold = 2**(meanlog + 2.5 * stdlog)
        spike_locations = np.argwhere(x > threshold).flatten()
        spike_regions = find_spiked_regions(spike_locations, dx=dx)
        spike_heights = get_spike_height(x, spike_regions)
        np.save(
            f"{outdir}/{sample.replace(".npz", "")}_spike_regions.npy", 
            spike_regions
        )
        np.save(
            f"{outdir}/{sample.replace(".npz", "")}_spike_heights.npy", 
            spike_heights
        )

        # Compute bins
        n = len(x)
        binsize = int(np.ceil(bin_frac * n))
        nbins = n // binsize + (n % binsize != 0)
        remainder = n % binsize
        if remainder > 0:
            x = np.concatenate([x, np.zeros(binsize - remainder)])
        xbinned = x.reshape([-1, binsize]).max(axis=1)

        #####################################################################
        ##  Plot data

        fig, ax = plt.subplots(1, 1, figsize=(10,6))

        # Plot base coverage
        ts = 1 + np.arange(0, len(xbinned))
        ax.plot(
            ts, logscale(xbinned), 
            linestyle='None',
            marker=".", 
            markersize=spike_marker_size,
        )

        ax.set_xlim(-ts.max() * 0.01, ts.max() * 1.05)
        # ax.set_ylim(-0.01, ax.get_ylim()[1])

        # Mark median and threshold
        xlims = ax.get_xlim()
        # ax.hlines(
        #     logscale(med), *xlims, 
        #     linestyle=':', linewidth=1, color="r",
        #     label=f"med",
        # )
        # ax.hlines(
        #     logscale(p75), *xlims, 
        #     linestyle="--", linewidth=1, color="r",
        #     label="75th pctl",
        # )
        ax.hlines(
            logscale(threshold), *xlims, 
            linestyle="-", linewidth=1, color="r",
            label="threshold",
        )
        ax.set_xlim(*xlims)

        # Mark spike regions
        ylims = ax.get_ylim()
        for reg, height in zip(spike_regions, spike_heights):
            if height >= threshold:
                ax.axvspan(
                    reg[0], reg[1],*ylims, 
                    color='orange', alpha=0.5,
                )
        ax.set_ylim(*ylims)
        
        # Add annotations
        annotations = REGION_MAPPING.get(scaffold_key, {})
        ylims = ax.get_ylim()
        labelheight = ylims[1]
        for gene in annotations:
            prev_start, prev_end = -1, -1
            regions = annotations[gene]
            for region in regions:
                start, end = region
                if prev_start != start or prev_end != end:
                    ax.axvline(start, *ylims, linestyle=":", color="k", zorder=1)
                    ax.axvline(end, *ylims, linestyle=":", color="k", zorder=1)
                    tbox = ax.text(
                        start, labelheight, gene.split("_")[1],
                        verticalalignment="top",
                    )
                    labelheight -= 0.025 * (ylims[1] - ylims[0])
                prev_start, prev_end = start, end
        
        ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        ax.set_xlabel(f"position")
        ax.set_ylabel("log2 read count")
        title = scaffold_key
        subtitle = sample.removeprefix("coverage_arrays_").removesuffix(".npz")
        ax.set_title(title + "\n" + subtitle)

        plt.savefig(f"{imgdir}/{sample.replace(".npz", "")}.png", 
                    bbox_inches="tight")
        plt.close()

        #####################################################################
        ##  Plot binned data
        
        fig, ax = plt.subplots(1, 1, figsize=(10,6))
        
        # Plot base coverage
        ts = 1 + np.arange(0, len(xbinned))
        ax.plot(
            ts, logscale(xbinned), 
            linestyle='None',
            marker=".", 
            markersize=spike_marker_size,
        )

        ax.set_xlim(-ts.max() * 0.01, ts.max() * 1.05)
        # ax.set_ylim(-0.01, ax.get_ylim()[1])

        # Mark median and threshold
        xlims = ax.get_xlim()
        # ax.hlines(
        #     logscale(med), *xlims, 
        #     linestyle=':', linewidth=1, color="r",
        #     label="med",
        # )
        # ax.hlines(
        #     logscale(p75), *xlims, 
        #     linestyle="--", linewidth=1, color="r",
        #     label="75th pctl",
        # )
        ax.hlines(
            logscale(threshold), *xlims, 
            linestyle="-", linewidth=1, color="r",
            label="threshold",
        )
        ax.set_xlim(*xlims)

        # Mark spike regions
        ylims = ax.get_ylim()
        for reg, height in zip(spike_regions, spike_heights):
            if height >= threshold:
                r0 = reg[0] // binsize
                r1 = reg[1] // binsize
                ax.axvspan(
                    r0, r1, *ylims, 
                    color='orange', alpha=0.5,
                )
        ax.set_ylim(*ylims)
        
        # Add annotations
        annotations = REGION_MAPPING.get(scaffold_key, {})
        ylims = ax.get_ylim()
        labelheight = ylims[1]
        for gene in annotations:
            prev_start, prev_end = -1, -1
            regions = annotations[gene]
            for region in regions:
                start, end = region
                start = 1 + start // binsize
                end = 1 + end // binsize
                if prev_start != start or prev_end != end:
                    ax.axvline(start, *ylims, linestyle=":", color="k", zorder=1)
                    ax.axvline(end, *ylims, linestyle=":", color="k", zorder=1)
                    tbox = ax.text(
                        start, labelheight, gene.split("_")[1],
                        verticalalignment="top",
                    )
                    labelheight -= 0.025 * (ylims[1] - ylims[0])
                prev_start, prev_end = start, end
        
        ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        ax.set_xlabel(f"bin (size {bin_frac*100:.2g}%={binsize}bp)")
        ax.set_ylabel("log2 read count")
        title = scaffold_key
        subtitle = sample.removeprefix("coverage_arrays_").removesuffix(".npz")
        ax.set_title(title + "\n" + subtitle)

        plt.savefig(f"{imgdir}/{sample.replace(".npz", "")}_binned.png", 
                    bbox_inches="tight")
        plt.close()


# #### Heatmaps



SCAFFOLD_KEYS = [
    "Ralstonia_solanacearum_strain_KACC_10722",
    "CP013959.1_Staphylococcus_aureus_strain_V605",
    "NC_003902.1_Xanthomonas_campestris_pv._campestris_str._ATCC_33913",
    "NC_003295.1_Ralstonia_solanacearum_GMI1000",
    "NC_012660.1_Pseudomonas_fluorescens_SBW25",
    "NC_002570.2_Bacillus_halodurans_C-125",
    "NC_002976.3_Staphylococcus_epidermidis_RP62A",
    "NC_007492.2_Pseudomonas_fluorescens_Pf0-1",
    "NC_008782.1_Acidovorax_sp._JS42",
    "NC_010002.1_Delftia_acidovorans_SPH-1",
    "NC_011071.1_Stenotrophomonas_maltophilia_R551-3",
    "NC_008702.1_Azoarcus_sp._BH72",
    "NC_010943.1_Stenotrophomonas_maltophilia_K279a,_strain_K279a",
    "NC_013446.2_Comamonas_testosteroni_CNB-2",
    "NC_014323.1_Herbaspirillum_seropedicae_SmR1",
    "NC_015563.1_Delftia_sp._Cs1-4",
    "NC_016830.1_Pseudomonas_fluorescens_F113",
    "NZ_AJVM01000001.1_Rhizobium_sp._AP16",
    "NC_018708.1_Acidovorax_sp._KKS102",
    "NC_020561.1_Sphingomonas_sp._MM-1",
    "NZ_KB906253.1_Curvibacter_lanceolatus_ATCC_14669",
    "NZ_KB976985.1_Acinetobacter_calcoaceticus_ANC_3811",
    "NZ_KI421535.1_Burkholderia_cepacia_ATCC_25416",
    "NZ_HG916852.1_Rhizobium_sp._LPU83",
    "NZ_CP007638.1_Pseudomonas_sp._WCS374",
    "NZ_CP008788.1_Klebsiella_oxytoca_KONIH1",
    "NZ_CP010536.1_Cupriavidus_basilensis_strain_4G11",
    "NZ_JXYQ01000001.1_Acidovorax_temperans_strain_KY4",
    "NZ_CP010516.1_Cupriavidus_gilardii_CR3",
    "NZ_CP010350.1_Acinetobacter_johnsonii_XBB1",
    "NZ_CP014540.1_Acinetobacter_baumannii_strain_XH857",
    "NC_016845.1_Klebsiella_pneumoniae_subsp._pneumoniae_HS11286",
    "NC_007606.1_Shigella_dysenteriae_Sd197",
    "NC_004337.2_Shigella_flexneri_2a_str._301",
    "Lee_A8Q_1_Ecoli_contig_1_polypolish",
]




dx = 200  # Padding for regions used in `find_spiked_regions`
bin_target = 500

CMAP = "plasma"
spike_color = "cyan"
spike_marker_size = 3

imgdir = f"{OUTDIR}/images/{SOIL_KEY.lower()}"
os.makedirs(imgdir, exist_ok=True)

progbar = tqdm.tqdm(
    SCAFFOLD_KEYS, total=len(SCAFFOLD_KEYS), 
    desc="scaffold", 
    leave=False,
    disable=DISABLE_PBAR,
)
for scaffold_key in progbar:
    progbar.set_description(scaffold_key)
    
    all_x = []
    all_medians = []
    all_meanlogs = []
    all_stdlogs = []
    all_p75s = []
    all_p25s = []
    all_iqrs = []
    all_thresholds = []
    all_spike_locations = []
    all_spike_regions = []
    all_spike_heights = []
    all_xbinned = []

    for sample in tqdm.tqdm(
        sample_list, 
        total=len(sample_list), 
        disable=DISABLE_PBAR,
        desc="Sample",
        leave=False,
    ):
        arr = np.load(f"{DATDIR}/{sample}", allow_pickle=True)
        x = arr[scaffold_key]

        #####################################################################
        ##  Compute and store spike information, etc.
        xpos = x[x > 0]
        med = np.median(xpos)  # store the median value
        p75 = np.percentile(xpos, 75)  # store the 75th percentile
        p25 = np.percentile(xpos, 25)  # store the 25th percentile
        meanlog = np.mean(np.log2(xpos))  # store the mean log2 of positives
        stdlog = np.std(np.log2(xpos))  # store the std log2 of positives
        iqr = p75 - p25
        # threshold = med + 1.5 * iqr
        threshold = 2**(meanlog + 2.5 * stdlog)
        spike_locations = np.argwhere(x > threshold).flatten()
        spike_regions = find_spiked_regions(spike_locations, dx=dx)
        spike_heights = get_spike_height(x, spike_regions)

        all_x.append(x)
        all_medians.append(med)
        all_p75s.append(p75)
        all_p25s.append(p25)
        all_meanlogs.append(meanlog)
        all_stdlogs.append(stdlog)
        all_iqrs.append(iqr)
        all_thresholds.append(threshold)
        all_spike_locations.append(spike_locations)
        all_spike_regions.append(spike_regions)
        all_spike_heights.append(spike_heights)
        
        n = len(x)
        bin_frac = bin_target / n
        binsize = int(np.ceil(bin_frac * n))
        nbins = n // binsize + (n % binsize != 0)
        remainder = n % binsize
        if remainder > 0:
            x = np.concatenate([x, np.zeros(binsize - remainder)])
        xbinned = x.reshape([-1, binsize]).max(axis=1)

        all_xbinned.append(xbinned)

    all_x = np.array(all_x)
    all_medians = np.array(all_medians)
    all_thresholds = np.array(all_thresholds)
    all_p75s = np.array(all_p75s)
    all_xbinned = np.array(all_xbinned)

    #####################################################################
    ##  Plot

    fig, ax = plt.subplots(1, 1, figsize=(13,7))
    
    sc = ax.imshow(
        all_xbinned, 
        aspect='auto', 
        norm="log",
        cmap=CMAP, 
    )

    divider = make_axes_locatable(ax)
    cax = divider.append_axes('right', size='5%', pad=0.1)
    cbar = fig.colorbar(sc, cax=cax)
    cbar.ax.set_title("reads", size=10)

    ticks = np.arange(len(sample_list))
    ticklabels = [str(i+1) for i in ticks]
    annot_space = 1/3
    ax.set_ylim(1/(1-annot_space) * ax.get_ylim()[0], ax.get_ylim()[1])
    ax.set_yticks(ticks, labels=ticklabels)

    # Add annotations
    annotations = REGION_MAPPING.get(scaffold_key, {})
    gene_order = sorted(
        list(annotations.keys()), key=lambda g: annotations[g][0][0], 
    )
    ylims = ax.get_ylim()
    labelheight0 = ylims[1] - (1 - annot_space) * (ylims[1] - ylims[0]) + 0.1
    labelheight = labelheight0
    for gene in gene_order:
        prev_start, prev_end = -1, -1
        regions = annotations[gene]
        for region in regions:
            start, end = region
            start = 1 + start // binsize
            end = 1 + end // binsize
            if prev_start != start or prev_end != end:
                ax.vlines(
                    start, ylims[0], (1 - annot_space) * ylims[0], 
                    linestyle="--", color="k", alpha=1, zorder=1,
                )
                tbox = ax.text(
                    start, labelheight, gene.split("_")[1],
                    verticalalignment="top", color='r', zorder=2,
                )
                labelheight -= 0.025 * (ylims[1] - ylims[0])
                if labelheight - 0.025 * (ylims[1] - ylims[0]) >= ylims[0]:
                    labelheight = labelheight0
            prev_start, prev_end = start, end
    
    # Add spike markers
    for sampidx, spike_regions in enumerate(all_spike_regions):
        spike_heights = all_spike_heights[sampidx]
        threshold = all_thresholds[sampidx]
        for spike_region, spike_height in zip(spike_regions, spike_heights):
            if spike_height < threshold:
                continue
            pos = np.mean(spike_region)
            binpos = pos // binsize
            ax.plot(
                binpos, sampidx, 
                marker='.',
                color=spike_color,
                markersize=spike_marker_size,
            )

    soilnames = [
        s.removeprefix("coverage_arrays_").removesuffix(".npz").replace("_", " ") 
        for s in sample_list
    ]
    ax.set_yticklabels(soilnames, fontsize=8, rotation=45)
    ax.set_xlabel(f"bin (size {bin_frac*100:.2g}%={binsize}bp)")
    ax.set_ylabel("sample")
    title = scaffold_key
    subtitle = SOIL_KEY
    ax.set_title(title + "\n" + subtitle)

    plt.savefig(f"{imgdir}/basecov_raw_{scaffold_key}.png")
    plt.close()






