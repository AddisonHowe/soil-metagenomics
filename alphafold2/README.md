# AlphaFold Protein Analysis

## Running AlphaFold

To predict protein structure, AlphaFold just needs an amino acid sequence in fasta format.
On Quest, we can run alphafold2 in two stages, the CPU and then the GPU stage.
Each stage takes on the order of a couple hours per sequence (for sequences around 1000 amino acids in length).
Multiple sequences can be processed at once, by specifying a comma-separated list of absolute paths as the `fasta_files` argument.
Each fasta file should contain only one sequence.

The scripts `alphafold2/run_alphafold_[cpu,gpu].sh` will run the relevant step, and take as input the single or comma-separated fasta file list (in absolute path) and the specified output directory (also as an absolute path).
A subdirectory will then be created for each fasta file, within the output directory.

The CPU process must finish before the GPU process can be started.

## Visualizing AlphaFold's results

We can visualize the results of AlphaFold using [pymol](https://www.pymol.org/).
We'll use the [open-source pymol](https://github.com/schrodinger/pymol-open-source) distribution, available through conda, that runs on Macs (Arm64 and Aarch64) and linux machines.
Create a conda environment and install pymol:

```bash
conda create -n pymol-env
conda activate pymol-env
conda install conda-forge::pymol-open-source
```

Once installed, from the terminal, run `pymol` and the GUI should appear.

### Using the pymol GUI

Pymol allows us to visualize the .pdb files that AlphaFold produces.
A useful tutorial can be found [here](https://rtguides.it.tufts.edu/bio/lectures/introduction-to-alphafold2/03-vizualize-predicted-structure.html).
We can load a file via the menu bar: file>open.

### Commands

#### Coloring

The color command allows us to color by the pLDDT confidence value, which is the `b` column in the .pdb output file.

```sh
spectrum b, red blue, minimum=0, maximum=100  # color from red (worst) to blue (best)
```

Or for a particular color for certain structures:

```sh
color cyan
```

#### Aligning proteins

After loading two pdb files, `<protein1.pdb>` `<protein2.pdb>`, we can align them using:

```sh
align <protein1>, <protein2>
```

#### Fetching proteins from the PDB

Use the `fetch` command to load existing proteins from the [PDB](https://www.rcsb.org), the [AlphaFold protein database](https://alphafold.ebi.ac.uk/), or [UniProt](https://www.uniprot.org/uniprotkb?query=*), based on PDB ID.

```sh
fetch <PDB_ID>
```

For proteins with multiple complexes, the remove command allows us to focus on only certain parts.
For example, using NarG:

```sh
fetch 1Q16
remove not chain A  # Chain A is the NarG subunit of the larger complex
```

### Running pymol scripts

We can run pymol scripts as follows (from [pymol's tutorial](https://pymol.org/tutorials/scripting/howtorunscripts.html)):

```bash
# run a PyMOL command script (opens PyMOL GUI)
pymol script.pml

# run a Python script
pymol script.py

# batch mode (no PyMOL GUI)
pymol -c script.pml

# use absolute path if "pymol" is not in $PATH
/opt/pymol-2.2.3/pymol -c script.pml
```

## Mounting remote directories (MacOS)

We can use `sshfs` to mount a remote directory, and access the data stored there without having to download anything to a local machine.
On Macs, a required software is [MacFuse](https://macfuse.github.io/), and full instructions can be found [here](Instructions can be found [here](https://phoenixnap.com/kb/sshfs-mac).)

Once installed, create a mount point directory somewhere convenient:

```bash
cd ~
mkdir mount_point
```

Then use sshfs to mount a remote directory, which will be accessible via the mount point:

```bash
sshfs <user>@<host>:</path/to/directory> mount_point
```

This will enable access to the directory on the remote machine, and running `ls mount_point` should show the contents of the remote directory.
This allows us to load the large .pdb files on Quest into pymol locally, so as to use the GUI capabilities.

To unmount the filesystem, eject the relevant macFUSE volume through finder or use the `umount` command.
