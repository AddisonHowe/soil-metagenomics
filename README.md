# Metagenomic soil analysis

## Installation

From the project directory, create a conda environment for the project as follows:

```bash
conda env create -p ./env -f environment.yml
```

Next install the project source code.
From the project directory, activate the environment and run:

```bash
conda activate env
python -m pip install -e '.[dev]'
```

Verify things have installed successfully by running:

```bash
pytest tests
```

## Links

- [NCBIFAM entry nitrate reductase subunit alpha](https://www.ebi.ac.uk/interpro/entry/ncbifam/TIGR01580/)
- [InterPro entry nitrate reductase, alpha subunit](https://www.ebi.ac.uk/interpro/entry/InterPro/IPR006468/)
- [Pfam entry PF00384 (Molybdopterin oxidoreductase)](https://www.ebi.ac.uk/interpro/entry/pfam/PF00384/)
- [Pfam entry PF01568 (Molydopterin dinucleotide binding domain)](https://www.ebi.ac.uk/interpro/entry/pfam/PF01568/)
