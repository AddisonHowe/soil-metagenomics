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
- [Nitrate reduction Ortholog Table](https://www.google.com/url?q=https://www.kegg.jp/kegg-bin/view_ortholog_table?md%3DM00529&source=gmail&ust=1756589228072000&usg=AOvVaw2vaH8k1EgoPI7un7sTLuIu)
- [Nitrite reduction Ortholog Table](https://www.kegg.jp/kegg-bin/view_ortholog_table?md=M00530)
- [Read this regarding enzymes](https://pmc.ncbi.nlm.nih.gov/articles/PMC4453514/)
- [KEGG listings](https://www.kegg.jp/kegg-bin/get_htext?ko00001+K02591)
