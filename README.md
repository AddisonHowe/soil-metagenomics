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
