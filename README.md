# PAMPHLET-LOCAL: PAM Prediction HomoLogous-Enhancement Toolkit

[PAMPHLET-LOCAL](https://github.com/GolshkovQ/PAMPHLET-LOCAL/tree/main) is a powerful tool designed for precise PAM (Protospacer Adjacent Motif) prediction in CRISPR-Cas systems. It utilizes homologous CRISPR-Cas searches and spacer expansion to enhance prediction accuracy. PAMPHLET-LOCAL offers flexibility and customization through various parameters to suit a wide range of research needs in CRISPR-Cas studies.

## Features

- **PAM Prediction**: Efficiently predicts PAM sequences in CRISPR-Cas systems.
- **Homologous Search**: Utilizes homologous CRISPR-Cas searches to improve accuracy.
- **Customizable Parameters**: Offers a wide array of parameters for fine-tuning the prediction process.
- **Flexible Output**: Outputs results in a user-specified directory for easy management.
- **Multi-threading Support**: Supports multi-threading to speed up BLAST searches.

## Installation

1. Clone the repository:

    ```bash
    git clone https://github.com/GolshkovQ/PAMPHLET-LOCAL/tree/main
    ```

2. Install required dependencies:

    ```bash
    pip install argparse biopython pyfaidx pickle urllib3 requests func_timeout
    ```

## Usage

Run PAMPHLET-LOCAL with the following required and optional arguments:

```bash
python pamphlet.py -s <spacer.fa> -r <repeat_sequence> -p <protein.fa> -o <output_directory> -bDB <bacteria_database> [options]
```

## Arguments

### Required Arguments:
- `-s, --spacer`: Spacer sequences (FASTA format), required.
- `-r, --repeat`: Repeat sequence for revising the flank sequence, required.
- `-p, --protein`: Class II Cas effector protein sequence (FASTA format), required.
- `-o, --outdir`: Output directory, required.
- `-bDB, --bacteriaDB`: Bacteria database for homologous CRISPR-Cas search and spacer expansion, required.

### Optional Arguments:
- `-pDB, --protoDB`: Protospacer source genome database for PAMPHLET-LOCAL, default is an empty string.
- `-cDB, --casDB`: Cas protein database for PAMPHLET-LOCAL, default is an empty string.
- `-L, --flanklen`: Length of flank sequence, default is `12`.
- `-b, --blastmode`: Spacer BLASTN mode, either `relax` or `common` (default is `common`).
- `--pcovs`: Minimum percent coverage of spacer sequence, default is `0.9`.
- `--pident`: Minimum percent identity of spacer sequence, default is `0.9`.
- `--rident`: Minimum percent identity of repeat sequence, default is `0.8`.
- `--reviseLen`: Revise length of spacer sequence (default is `False`).
- `-f, --freqmode`: Base frequency calculation mode, either `sigmoid` or `linear` (default is `sigmoid`).
- `-O, --orientation`: Orientation of the repeat sequence, either `positive` or `negative` (default is `positive`).
- `-bt, --BlastThreads`: Number of threads for BLASTN, default is `1`.

## Example

```bash
python pamphlet.py \
    -p /path/to/protein.fa \
    -o /path/to/output_directory \
    -r GTTTTGGTAGCATTCAAAATAACATAGCTCTAAAAC \
    -s /path/to/spacer.fa \
    -bDB /path/to/bacteria_database \
    -cDB /path/to/cas_database \
    -pDB /path/to/protospacer_database \
    -b relax --pcovs 0.9 --pident 0.9 --rident 0.9 -f sigmoid -O positive -bt 10
```

## Citation

If you use PAMPHLET-LOCAL in your research, please cite the following article:

**PAMPHLET: PAM Prediction HomoLogous-Enhancement Toolkit for Precise PAM Prediction in CRISPR-Cas Systems**  
Available on bioRxiv: [https://www.biorxiv.org/content/10.1101/2024.04.09.587696v2](https://www.biorxiv.org/content/10.1101/2024.04.09.587696v2)  
DOI: [10.1101/2024.04.09.587696v2](https://doi.org/10.1101/2024.04.09.587696v2)

