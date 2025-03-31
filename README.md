# XML2FASTA

A Python tool designed to extract all variant sequences from UniProt (variant) XML files to FASTA files.  

## Dependencies

This script requires Python 3.7+ and the following packages (see installation section on how to install).

- `Bio`
- `biopython`
- `numpy`
- `Requests`

## Installation

1. Clone this repository:

```
git clone https://github.com/leontwntfr/XML2FASTA.git
cd XML2FASTA
```

2. Install the required dependencies using ``pip``:

```
pip install -r requirements.txt
```

3. Run the tool:

```
python XML.py
```

## Usage

The tool provides an user interface where input files (multiple are possible) and an output directory can be selected. The following input files are allowed:
- EBL variant XML files (e.g. obtained from UniProt variant viewer)
- UniProt entry XML files (also including multiple entries)

The "Submit" button starts the process and the generated FASTA files are saved in the indicated output directory.

In the case of providing UniProt entry XML files, an internet connection is required as the associated variant XML files are being downloaded first when pressing "Submit".
