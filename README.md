# EZmito2 🧬

<div align="center">

![Python](https://img.shields.io/badge/python-3.7+-blue.svg)
![License](https://img.shields.io/badge/license-MIT-green.svg)
![Bioinformatics](https://img.shields.io/badge/field-Bioinformatics-purple.svg)
![Status](https://img.shields.io/badge/status-active-success.svg)

**A comprehensive toolkit for mitochondrial genome analysis**

[Features](#features) • [Installation](#installation) • [Usage](#usage) • [Web Server](#web-server) • [Citation](#citation)

</div>

---

## 🌟 Overview

**EZmito2** is a simple and fast command-line tool designed for comprehensive mitochondrial genome analyses. It provides a suite of utilities for processing, analyzing, and visualizing mitochondrial DNA sequences with a focus on ease of use and efficiency.

> 💻 **Web Server Available!** A parallel web-based version is available at [ezmito.unisi.it](http://ezmito.unisi.it)

---

## ✨ Features

EZmito2 includes seven powerful subcommands for different mitochondrial genome analyses:

### 🔄 **EZcircular**
Circularize mitochondrial sequences starting from a specific gene
- Supports both circular and linear genomes
- Customizable starting gene position
- Handles BED annotation files

### 🧪 **EZcodon**
Analyze codon usage patterns in heavy (J) and light (N) strands
- Calculate RSCU (Relative Synonymous Codon Usage)
- Amino acid frequency analysis
- Beautiful visualization plots
- Support for multiple genetic codes

### 🗺️ **EZmap**
Create custom circular or linear maps of mitochondrial genomes
- Circular and linear plot generation
- Customizable color schemes
- GFF3 annotation support
- Publication-ready figures

### 🔍 **EZmix**
Detect potential chimeras in mitochondrial assemblies
- BLAST-based similarity detection
- Adjustable identity and length thresholds
- Visual representation of homologous regions

### 🧬 **EZpipe**
Prepare mitochondrial sequences for phylogenetic analysis
- Automated alignment with MAFFT
- Gblocks trimming
- Codon position filtering (2nd/3rd position removal)
- PartitionFinder configuration generation
- NEXUS and PHYLIP output formats

### 📊 **EZskew**
Analyze nucleotide skew biases in mitochondrial strands
- AT, CG, AC, and GT skew calculations
- Position-specific analysis (1st, 2nd, 3rd codon positions)
- Comprehensive visualization
- Support for multiple genetic codes

### ✂️ **EZsplit**
Extract individual protein-coding genes from complete mitogenomes
- Batch processing of multiple genomes
- GFF3-based gene extraction
- Automatic gene name standardization
- Missing gene detection

- Splits sequences into transmembrane (TM), matrix (MA), and inner membrane (IM) domains following the  pipeline
### 🔀 **EZtrampo**
Partition aligned mitochondrial PCGs into transmembrane, matrix, and inner membrane domains, following the [TRAMPO](https://github.com/dbajpp0/TRAMPO) pipeline
- Domain boundary inference using reference models
- Support for custom reference organisms via user-supplied sequences and domain tables
- Nexus partition files by gene, codon position, and structural domain
- Composition plots and tables for skews, RSCU, and physicochemical properties

---

## 📋 Requirements

### System Requirements
- Linux/Unix-based operating system
- Conda (Miniconda or Anaconda)
- 4GB RAM minimum
- 2GB disk space

### Software Dependencies
All dependencies are automatically installed via the provided conda environment (some examples here, a better description in the YML file):
- Python 3.7+
- Biopython
- pandas
- matplotlib
- numpy
- pyfiglet
- pycirclize
- BCBio
- itaxotools.pygblocks
- MAFFT
- BLAST+

---

## 🚀 Installation

### Quick Install (Recommended)

EZmito2 provides an automated installation script that sets up everything for you:
```bash
# Clone the repository
git clone https://github.com/yourusername/ezmito.git
cd ezmito

# Run the installer
bash install.sh
```

The installer will:
1. ✅ Check if Conda is installed (and offer to install it if missing)
2. ✅ Create the `ezmito_env` conda environment
3. ✅ Install all required dependencies
4. ✅ Set up executable permissions

### Manual Installation

If you prefer to install manually:
```bash
# Clone the repository
git clone https://github.com/yourusername/ezmito.git
cd ezmito

# Make files executable
chmod 775 ezmito.yml
chmod 775 ezmito.py

# Create conda environment
conda env create -f ezmito_env.yml

# Activate the environment
conda activate ezmito_env
```

### Verify Installation

After installation, verify everything is working:
```bash
# Activate the environment
conda activate ezmito_env

# Check EZmito version and help
python ezmito.py --help
```

---

## 💻 Usage

### Activating the Environment

Before running EZmito2, always activate the conda environment:
```bash
conda activate ezmito_env
```

### General Syntax
```bash
python ezmito.py <command> [options]
```


### Available Commands

#### 🔄 EZcircular
```bash
python ezmito.py ezcircular -i input.fasta -b input.bed -s cox1 -f circular -o output_dir
```

**Options:**
- `-i, --input`: Input FASTA file (required)
- `-b, --bed`: Input BED file (required)
- `-s, --start`: Starting gene (default: cox1)
- `-f, --feature`: Genome feature - circular or linear (default: circular)
- `-o, --outdir`: Output directory (default: outdir)

---

#### 🧪 EZcodon
```bash
python ezmito.py ezcodon -J heavy_genes/ -N light_genes/ -c 2 -o output_dir
```

**Options:**
- `-J, --heavy`: Path to heavy chain FASTA files
- `-N, --light`: Path to light chain FASTA files
- `-c, --code`: Genetic code (required)
- `-o, --outdir`: Output directory (default: outdir)

**Supported Genetic Codes:**
- `2` - Vertebrate Mitochondrial Code
- `4` - Mold, Protozoan, and Coelenterate Mitochondrial Code
- `5` - Invertebrate Mitochondrial Code
- `9` - Echinoderm and Flatworm Mitochondrial Code
- `13` - Ascidian Mitochondrial Code
- `14` - Alternative Flatworm Mitochondrial Code
- `21` - Trematode Mitochondrial Code
- `33` - Cephalodiscidae Mitochondrial UAA-Tyr Code

---

#### 🗺️ EZmap
```bash
python ezmito.py ezmap -g input.gff -f circular -colJ '#add8e6' -colN '#B22222' -o output_dir
```

**Options:**
- `-g, --gff`: Input GFF3 file (required)
- `-f, --feature`: Genome feature - circular or linear (default: circular)
- `-colJ, --colorJ`: Heavy strand color (default: #add8e6)
- `-colN, --colorN`: Light strand color (default: #B22222)
- `-o, --outdir`: Output directory (default: outdir)

---

#### 🔍 EZmix
```bash
python ezmito.py ezmix -i sequences.fasta -id 0.95 -len 200 -o output_dir
```

**Options:**
- `-i, --input`: MultiFASTA input file (required)
- `-id, --identity`: Identity threshold 0.5-1 (default: 0.95)
- `-len, --length`: Length threshold in bp (default: 200)
- `-bn, --blastn`: Path to blastn executable (default: uses conda installation)
- `-o, --outdir`: Output directory (default: outdir)

---

#### 🧬 EZpipe
```bash
python ezmito.py ezpipe -i genes/ -c 5 -p 3 -o output_dir
```

**Options:**
- `-i, --input`: Path to FASTA files (required)
- `-c, --code`: Genetic code (required)
- `-p, --positions`: Number of codon positions - 2 or 3 (default: 3)
- `-o, --outdir`: Output directory (default: outdir)

---

#### 📊 EZskew
```bash
python ezmito.py ezskew -J heavy_genes/ -N light_genes/ -c 2 -o output_dir
```

**Options:**
- `-J, --heavy`: Path to heavy chain FASTA files
- `-N, --light`: Path to light chain FASTA files
- `-c, --code`: Genetic code (required)
- `-o, --outdir`: Output directory (default: outdir)

---

#### ✂️ EZsplit
```bash
python ezmito.py ezsplit -i mitogenomes.fasta -g annotations.gff -o output_dir
```

**Options:**
- `-i, --input`: MultiFASTA input file of complete mitogenomes (required)
- `-g, --gff`: GFF3 input file (required)
- `-o, --outdir`: Output directory (default: outdir)


#### 🔀 EZtrampo
```bash
python ezmito.py eztrampo -p genes/ -c 5 -m dme -g panc -n 4 -o output_dir
```
**Options:**
- `-p, --path`: Path to FASTA files directory (required)
- `-c, --code`: Genetic code (required)
- `-m, --model`: Model organism nickname (required)
- `-g, --gene_order`: Gene order model nickname
- `-s, --sequence`: User's model organism genes in FASTA (amino acid) format
- `-t, --tables`: User's model organism table files in TMHMM format
- `-n, --threads`: Number of threads for MAFFT (default: 1)
- `-o, --outdir`: Output directory (default: outdir)

**Model Organisms:**
- `hsa` - *Homo sapiens* (Chordata)
- `ppe` - *Patiria pectinifera* (Echinodermata)
- `dme` - *Drosophila melanogaster* (Pancrustacea+Chelicerata)
- `aca` - *Albinaria caerulea* (Mollusca)
- `lte` - *Lumbricus terrestris* (Annelida)
- `cel` - *Caenorhabditis elegans* (Nematoda)
- `mse` - *Metridium senile* (Cnidaria)
- `user` - Custom model (requires `-s` and `-t`)

**Gene Order Models:**
- `vert` - All vertebrates
- `panc` - Arthropods
- `ances` - *Lumbricus terrestris*, *Caenorhabditis elegans*, *Metridium senile*
- `albin` - *Albinaria caerulea*
- `meta` - *Metacrangonyx*

**Supported Genetic Codes:**
- `2` - Vertebrate Mitochondrial Code
- `3` - Yeast Mitochondrial Code
- `4` - Mold, Protozoan, and Coelenterate Mitochondrial Code
- `5` - Invertebrate Mitochondrial Code
- `9` - Echinoderm and Flatworm Mitochondrial Code
- `13` - Ascidian Mitochondrial Code
- `14` - Alternative Flatworm Mitochondrial Code
- `16` - Chlorophycean Mitochondrial Code
- `21` - Trematode Mitochondrial Code
- `22` - Scenedesmus obliquus Mitochondrial Code
- `23` - Thraustochytrium Mitochondrial Code
- `33` - Cephalodiscidae Mitochondrial UAA-Tyr Code
- 
---

## 🌐 Web Server

Don't want to use the command line? Try our **user-friendly web interface** at:

### [🔗 ezmito.unisi.it](http://ezmito.unisi.it)

The web server provides:
- ✅ No installation required
- ✅ Intuitive graphical interface
- ✅ Example datasets
- ✅ Direct download of results
- ✅ Job management system

---

## 📊 Output Files

Each tool generates specific outputs:

- **EZcircular**: Rearranged FASTA and BED files
- **EZcodon**: RSCU tables, amino acid frequency plots, codon usage plots
- **EZmap**: Publication-ready circular or linear genome maps (PDF)
- **EZmix**: Similarity plots showing potential chimeric regions (PDF)
- **EZpipe**: Aligned sequences, partitioned datasets, PartitionFinder config
- **EZskew**: Skew analysis tables and plots (CSV, PDF)
- **EZsplit**: Individual gene FASTA files, missing genes report
- **EZtrampo**: Partition files, plots (PDF, HTML, CSV), statistic tables, alignment

---

## 🛠️ Troubleshooting

### Common Issues

**Issue: "Conda not found"**
```bash
# The installer will offer to install Conda automatically
# Or install manually from: https://docs.conda.io/en/latest/miniconda.html
```

**Issue: "Environment already exists"**
```bash
# Remove the existing environment
conda env remove -n ezmito_env

# Reinstall
conda env create -f ezmito_env.yml
```

**Issue: "Permission denied"**
```bash
# Make sure files are executable
chmod 775 install.sh
chmod 775 ezmito.py
```

---

## 📖 Citation

If EZmito2 helps your research, while eaiting for the new manuscript publication, please cite:

> **Cucini C., Leo C., Iannotti N., Boschi S., Brunetti C., Pons J., Fanciulli P. P., Frati F., Carapelli A., & Nardi F. (2021)**  
> *EZmito: a simple and fast tool for multiple mitogenome analyses*  
> Mitochondrial DNA Part B, 6(3), 1101-1109.  
> DOI: [10.1080/23802359.2021.1899865](https://doi.org/10.1080/23802359.2021.1899865)

---

## 📝 License

This project is licensed under the MIT License - see the LICENSE file for details.

---

## 👥 Authors

- Cucini C., Leo C., Iannotti N., Boschi S., Brunetti C., Pons J., Fanciulli P. P., Frati F., Carapelli A., & Nardi F.

---

## 🐛 Issues & Support

Found a bug or need help? Please open an issue on our [GitHub Issues page](https://github.com/yourusername/ezmito/issues).

---

## 📚 Documentation

For more detailed documentation, visit the tutorial at:
- 🎓 [Tutorial Videos](#) (coming soon)

---

<div align="center">

⭐ Star us on GitHub if you find EZmito2 useful!

</div>
