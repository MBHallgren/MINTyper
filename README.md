Read more about mintyper here:
> https://academic.oup.com/biomethods/article/6/1/bpab008/6243723

# Table of contents

* [Introduction](#introduction)
* [Installation](#installation)
* [Standard usage](#standard-options)
* [Output](#output)
* [Licens](#licens)
    
# Introduction
Mintyper is a tool designed to quickly estimate a distance matrix from a set of input sequencing reads. It can take both short read from Illumina's sequencing platforms and long reads from Oxford Nanopore's platform. 

# Installation

## Dependencies

Mintyper requires the following dependencies:

* Python 3.6 or higher
* KMA (https://bitbucket.org/genomicepidemiology/kma/src/master/)
* CCphylo (https://bitbucket.org/genomicepidemiology/ccphylo/src/master/)

### Dependencies installation

#### Conda (works on linux only due to KMA conda package)

```bash
git clone https://github.com/mbhallgren/mintyper
cd mintyper
conda env create -f environment.yml -n mintyper
conda activate mintyper
```

#### Manual installation (Linus and MacOS, requires libz-dev installed)

```bash
git clone https://github.com/MBHallgren/MINTyper.git
cd mintyper
python3 scripts/install_dependencies.py
```

## Mintyper installation

### Using pip

The easiest way to install mintyper is using pip:

```bash
pip install mintyper
```

### From source

```bash
git clone https://github.com/MBHallgren/MINTyper.git
cd mintyper
python3 setup.py install
```

## Database

Mintyper provides two options: Either the user can give a single fasta file as reference sequence by using the -ref argument, or a whole KMA indexed database can be used with the -db argument.

A guide for indexing a given database correctly can be found at https://bitbucket.org/genomicepidemiology/kma/src/master/

A database of complete bacterial genomes can be downloaded with: 
```bash
wget https://cge.food.dtu.dk/services/MINTyper/bac_species_db.tar.gz .
```

# Standard usage

```bash
mintyper --illumina /path/to/illlumina_files/* --nanopore /path/to/nanopore_files/* --db /path/to/database --output output_folder_name
```

# Output

The output consists of a distance matrix in relaxed phylip format (distmatrix.phy), a list of template candidates (Template_kma_results.spa), the chosen template sequence (template_sequence.fasta), the combined VCF file of all analyzed reads (Combined.vcf.gz), a newick file of the generated phylogenetic tree (outtree.newick), a cluster file (Optional, only outputted if -cluster_length > 0) and a logfile. Additionally, a data folder containing the consensus sequences, the individual VCFs and alignment files is created.

# Citation #
1. Malte B Hallgren, Søren Overballe-Petersen, Ole Lund, Henrik Hasman, Philip T L C Clausen, mintyper: an outbreak-detection method for accurate and rapid SNP typing of clonal clusters with noisy long reads, Biology Methods and Protocols, Volume 6, Issue 1, 2021, bpab008, https://doi.org/10.1093/biomethods/bpab008

# Licens
```bash
Copyright (c) 2021, Malte Hallgren, Technical University of Denmark All rights reserved.

Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License. You may obtain a copy of the License at

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.
```

# Additional information
As of April 2023 IQtree and FastTree are no longer supported automatically. If you wish to use these programs, please install them manually and use the multiple alignment as input.
