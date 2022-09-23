Read more about mintyper here:
> https://academic.oup.com/biomethods/article/6/1/bpab008/6243723

# Table of contents

* [Introduction](#introduction)
* [Requirements](#requirements)
* [Installation](#installation)
* [Database](#database)
* [Standard usage](#standard-options)
* [Output](#output)
* [Licens](#licens)
    
# Introduction
mintyper is a tool designed to quickly estimate a distance matrix from a set of input sequencing reads. It can take both short read from Illumina's sequencing platforms and long reads from Oxford Nanopore's platform. 

# Requirements

* Unix based system
* [Python3](https://www.python.org/) Any version of python3
* [KMA](https://bitbucket.org/genomicepidemiology/kma) Installation requirements for KMA installation can be found at the bitbucket page.
* [ccphylo](https://bitbucket.org/genomicepidemiology/ccphylo/src/master/) Installation requirements for ccphylo installation can be found at the bitbucket page.
* (optional) IQtree installed acording to http://www.iqtree.org/doc/Quickstart. Using IQTree requires you to have an executable file for IQtree.
* (optional) Fasttree installed acording to http://www.microbesonline.org/fasttree/#Install. Using FastTree required a root installation and NOT an executable file.

# Installation
The following commands will install the latest version of mintyper if python3, a C-compiler and zlib development files are installed (See KMA bitbucket).:

```bash
#For a system-wide installation:
git clone https://bitbucket.org/genomicepidemiology/mintyper.git
cd mintyper
python3 setup.py build
python3 setup.py install

#If setuptools installation (the above) is not an option, instead run:
git clone https://bitbucket.org/genomicepidemiology/mintyper.git
cd mintyper
python3 local_install.py
```


# Database

mintyper provides two options: Either the user can give a single fasta file as reference sequence by using the -ref argument, or a whole KMA indexed database can be used with the -db argument.

A guide for indexing a given database correctly can be found at https://bitbucket.org/genomicepidemiology/kma/src/master/

A database of complete bacterial genomes can be downloaded with: 
```bash
wget ftp://ftp.cbs.dtu.dk/public/CGE/databases/KmerFinder/version/latest/bacteria.tar.gz .
```
# Standard usage and testrun

To make sure the installation has been completed succesfully, run mintyper on the assembled MinION data found in the testrun folder:
```bash
"mintyper -i_assemblies testrun/data/* -ref testrun/template_sequence -o output"
#Or for the locally installed version:
"python3 mintyper_local.py -i_assemblies testrun/data/* -ref testrun/template_sequence -exe_path /path/to/installationdir/ -o output"
```



Make sure to index your reference database with "kma index -Sparse ATG" to make the reference finding faster! (See https://bitbucket.org/genomicepidemiology/kma)

Use the dcmMethylations file as input to mask DCM motifs (CC(A/T)GG). Additionally, feel free to add other motifs in a multifasta format in the same file, if they are to be masked. A masking file for multiple motifs should look like:
```bash
>Header1
motif_1_sequence
>Header2
motif_2_sequence
>Header3
motif_3_sequence
```
Make sure to move all your isolate files to an empty directory. Nanopore and Illumina files should be placed in seperated directories, because the KMA alignments are different for the two technologies.

An example  of a standard usage could be:
```bash
mintyper -i_illumina /home/usr/illuminaDirectory/ -i_nanopore /home/usr/NanoporeDirectory/ -masking_scheme /home/usr/currentDir/dcmFile -prune_distance 5 -db /home/user/databases/Bacteria.ATG -thread 6 -o output"
```

Run `mintyper --help` to see the program's usage options.

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
