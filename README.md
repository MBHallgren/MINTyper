Read more about MINTyper here:
> PUBLICATION TBA

# Table of contents

* [Introduction](#introduction)
* [Requirements](#requirements)
* [Installation](#installation)
* [Database](#database)
* [Standard usage](#standard-options)
    
# Introduction
MINTyper is a tools designed to quickly a distance matrix from a set of input sequencing reads. It can take both short read from Illumina's sequencing platforms and long reads from Oxford Nanopore's platform. 

# Requirements

* Linux or macOS
* [Python3](https://www.python.org/) Any version of python3
* [KMA](https://bitbucket.org/genomicepidemiology/kma) Installation requirements for KMA installation can be found at the bitbucket page.
* [ccphylo](https://bitbucket.org/genomicepidemiology/ccphylo/src/master/) Installation requirements for ccphylo installation can be found at the bitbucket page.


# Installation
The following commands will install the latest version of MINTyper if python3, a C-compiler and zlib development files are installed (See KMA bitbucket).:

```bash
git clone https://github.com/s153002/MINTyper.git
cd MINTyper
python3 MINTyperInstall.py
```


# Database

MINTyper provides two options: Either the user can give a single fasta file as reference sequence by using the -ref argument, or a whole KMA indexed database can be used with the -db argument.

A guide for indexing a given database correctly can be found at https://bitbucket.org/genomicepidemiology/kma/src/master/

A database of complete bacterial genomes can be found at: http://www.cbs.dtu.dk/public/CGE/databases/KmerFinder/version/latest/bacteria.tar.gz 

# Standard usage

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
"python3 MINTyper.py -i_path_illumina /home/usr/illuminaDirectory/ -i_path_nanopore /home/usr/NanoporeDirectory/ -masking_scheme /home/usr/currentDir/dcmFile -prune_distance 5 -db /home/user/databases/Bacteria.ATG -thread 6 -o output"
```

Run `MINTyper.py --help` to see the program's usage options:

```
usage: MINTyper.py [-h] [-i_path_illumina I_PATH_ILLUMINA]
                   [-i_path_nanopore I_PATH_NANOPORE] [-pe]
                   [-masking_scheme MASKING_SCHEME]
                   [-prune_distance PRUNE_DISTANCE] [-bc BC]
                   [-db REF_KMA_DATABASE] [-thread MULTI_THREADING]
                   [-ref REFERENCE] [-version] [-exepath EXEPATH]
                   [-o OUTPUT_NAME]

.

optional arguments:
  -h, --help            show this help message and exit
  -i_path_illumina I_PATH_ILLUMINA
                        The path to the directory containing ONLY the input
                        illumina files. Should be used when analyzing >5 read-
                        files at a time.
  -i_path_nanopore I_PATH_NANOPORE
                        The path to the directory containing ONLY the input
                        nanopore files. Should be used when analyzing >5 read-
                        files at a time.
  -pe                   If paipred ends are used give input as True (-pe
                        True). If Paired-ends are used, it is important that
                        the files are written in the correct order, such as:
                        sample1_1.fasta sample1_2.fasta sample2_1.fasta
                        sample2_1.fasta
  -masking_scheme MASKING_SCHEME
                        Give a fasta file containing a motof that you wish to
                        mask in the aligned concensus files.
  -prune_distance PRUNE_DISTANCE
                        X lenght that SNPs can be located between each other.
                        Default is 10. If two SNPs are located within X lenght
                        of eachother, everything between them as well as X
                        lenght on each side of the SNPs will not be used in
                        the alignments to calculate the distance matrix.
  -bc BC                Base calling parameter for nanopore KMA mapping.
                        Default is 0.7
  -db REF_KMA_DATABASE  Comeplete path for the ref_kma_database for KMA
                        mapping
  -thread MULTI_THREADING
                        Set this parameter to x-number of threads that you
                        would like to use during KMA-mapping.
  -ref REFERENCE        KMA will by default determine the best template
                        against the given database. However, if you want to
                        align your query sequences against a reference of your
                        own choice, use this function. If this is left blank,
                        KMA will determine the optimal reference.
  -version              current version of MINTyper
  -exepath EXEPATH      Complete path to the MINTyper repo that you cloned,
                        in which the executables are located
  -o OUTPUT_NAME        Name that you would like the output directory to be
                        called.
```
