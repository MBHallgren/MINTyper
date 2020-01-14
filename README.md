Read more about Mint3 here:
> PUBLICATION TBA

# Table of contents

* [Introduction](#introduction)
* [Requirements](#requirements)
* [Installation](#installation)
* [Options and usage](#options-and-usage)
    * [Functions](#functions)
    * [Standard usage](#standard-options)
# Introduction
Mint3 is a tools designed to quickly a distance matrix from a set of input sequencing reads. It can take both short read illumina sequences and long read Oxford Nanopore sequences. 

# Requirements

* Linux or macOS
* [Python](https://www.python.org/) 3.4 or later
* [KMA](https://bitbucket.org/genomicepidemiology/kma) As new as possible, but most should work
* [ccphylo](https://bitbucket.org/genomicepidemiology/ccphylo/src/master/) As new as possible, but most should work


# Installation
The following instructions will install the latest version of Mint3:

First, clone the github repository:

"git clone https://github.com/s153002/Mint3.git"

Then run the  installer:

"python3 mintinstall.py"

### Standard usage

Make sure to index your reference database with "kma index -Sparse ATG" to make the reference finding faster! (See https://bitbucket.org/genomicepidemiology/kma)

Use the dcmMethylations file as input to mask DCM motifs (CC(A/T)GG). Additionally, feel free to add other motifs in a multifasta format in the same file, if they are to be masked.

Make sure to move all your isolate files to an empty directory. Nanopore and Illumina files should be placed in seperated directories.

An example  os a standard usage could be:

"Mint3 -i_path_illumina /home/usr/illuminaDirectory/ -i_path_nanopore /home/usr/NanoporeDirectory/ -dcm /home/usr/currentDir/dcmFile -o run1 -db /home/user/databases/Bacteria.ATG -thread 6"

