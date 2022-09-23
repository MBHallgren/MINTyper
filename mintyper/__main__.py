"""
__main__ for mintyper
"""

#!/usr/bin/env python3
import argparse
import sys
from mintyper_functions import mintyper


def main():
    """
    main call
    :return:
    """
    args = parse_args(sys.argv[1:])
    mintyper(args)

def parse_args(args):
    """
    Arguments parser for running mintyper locally
    :param args:
    :return:
    """
    description = 'mintyper: an outbreak-detection method for accurate \
     and rapid SNP typing of clonal clusters with noisy long reads.'
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-i_illumina', action="store", type=str, dest='i_illumina', nargs="+",
                        default=[], help='Illumina input files')
    parser.add_argument('-i_nanopore', action="store", type=str, dest='i_nanopore', nargs="+",
                        default=[], help='Nanopore input files')
    parser.add_argument('-i_assemblies', action="store", type=str, dest='i_assemblies', nargs="+",
                        default=[],help='Assembly input files. Assemblies should not be'
                                        ' mixed with raw reads in the same analysis.')
    parser.add_argument("-pe", action="store_true", dest="paired_end", default=False,
                        help="If paipred ends are used for illumina reads,"
                             " otherwise they will be considered at individual isolates.")
    parser.add_argument("-masking_scheme", type=str, action="store", dest="masking_scheme",
                        default="", help="Give a fasta file containing a motif"
                                         ' that is to be masked in the aligned concensus files.')
    parser.add_argument("-prune_distance", type=int, action="store", dest="prune_distance",
                        default=10, help="X length that SNPs can be located between each other."
                                         " Default is 10."
                                         " If two SNPs are located within X length of each other,"
                                         " everything between them as well as X length on each side"
                                         " of the SNPs will not be used in the alignments"
                                         " to calculate the distance matrix.")
    parser.add_argument("-insig_prune", action="store_true", dest="insig_prune", default=False,
                        help="By default insignificant bases are included in the pruning process."
                             " Use this flag, if you wish to NOT include them in the pruning.")
    parser.add_argument("-bc", action="store", type=float, default=0.7, dest="bc",
                        help="Base calling parameter for nanopore KMA mapping. Default is 0.7.")
    parser.add_argument("-db", action="store", type=str, default="", dest="ref_kma_database",
                        help="Comeplete path to a KMA indexed reference database.")
    parser.add_argument("-thread", action="store", default=1, dest="multi_threading",
                        help="Set this parameter to x-number of threads KMA should use.")
    parser.add_argument("-ref", action="store", type=str, default="", dest="reference",
                        help="To align your query sequences against a reference"
                             " of your own choice, use this argument.")
    parser.add_argument("-cluster_length", type=int, action="store", dest="cluster_length",
                        default=-1, help="Maximum distance within the same cluster."
                                         " Default is -1, thus not estimating any clusters."
                                         " Set to positive integer to generate a cluster file.")
    parser.add_argument("-iqtree", action="store_true", default=False, dest="iqtree",
                        help="To use iqtree if it is installed."
                             " Default clustering is ccphylo (preinstalled).")
    parser.add_argument("-fast_tree", action="store_true", dest="fast_tree", default=False,
                        help="To use FastTree if it is installed."
                             " Default clustering is ccphylo (preinstalled)..")
    parser.add_argument('-version', action='version', version='mintyper 1.1.0',
                        help="current version of mintyper")
    parser.add_argument("-exe_path", action="store", dest="exe_path", default="", required=True,
                        help="Complete path to the mintyper repo that you cloned."
                             " It should contain a kma folder and a ccphylo folder.")
    parser.add_argument("-o", action="store", dest="output_name",
                        help="Output directory name.")
    parser.add_argument("-cge", action="store_true", dest="cge",
                        help="Does nothing. Internal webserver call for CGE.")
    args = parser.parse_args()

    return args


if __name__ == '__main__':
    main()
