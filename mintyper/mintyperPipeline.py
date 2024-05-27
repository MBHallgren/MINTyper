import os
import sys
import logging
import multiprocessing
import time
import pandas as pd


sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)),'..')] + sys.path

from mintyper.validateInput import validate_input
from mintyper.kma import KMARunner
from mintyper.findBestTemplate import find_best_template_from_spa_file
import mintyper.ccphylo as ccphylo
def mintyper_pipeline(arguments):
    if not os.path.exists(arguments.output):
        os.makedirs(arguments.output)

    logging.basicConfig(
        format='%(asctime)s %(message)s',
        filename=arguments.output + '/mintyper.log',
        level=logging.INFO)
    validate_input(arguments)
    logging.info('Input validated. Mintyper is running with the following parameters:')
    logging.info(arguments)

    all_input_files_string = ' '.join(arguments.illumina + arguments.nanopore + arguments.iontorrent)

    threads = int(multiprocessing.cpu_count()/2)

    if arguments.reference != None:
        arguments.database = arguments.output + '/tmp_db'
        template_number = 1
        logging.info('Reference provided: {}'.format(arguments.reference))
        arguments.reference = concat_reference(arguments.reference, arguments.output)
        KMARunner(arguments.reference,
                  arguments.database,
                  '',
                  arguments.threads,
                  '').index()
        with open(arguments.output + '/tmp_db.name') as f:
            reference_header_text = f.readline().strip()
    else:
        KMARunner(all_input_files_string,
                  arguments.output + '/read_mapping',
                  arguments.database,
                  arguments.threads,
                  '-mem_mode -Sparse -ss c -t {}'.format(threads)).run()
        template_number, template_score, reference_header_text = find_best_template_from_spa_file(arguments.output + '/read_mapping.spa', arguments.database)
        logging.info('Best template found: {}'.format(reference_header_text))

    os.system('mkdir {}/alignments'.format(arguments.output))

    if arguments.nanopore != []:
        for item in arguments.nanopore:
            prefix = item.split('/')[-1].split('.')[0]
            KMARunner(item,
                      arguments.output + '/alignments/' + prefix,
                      arguments.database,
                      arguments.threads,
                      '-mint3 -Mt1 {} -t {} -vcf'.format(template_number, threads)).run()
    if arguments.iontorrent != []:
        for item in arguments.iontorrent:
            prefix = item.split('/')[-1].split('.')[0]
            KMARunner(item,
                      arguments.output + '/alignments/' + prefix,
                      arguments.database,
                      arguments.threads,
                      '-mint2 -Mt1 {} -t {} -vcf'.format(template_number, threads)).run()

    if arguments.illumina != []:
        for i in range(0, len(arguments.illumina), 2):
            prefix = arguments.illumina[i].split('/')[-1].split('.')[0]
            KMARunner(arguments.illumina[i] + ' ' + arguments.illumina[i+1],
                      arguments.output + '/alignments/' + prefix,
                      arguments.database,
                      arguments.threads,
                      '-mint2 -Mt1 {} -t {} -vcf'.format(template_number, threads)).run()

    time.sleep(3) #CCphylo might crash unless this. not sure why.

    ccphylo_flag = 1
    if arguments.pairwise == True:
        ccphylo_flag = 2
    if arguments.insig_prune == True:
        ccphylo_flag += 32

    ccphylo.CcphyloDist(arguments.output,
                        reference_header_text,
                        ccphylo_flag,
                        arguments.prune_distance,
                        arguments.masking_scheme).run()

    if arguments.cluster_length > 0:
        ccphylo.CcphyloDBSCAN(arguments.output, arguments.cluster_length).run()

    ccphylo.CcphyloTree(arguments.output).run()

    coord_to_filename = load_matrix_file(arguments.output + '/distmatrix.txt')
    update_variant_file(arguments.output + '/matrix_SNVs.txt', coord_to_filename)

    cmd = "cat {}/alignments/*.vcf.gz > {}/combined.vcf.gz" \
        .format(arguments.output, arguments.output)
    os.system(cmd)



def concat_reference(reference, output):
    sequence = ''
    header = '>concatenated_draft_reference'
    with open(reference, 'r') as f:
        for line in f:
            if not line.startswith('>'):
                line = line.strip()
                sequence += line

    with open('{}/{}'.format(output, 'draft_genome.fasta'), 'w') as f:
        print (header, file = f)
        print (sequence, file = f)

    return '{}/{}'.format(output, 'draft_genome.fasta')


def load_matrix_file(matrix_file):
    coord_to_filename = {}
    t = 0
    with open(matrix_file, 'r') as f:
        for line in f:
            if t != 0:  # Skip the first line
                coord = t - 1  # Adjust to start coordinates from 0
                filename = line.strip().split()[0]
                coord_to_filename[coord] = filename
            t += 1

    return coord_to_filename


def update_variant_file(variant_file, coord_to_filename):
    output = []
    with open(variant_file, 'r') as f:
        for line in f:
            line = line.strip().split('\t')
            coord_1 = int(line[0].split(', ')[0][1:])
            coord_2 = int(line[0].split(', ')[1][:-1])
            mutation = line[1]
            output_string = '({}, {})\t{}'.format(coord_to_filename[coord_1], coord_to_filename[coord_2], mutation)
            output.append(output_string)

    with open(variant_file, 'w') as f:
        for line in output:
            print(line, file=f)


