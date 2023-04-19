import os
import sys
import logging
import multiprocessing

sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)),'..')] + sys.path

from mintyper.validateInput import validate_input
from mintyper.kma import KMARunner
from mintyper.findBestTemplate import find_best_template_from_spa_file
from mintyper.ccphylo import CcphyloRunner

def mintyper_pipline(arguments):
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
    print (all_input_files_string)

    if arguments.reference != None:
        arguments.database = arguments.output + '/tmp_db'
        template_number = 1
        info.logging('Reference provided: {}'.format(arguments.reference))
        KMARunner(arguments.reference,
                  arguments.database,
                  '',
                  '').index()
        with open(arguments.output + '/tmp_db.name') as f:
            reference_header_text = f.readline()
    else:
        KMARunner(all_input_files_string,
                  arguments.output + '/read_mapping',
                  arguments.database,
                  '-mem_mode -Sparse -ss c').run()
        template_number, template_score, reference_header_text = find_best_template_from_spa_file(arguments.output + '/read_mapping.spa', arguments.database)
        logging.info('Best template found: {}'.format(reference_header_text))

    threads = int(multiprocessing.cpu_count()/2)

    os.system('mkdir {}/alignments'.format(arguments.output))

    if arguments.nanopore != []:
        for item in arguments.nanopore:
            prefix = item.split('/')[-1].split('.')[0]
            KMARunner(item,
                      arguments.output + '/alignments/' + prefix,
                      arguments.database,
                      '-mint3 -Mt1 {} -t {}'.format(template_number, threads)).run()
    if arguments.iontorrent != []:
        for item in arguments.iontorrent:
            prefix = item.split('/')[-1].split('.')[0]
            KMARunner(item,
                      arguments.output + '/alignments/' + prefix,
                      arguments.database,
                      '-mint2 -Mt1 {} -t {}'.format(template_number, threads)).run()

    if arguments.illumina != []:
        for i in range(0, len(arguments.illumina), 2):
            prefix = arguments.illumina[i].split('/')[-1].split('.')[0]
            KMARunner(arguments.illumina[i] + ' ' + arguments.illumina[i+1],
                      arguments.output + '/alignments/' + prefix,
                      arguments.database,
                      '-mint2 -Mt1 {} -t {}'.format(template_number, threads)).run()

    if arguments.assemblies != []:
        for item in arguments.assemblies:
            prefix = item.split('/')[-1].split('.')[0]
            KMARunner(item,
                      arguments.output + '/alignments/' + prefix,
                      arguments.database,
                      '-mint1 -Mt1 {} -t {}'.format(template_number, threads)).run()

    time.sleep(3) #CCphylo might crash unless this. not sure why.

    ccphylo_flag = 1
    if mintyper_input.assemblies != []:
        ccphyloflag = 10
    if mintyper_input.insig_prune == True:
        ccphyloflag = 32

    CcphyloRunner(arguments.output,
                  reference_header_text,
                  ccphylo_flag,
                  arguments.cluster_length).dist()

    if arguments.cluster_length > 0:
        CcphyloRunner(arguments.output,
                      reference_header_text,
                      ccphylo_flag,
                      arguments.cluster_length).dbscan()

    CcphyloRunner(arguments.output,
                  reference_header_text,
                  ccphylo_flag,
                  arguments.cluster_length).tree()





