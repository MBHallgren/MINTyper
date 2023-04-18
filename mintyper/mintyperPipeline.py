import os
import sys
import logging

sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)),'..')] + sys.path

from mintyper.validateInput import validate_input
from mintyper.kma import KMARunner
from mintyper.findBestTemplate import find_best_template_from_spa_file

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
        KMARunner(arguments.reference,
                  arguments.database,
                  '',
                  '').index()
    else:
        KMARunner(all_input_files_string,
                  arguments.output + '/read_mapping',
                  arguments.database,
                  '-mem_mode -Sparse -ss c').run()
        template_number, template_score, reference_header_text = find_best_template_from_spa_file(arguments.output + '/read_mapping.spa', arguments.database)
        logging.info('Best template found: {}'.format(reference_header_text))

    sys.exit()


    if arguments.illumina != []:
        KMARunner(arguments.illumina,
                  arguments.output + '/illumina_alignment',
                  arguments.database,
                  arguments.kma_arguments).run()

