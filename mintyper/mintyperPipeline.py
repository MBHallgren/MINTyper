import os
import sys
import logging

sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)),'..')] + sys.path

from mintyper.validateInput import validate_input
from mintyper.kma import KMARunner
from mintyper.findBestTemplate import find_best_template

def mintyper_pipline(arguments):
    logging.basicConfig(
        format='%(asctime)s %(message)s',
        filename=arguments.output + '/mintyper.log',
        level=logging.INFO)
    validate_input(arguments)
    logging.info('Input validated. Mintyper is running with the following parameters:')
    logging.info(arguments)

    all_input_files_string = ' '.join(arguments.illumina + arguments.nanopore + arguments.iontorrent)

    if args.reference != None:
        args.database = arguments.output + '/tmp_db'
        KMARunner(args.reference,
                  args.database,
                  '',
                  '').index()
    else:
        KMARunner(all_input_files_string,
                  arguments.output + '/read_mapping',
                  arguments.database,
                  'mem_mode -Sparse -ss c').run()
        template_number, template_score, reference_header_text, find_best_template(arguments.output + '/read_mapping.res', arguments.database)
        logging.info('Best template found: {}'.format(reference_header_text))

    sys.exit()


    if arguments.illumina != []:
        KMARunner(arguments.illumina,
                  arguments.output + '/illumina_alignment',
                  arguments.database,
                  arguments.kma_arguments).run()

