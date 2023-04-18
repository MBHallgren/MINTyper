import os
import sys
import logging

sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)),'..')] + sys.path

def validate_input(arguments):
    check_illumina_pe(arguments.illumina)
    check_database(arguments.database)

def check_illumina_pe(files):
    if len(files) % 2 != 0:
        logging.info('Error: Illumina files must be paired-end. Exiting.')
        sys.exit('Error: Illumina files must be paired-end. Exiting.')

def check_database(database):
    db_list = ['.name', '.comp.b', '.length.b', '.seq.b']
    for db in db_list:
        if not os.path.isfile(database + db):
            logging.info('Error: KMA database not found. Exiting.')
            sys.exit('Error: KMA database not found. Exiting.')
