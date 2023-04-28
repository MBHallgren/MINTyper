import os
import sys
import logging

sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)),'..')] + sys.path

def validate_input(arguments):
    check_illumina_pe(arguments.illumina)
    check_database(arguments.database, arguments.reference)

def check_illumina_pe(files):
    if len(files) % 2 != 0:
        logging.info('Error: Illumina files must be paired-end. Exiting. An odd number of files was provided.')
        sys.exit('Error: Illumina files must be paired-end. Exiting.')

def check_database(database, reference):
    if reference == None:
        db_list = ['.name', '.comp.b', '.length.b', '.seq.b']
        for db in db_list:
            if not os.path.isfile(database + db):
                logging.info('Error: KMA database not found. Exiting.')
                sys.exit('Error: KMA database not found. Exiting.')

def check_for_ccphylo():
    """Checks if ccphylo is installed"""
    try:
        subprocess.call(["ccphylo"], stdout=open(os.devnull, 'wb'))
    except FileNotFoundError:
        self.logger.info("ccphylo is not installed correctly directly in the PATH.")
        sys.exit(1)

def check_for_kma():
    """Checks if kma is installed"""
    try:
        subprocess.call(["kma"], stdout=open(os.devnull, 'wb'))
    except FileNotFoundError:
        self.logger.info("kma is not installed correctly directly in the PATH.")
        sys.exit(1)
