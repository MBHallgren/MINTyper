import os
import sys
import logging
import subprocess

class KMARunner():
    def __init__(self, input, output, reference_database, argument_string):
        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(logging.INFO)
        self.check_for_kma()
        self.input = input
        self.output = output
        self.reference_database = reference_database
        self.argument_string = argument_string


    def check_for_kma(self):
        """Checks if kma is installed"""
        try:
            subprocess.call(["kma"], stdout=open(os.devnull, 'wb'))
        except FileNotFoundError:
            self.logger.info("kma is not installed correctly directly in the PATH.")
            sys.exit(1)

    def run(self):
        """runs kma"""
        kma_cmd = "kma -t_db {} -i {} -o {} {}".format(self.reference_database, self.input, self.output, self.argument_string)
        self.logger.info("Running KMA with the following command: {}".format(kma_cmd))
        os.system(kma_cmd)

    def index(self):
        """indexes kma database"""
        kma_cmd = "kma index -i {} -o {}".format(self.input, self.output)
        self.logger.info("Running KMA index with the following command: {}".format(kma_cmd))
        os.system(kma_cmd)