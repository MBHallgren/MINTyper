import os
import sys
import logging
import subprocess

class CcphyloDist():
    def __init__(self, target_dir, reference_header_text, ccphylo_flag):
        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(logging.INFO)
        self.check_for_ccphylo()
        self.target_dir = target_dir
        self.reference_header_text = reference_header_text
        self.ccphylo_flag = ccphylo_flag
        self.prepare_list_of_alignment_files()

    def prepare_list_of_alignment_files(self):
        """Returns a list of alignment files"""
        run_list = []
        fsa_list = os.listdir(self.target_dir + '/alignments/')
        for item in fsa_list:
            if item.endswith(".fsa"):
                if os.path.getsize(self.target_dir + '/alignments/' + item) > 0:  # non empty alignment
                    run_list.append(self.target_dir + '/alignments/' + item)
                else:
                    logging.info(
                        'The alignment file {} is empty and therefore it was excluded from the analysis'.format(item))
        self.alignment_string = " ".join(run_list)

    def check_for_ccphylo(self):
        """Checks if ccphylo is installed"""
        try:
            subprocess.call(["ccphylo"], stdout=open(os.devnull, 'wb'))
        except FileNotFoundError:
            self.logger.info("ccphylo is not installed correctly directly in the PATH.")
            sys.exit(1)

    def run(self):
        """runs ccphylo"""
        cmd = "ccphylo dist --input {} --output {}/distmatrix.txt --reference \"{}\"" \
              " --min_cov 1 --normalization_weight 0" \
            .format(self.alignment_string, self.target_dir, self.reference_header_text)
        if self.ccphylo_flag != 1:
            cmd += " -f {} 2>&1".format(self.ccphylo_flag)
        else:
            cmd += " 2>&1"
        self.logger.info("Running ccphylo with the following command: {}".format(cmd))
        os.system(cmd)


class CcphyloTree():
    def __init__(self, target_dir):
        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(logging.INFO)
        self.check_for_ccphylo()
        self.target_dir = target_dir
        self.prepare_list_of_alignment_files()

    def check_for_ccphylo(self):
        """Checks if ccphylo is installed"""
        try:
            subprocess.call(["ccphylo"], stdout=open(os.devnull, 'wb'))
        except FileNotFoundError:
            self.logger.info("ccphylo is not installed correctly directly in the PATH.")
            sys.exit(1)

    def prepare_list_of_alignment_files(self):
        """Returns a list of alignment files"""
        run_list = []
        fsa_list = os.listdir(self.target_dir + '/alignments/')
        for item in fsa_list:
            if item.endswith(".fsa"):
                if os.path.getsize(self.target_dir + '/alignments/' + item) > 0:  # non empty alignment
                    run_list.append(self.target_dir + '/alignments/' + item)
                else:
                    logging.info(
                        'The alignment file {} is empty and therefore it was excluded from the analysis'.format(item))
        self.alignment_string = " ".join(run_list)

    def run(self):
        """runs ccphylo tree"""
        cmd = "ccphylo tree --input {}/distmatrix.txt --output {}/tree.newick" \
            .format(self.target_dir, self.target_dir)
        self.logger.info("Running ccphylo with the following command: {}".format(cmd))
        os.system(cmd)

class CcphyloDBSCAN():
    def __init__(self, target_dir, cluster_length):
        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(logging.INFO)
        self.check_for_ccphylo()
        self.target_dir = target_dir
        self.cluster_length = cluster_length
        self.prepare_list_of_alignment_files()

    def check_for_ccphylo(self):
        """Checks if ccphylo is installed"""
        try:
            subprocess.call(["ccphylo"], stdout=open(os.devnull, 'wb'))
        except FileNotFoundError:
            self.logger.info("ccphylo is not installed correctly directly in the PATH.")
            sys.exit(1)

    def prepare_list_of_alignment_files(self):
        """Returns a list of alignment files"""
        run_list = []
        fsa_list = os.listdir(self.target_dir + '/alignments/')
        for item in fsa_list:
            if item.endswith(".fsa"):
                if os.path.getsize(self.target_dir + '/alignments/' + item) > 0:  # non empty alignment
                    run_list.append(self.target_dir + '/alignments/' + item)
                else:
                    logging.info(
                        'The alignment file {} is empty and therefore it was excluded from the analysis'.format(item))
        self.alignment_string = " ".join(run_list)

    def run(self):
        """runs ccphylo dbscan"""
        cmd = "ccphylo dbscan --max_distance {} --input {}{} --output {}{}" \
            .format(cluster_length, self.target_dir, '/distmatrix.txt', self.target_dir, '/clusters.txt')
        self.logger.info("Running ccphylo with the following command: {}".format(cmd))
        os.system(cmd)



