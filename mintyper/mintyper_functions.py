"""
Functions included in mintyper
"""

#Import Libraries
import sys
import os
import time
import subprocess

def eval_pe(mintyper_input):
    mintyper_input.i_illumina.sort()
    if mintyper_input.i_illumina[0].split("R1")[0] == mintyper_input.i_illumina[1].split("R2")[0]:
        mintyper_input.paired_end = True
        print ("# Paired-end illumina input not given but determined by the eval_pe function", file=mintyper_input.logfile)

    #hits = 0
    #for i in range(0, len(mintyper_input.i_illumina), 2):
    #    if mintyper_input.i_illumina[i].split("R1")[0] == mintyper_input.i_illumina[i+1].split("R2")[0]:
    #        hits += 1
    #rate = (hits*2)/len(mintyper_input.i_illumina)
    #if rate >= 0.5: #At least 50% are paired-end, consider all as paired-end
    #    mintyper_input.paired_end = True
    #    print ("# Paired-end illumina input not given but determined by the eval_pe function", file=mintyper_input.logfile)
    return mintyper_input

def mintyper(args):
    """
    Pipeline for mintyper. This is the main func which is called my __main__.py or mintyper_local.py
    :param args:
    :return:
    """

    args = cge_server_input(args)

    mintyper_input = MintyperHandler(args)

    if mintyper_input.cge: #Reevaluates CGE input
        joined_list = mintyper_input.i_illumina + mintyper_input.assemblies + mintyper_input.i_nanopore
        nanopore_list = []
        assembly_list = []
        illumina_list = []
        with open(mintyper_input.target_dir + "fingerReport.tsv", 'r') as infile:
            for line in infile:
                if line[0] != "#":
                    line = line.rstrip().split("\t")
                    if line[1] == "Illumina":
                        illumina_list.append(line[0])
                    elif line[1] == "Nanopore":
                        nanopore_list.append(line[0])
                    elif line[1] == "fastA":
                        assembly_list.append(line[0])
        mintyper_input.i_illumina = illumina_list
        mintyper_input.i_nanopore = nanopore_list
        mintyper_input.assemblies = assembly_list



    start_time = time.time()
    print("# Running mintyper 1.1.0 with following input conditions:", file=mintyper_input.logfile)
    print (args, file=mintyper_input.logfile)

    best_template, template_name = find_template(mintyper_input)
    mintyper_input.best_template = str(best_template)
    mintyper_input.template_name = str(template_name)

    if mintyper_input.i_illumina != []:
        mintyper_input = eval_pe(mintyper_input)

        if mintyper_input.paired_end == True:
            illumina_alignment_pe(mintyper_input)
        else:
            illumina_alignment_se(mintyper_input)
    if mintyper_input.i_nanopore != []:
        nanopore_alignment(mintyper_input)
    if mintyper_input.assemblies != []:
        print ("Assemblies are currently not supported", file=mintyper_input.logfile)
        print ("Assemblies are currently not supported")
        #assembly_alignment(mintyper_input)

    print ("calculating distance matrix")

    time.sleep(3)

    run_ccphylo(mintyper_input)

    end_time = time.time()
    delta_time = end_time - start_time

    print("mintyper total runtime: " + str(delta_time) + " seconds", file=mintyper_input.logfile)
    mintyper_input.logfile.close()
    print ("mintyper has completed")

    cleanUp(mintyper_input)

    cmd = "cat {}data_files/*.vcf.gz > {}combined.vcf.gz"\
        .format(mintyper_input.target_dir, mintyper_input.target_dir)
    os.system(cmd)

class MintyperHandler:
    """
    This class handles the input arguments and stores them.
    It also handles input/out pathing, merges input to fit KMA
    and it validates that the input is the correct format.
    """
    def __init__(self, args):
        self.i_illumina = args.i_illumina
        self.i_nanopore = args.i_nanopore
        self.paired_end = args.paired_end
        self.masking_scheme = args.masking_scheme
        self.prune_distance = args.prune_distance
        self.bc = args.bc
        self.ref_kma_database = args.ref_kma_database
        self.multi_threading = args.multi_threading
        self.reference = args.reference
        self.output_name = args.output_name
        self.exe_path = args.exe_path
        self.assemblies = args.i_assemblies
        self.insig_prune = args.insig_prune
        self.iqtree = args.iqtree
        self.fast_tree = args.fast_tree
        self.cluster_length = args.cluster_length
        self.cge = args.cge
        self._validation = self._validate_input()
        self.type = self._input_type()
        self.target_dir, self.logfile = self._handle_output()
        self.file_string = self._combine_input_files()
        self.best_template = ""
        self.tempalte_name = ""

    def _combine_input_files(self):
        self.i_illumina.sort()
        self.i_nanopore.sort()
        self.assemblies.sort()
        total_input_files = []
        if self.i_illumina != []:
            total_input_files.extend(self.i_illumina)
        if self.i_nanopore != []:
            total_input_files.extend(self.i_nanopore)
        if self.assemblies != []:
            total_input_files.extend(self.assemblies)
        total_input_files = " ".join(total_input_files)
        return total_input_files

    def _validate_input(self):
        if self.i_illumina == [] and self.i_nanopore == [] and self.assemblies == []:
            sys.exit("No input sequences was giving")

        if self.iqtree == True and self.FastTree == True:
            sys.exit("You selected both iqtree and fasttree."
                     " Please choose only one or neither for CCphylo.")

        if 3 > (len(self.i_nanopore) + len(self.i_illumina) + len(self.assemblies)):
            sys.exit("Less than 3 input files were given.")
        return True

    def _input_type(self):
        type = ""
        if self.i_illumina != []:
            type += "i"
        if self.i_nanopore != []:
            type += "n"
        if self.assemblies != []:
            type += "a"
        return type

    def _handle_output(self):
        if self.output_name[0] == "/":
            target_dir = self.output_name + "/"
        else:
            current_path = os.getcwd()
            mintyper_input.target_dir = current_path + "/" + self.output_name + "/"
        os.system("mkdir {}".format(target_dir))
        logfilename = target_dir + "logfile"
        logfile = open(logfilename, 'w')
        os.system("mkdir {}data_files".format(target_dir))
        if self.exe_path[-1] != "/":
            self.exe_path += "/"

        return target_dir, logfile

def check_draft_assembly(reference_file):
    headers = 0
    header_text = ""
    sequence = ""
    if reference_file.endswith(".gz"):
        with gzip.open(reference_file, 'rt') as reference:
            for line in reference:
                line = line.strip()
                if line.startswith(">"):
                    headers += 1
                    header_text = line  # Save header text
                else:
                    sequence += line
    else:
        with open(reference_file, 'r') as reference:
            for line in reference:
                line = line.strip()
                if line.startswith(">"):
                    headers += 1
                else:
                    sequence += line
    if headers >= 2:
        header_text = ">template_sequence"
        return True, header_text, sequence
    else:
        return False, header_text, sequence

def find_template(mintyper_input):
    if mintyper_input.reference != "":
        #Check draftassembly
        draft_genome, header_text, sequence = check_draft_assembly(mintyper_input.reference)
        if draft_genome:
            # Concatenate contigs
            with open("concatenated_genome_{}".format(mintyper_input.reference), 'w') as draft_genome_output:
                print(header, file=draft_genome_output)
                print(sequence, file=draft_genome_output)
            mintyper_input.reference = "concatenated_genome_{}".format(mintyper_input.reference)
            print("# Input: draft genome.", file=mintyper_input.logfile)
        else:
            print ("# Input reference: %s" % mintyper_input.reference, file=mintyper_input.logfile)
        print("#Making articial DB", file=mintyper_input.logfile)
        cmd = "{} index -i {}  -o {}tmp_db.ATG -Sparse ATG"\
            .format(mintyper_input.exe_path + "kma/kma",
                    mintyper_input.reference, mintyper_input.target_dir)
        os.system(cmd)
        print("# Mapping reads to template", file=mintyper_input.logfile)
        cmd = "{} -i {} -o {}template_kma_results -t_db {}tmp_db.ATG -Sparse -mp 20 -ss d".\
            format(mintyper_input.exe_path + "kma/kma", mintyper_input.file_string,
                   mintyper_input.target_dir, mintyper_input.target_dir)
        os.system(cmd)

        try:
            infile_template = open(mintyper_input.target_dir + "template_kma_results.spa", 'r')
            line = infile_template.readlines()[1]
            best_template = line.split("\t")[1]
            template_name = line.split("\t")[0]
            infile_template.close()
        except IndexError as error:
            print ("The given template does not match the given input."
                   "No analysis can be carried out.")
            sys.exit()
        print("# Best template found was " + template_name, file=mintyper_input.logfile)
        print("#Template number was: " + str(best_template), file=mintyper_input.logfile)
        cmd = "{} seq2fasta -t_db {}tmp_db.ATG -seqs {} > {}template_sequence.fasta"\
            .format(mintyper_input.exe_path + "kma/kma", mintyper_input.target_dir,
                    str(best_template), mintyper_input.target_dir)
        os.system(cmd)
        print("# Mapping reads to template", file=mintyper_input.logfile)
        return best_template, template_name
    else:
        print("# Finding best template", file=mintyper_input.logfile)
        cmd = "{} -i {} -o {}template_kma_results -ID 50 -t_db {} -Sparse -mp 20 -ss c"\
            .format(mintyper_input.exe_path + "kma/kma", mintyper_input.file_string,
                    mintyper_input.target_dir, mintyper_input.ref_kma_database)
        os.system(cmd)
        try:
            infile_template = open(mintyper_input.target_dir + "template_kma_results.spa", 'r')
            line = infile_template.readlines()[1]
            best_template = line.split("\t")[1]
            template_name = line.split("\t")[0]
            infile_template.close()
        except IndexError as error:
            sys.exit("Never found a template. Exiting. Check your SPA file.")
        print("# Best template found was " + template_name, file=mintyper_input.logfile)
        print("# Template number was: " + str(best_template), file=mintyper_input.logfile)
        cmd = "{} seq2fasta -t_db {} -seqs {} > {}template_sequence.fasta"\
            .format(mintyper_input.exe_path + "kma/kma", mintyper_input.ref_kma_database,
                    str(best_template), mintyper_input.target_dir)
        os.system(cmd)
        print("# Mapping reads to template", file=mintyper_input.logfile)
        return best_template, template_name

def assembly_alignment(mintyper_input):
    print("Assembly input was given")
    for item in mintyper_input.assemblies:
        cmd = "{} -i {} -o {}{}_alignment -ref_fsa -ca -dense -cge -vcf -bc90 -Mt1 {} -t {}"\
            .format(mintyper_input.exe_path + "kma/kma", item, mintyper_input.target_dir,
                    item.split("/")[-1], mintyper_input.best_template, mintyper_input.tempalte_name)
        if mintyper_input.reference != "":
            cmd += " -t_db {}tmp_db.ATG".format(mintyper_input.target_dir)
        else:
            cmd += " -t_db {}".format(mintyper_input.ref_kma_database)
        os.system(cmd)
    print ("# Alignment completed succesfully", file=mintyper_input.logfile)

def illumina_alignment_se(mintyper_input):
    print("Single end illumina input was given")
    for item in mintyper_input.i_illumina:
        cmd = "{} -i {} -o {}{}_alignment -ref_fsa -ca -dense -cge -vcf -bc90 -Mt1 {} -t {}"\
            .format(mintyper_input.exe_path + "kma/kma", item, mintyper_input.target_dir,
                    item.split("/")[-1], mintyper_input.best_template, mintyper_input.tempalte_name)
        if mintyper_input.reference != "":
            cmd += " -t_db {}tmp_db.ATG".format(mintyper_input.target_dir)
        else:
            cmd += " -t_db {}".format(mintyper_input.ref_kma_database)
        os.system(cmd)
    print ("# Alignment completed succesfully", file=mintyper_input.logfile)

def illumina_alignment_pe(mintyper_input):
    print("Paired end illumina input was given")
    for i in range(0, len(mintyper_input.i_illumina), 2):
        cmd = "{} -ipe {} {} -o {}{}_alignment -ref_fsa -ca -dense -cge" \
              " -vcf -bc90 -Mt1 {} -t {}"\
            .format(mintyper_input.exe_path + "kma/kma", mintyper_input.i_illumina[i],
                    mintyper_input.i_illumina[i+1],
                    mintyper_input.target_dir, mintyper_input.i_illumina[i].split("/")[-1],
                    mintyper_input.best_template, mintyper_input.tempalte_name)
        print (cmd, file=mintyper_input.logfile)
        if mintyper_input.reference != "":
            cmd += " -t_db {}tmp_db.ATG".format(mintyper_input.target_dir)
        else:
            cmd += " -t_db {}".format(mintyper_input.ref_kma_database)
        os.system(cmd)

    print ("# Alignment completed succesfully", file=mintyper_input.logfile)

def nanopore_alignment(mintyper_input):
    print("Nanopore input")

    for item in mintyper_input.i_nanopore:
        cmd = "{} -i {} -o {}{}_alignment -mp 20 -1t1 -dense -vcf -ref_fsa" \
              " -ca -bcNano -Mt1 {} -t {} -bc {}" \
            .format(mintyper_input.exe_path + "kma/kma", item, mintyper_input.target_dir,
                    item.split("/")[-1],mintyper_input.best_template,
                    mintyper_input.tempalte_name, mintyper_input.bc)
        if mintyper_input.reference != "":
            cmd += " -t_db {}tmp_db.ATG".format(mintyper_input.target_dir)
        else:
            cmd += " -t_db {}".format(mintyper_input.ref_kma_database)
        os.system(cmd)

    print ("# Alignment completed succesfully", file=mintyper_input.logfile)

def cleanUp(mintyper_input):
    if mintyper_input.reference != "":
        os.system("rm {}*tmp_db*".format(mintyper_input.target_dir))
    os.system("mv {}*alignment* {}data_files/."
              .format(mintyper_input.target_dir, mintyper_input.target_dir))


def checkOutputName(mintyper_input):
    # if used on server and output path is provided:
    if output_name[0] == "/":
        mintyper_input.target_dir = output_name + "/"
        cmd = "mkdir " + mintyper_input.target_dir
        os.system(cmd)
        cmd = "chmod 775 " + mintyper_input.target_dir
        os.system(cmd)
        cmd = "chmod 775 " + mintyper_input.target_dir + "data_files/"
        os.system(cmd)
        logfilename = mintyper_input.target_dir + "logfile"
        logfile = open(logfilename, 'w')
    else:
        current_path = os.getcwd()
        mintyper_input.target_dir = current_path + "/" + output_name + "/"
        cmd = "mkdir " + output_name
        os.system(cmd)

        logfilename = mintyper_input.target_dir + "logfile_" + output_name
        logfile = open(logfilename, 'w')

    return mintyper_input.target_dir, logfile

def run_ccphylo(mintyper_input):
    run_list = []
    fsa_list = os.listdir(mintyper_input.target_dir)
    for item in fsa_list:
        if item.endswith(".fsa"):
            if os.path.getsize(mintyper_input.target_dir + item) > 0:
                run_list.append(mintyper_input.target_dir + item)
            else:
                print ("Could not produce an alignment with {} and therefore it was excluded from the analysis".format(item))
    fsa_string = " ".join(run_list)

    if mintyper_input.iqtree:
        ccphyloflag = 1
        cmd = "{} trim --input {} --reference \"{}\" > {}trimmedalign.fsa"\
            .format(mintyper_input.exe_path + "ccphylo/ccphylo",
                    fsa_string, mintyper_input.template_name,
                    intyper_input.target_dir)
        if mintyper_input.assemblies != []:
            ccphyloflag += 8
        if mintyper_input.insig_prune == True:
            ccphyloflag += 32
        cmd = cmd + " -f {}".format(ccphyloflag)
        # Only include max flag value
        if mintyper_input.prune_distance != 0:
            cmd = cmd + " -pr {}".format(mintyper_input.prune_distance)
        if mintyper_input.masking_scheme != "":
            cmd = cmd + " -m {}".format(mintyper_input.masking_scheme)
        os.system(cmd)

        cmd = "iqtree -s {}trimmedalign.fsa > {}iqtree"\
            .format(mintyper_input.target_dir, mintyper_input.target_dir)
        os.system(cmd)
    elif mintyper_input.fast_tree:
        ccphyloflag = 1
        cmd = "{} trim --input {} --reference \"{}\" > {}trimmedalign.fsa"\
            .format(mintyper_input.exe_path + "ccphylo/ccphylo", fsa_string,
                    mintyper_input.template_name, mintyper_input.target_dir)
        if mintyper_input.assemblies != []:
            ccphyloflag += 8
        if mintyper_input.insig_prune == True:
            ccphyloflag += 32
        if mintyper_input.prune_distance != 0:
            cmd = cmd + " -pr {}".format(mintyper_input.prune_distance)
        if mintyper_input.masking_scheme != "":
            cmd = cmd + " -m {}".format(mintyper_input.masking_scheme)
        os.system(cmd)

        cmd = cmd + " -f {}".format(ccphyloflag)
        os.system(cmd)
        cmd = "FastTree -nt -gtr {}trimmedalign.fsa > {}fasttree"\
            .format(mintyper_input.target_dir, mintyper_input.target_dir)
        os.system(cmd)
    else:
        ccphyloflag = 1
        cmd = "{} dist --input {} --output {}{} --reference \"{}\"" \
              " --min_cov 1 --normalization_weight 0 2>&1" \
            .format(mintyper_input.exe_path + "ccphylo/ccphylo", fsa_string,
                    mintyper_input.target_dir, "distmatrix.txt", mintyper_input.template_name)
        if mintyper_input.assemblies != []:
            ccphyloflag = 10
        if mintyper_input.insig_prune == True:
            ccphyloflag = 32
        cmd = "{} dist --input {} --output {}{} --reference \"{}\"" \
              " --min_cov 1 --normalization_weight 0 -f {} 2>&1"\
            .format(mintyper_input.exe_path + "ccphylo/ccphylo", fsa_string,
                    mintyper_input.target_dir, "distmatrix.txt", mintyper_input.template_name, ccphyloflag)
        proc = subprocess.Popen(cmd, shell=True,
                                stdout=subprocess.PIPE, )
        output = proc.communicate()[0].decode()

        print(output, file=mintyper_input.logfile)

        time.sleep(2)
        if os.path.getsize(mintyper_input.target_dir + "distmatrix.txt") == 0:
            print("Error: Could not produce a distance matrix with ccphylo. Please check your input files. Check the logfile for Errors.")

        if mintyper_input.cluster_length > 0:
            cmd = "{} dbscan --max_distance {} --input {}{} --output {}{}"\
                .format(mintyper_input.exe_path + "ccphylo/ccphylo", mintyper_input.cluster_length,
                        mintyper_input.target_dir, "distmatrix.txt", mintyper_input.target_dir,
                        "cluster.dbscan")
            os.system(cmd)

        os.system("{} tree --input {} --output {}".format(mintyper_input.exe_path + "ccphylo/ccphylo", mintyper_input.target_dir + "distmatrix.txt", mintyper_input.target_dir + "outtree.newick"))

def cge_server_input(args):
    #Quick, bad code for easy fix. Rewrite when return from vacation.
    if args.cge:
        if args.i_illumina != []:
            illumina_files = os.listdir(args.i_illumina[0])
            illumina_list = []
            for item in illumina_files:
                illumina_list.append(args.i_illumina[0] + item)
            args.i_illumina = illumina_list
        if args.i_nanopore != []:
            nanopore_files = os.listdir(args.i_nanopore[0])
            nanopore_list = []
            for item in nanopore_files:
                nanopore_list.append(args.i_nanopore[0] + item)
            args.i_nanopore = nanopore_list
        if args.i_assemblies != []:
            assemblies_files = os.listdir(args.i_assemblies[0])
            assemblies_list = []
            for item in assemblies_files:
                assemblies_list.append(args.i_assemblies[0] + item)
            args.i_assemblies = assemblies_list
    return args

