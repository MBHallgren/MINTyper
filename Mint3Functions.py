# Copyright (c) 2019, Malte Bjørn Hallgren Technical University of Denmark
# All rights reserved.
#

#Import Libraries
import sys
import os
import argparse
import time
import gc
import numpy as np
import array
from optparse import OptionParser
from operator import itemgetter
import re
import json
import sqlite3
import json

def findTemplateResearch(total_filenames, target_dir, kma_database_path, logfile, reference):
    #NON-LEAN
    no_template_found = False
    best_template = ""
    best_template_score = 0.0
    templatename = ""
    if reference != "":
        best_template = reference
        print("# The reference given by the user was: " + best_template, file=logfile)
        print("#Making articial DB", file=logfile)
        cmd = "kma_index -i " + best_template + " -o " + target_dir + "temdb.ATG -Sparse ATG"
        os.system(cmd)
        print("# Mapping reads to template", file=logfile)
        cmd = "kma -i " + total_filenames + " -o " + target_dir + "template_kma_results" + " -t_db " + target_dir + "temdb.ATG" + " -Sparse -mp 20"
        os.system(cmd)

        print ("midway complete")
        try:
            infile_template = open(target_dir + "template_kma_results.spa", 'r')
            line = infile_template.readlines()[1]
            best_template = line.split("\t")[1]
            templatename = line.split("\t")[0]
            infile_template.close()
        except IndexError as error:
            print ("The template you have stated as a reference does not match the input reads to a good enough degree to make any type of analysis.")
            sys.exit()
        print("# Best template found was " + templatename, file=logfile)
        print("#Template number was: " + str(best_template), file=logfile)
        cmd = "kma seq2fasta -t_db " + target_dir + "temdb.ATG -seqs " + str(best_template) + " > " + target_dir + "template_sequence"
        os.system(cmd)
        print("# Mapping reads to template", file=logfile)
        return best_template, templatename
    else:
        print("# Finding best template", file=logfile)
        cmd = "kma -i " + total_filenames + " -o " + target_dir + "template_kma_results" + " -ID 50 -t_db " + kma_database_path + " -Sparse -mp 20"
        os.system(cmd)
        try:
            infile_template = open(target_dir + "template_kma_results.spa", 'r')
            line = infile_template.readlines()[1]
            best_template = line.split("\t")[1]
            templatename = line.split("\t")[0]
            infile_template.close()
        except IndexError as error:
            sys.exit("Never found a template. Exiting. Check your SPA file.")
        print("# Best template found was " + templatename, file=logfile)
        print("# Template number was: " + str(best_template), file=logfile)
        cmd = "kma seq2fasta -t_db " + kma_database_path + " -seqs " + str(best_template) + " > " + target_dir + "template_sequence"
        os.system(cmd)
        print("# Mapping reads to template", file=logfile)
        return best_template, templatename

#Mapping
def illuminaMappingForward(complete_path_illumina_input, illumina_input, best_template, target_dir, kma_database_path, logfile, multi_threading, reference):
    if reference != "":
        kma_database_path = target_dir + "temdb.ATG"

    # Illumina input
    if illumina_input != "":
        for i in range(len(illumina_input)):
            cmd = "kma -i {} -o {}{}_mapping_results -t_db {} -ref_fsa -ca -dense -cge -vcf -bc90 -Mt1 {} -t {}".format(complete_path_illumina_input[i], target_dir, illumina_input[i], kma_database_path, str(best_template), str(multi_threading))
            os.system(cmd)
        print ("# Illumina mapping completed succesfully", file=logfile)

def illuminaMappingPE(complete_path_illumina_input, illumina_input, best_template, target_dir, kma_database_path, logfile, multi_threading, reference):
    if reference != "":
        kma_database_path = target_dir + "temdb.ATG"

    # Illumina input
    if illumina_input != "":
        for i in range(0, len(illumina_input), 2):
            cmd = "kma -ipe {} {} -o {}{}_mapping_results -t_db {} -ref_fsa -ca -dense -cge -vcf -bc90 -Mt1 {} -t {}".format(complete_path_illumina_input[i], complete_path_illumina_input[i+1], target_dir, illumina_input[i], kma_database_path, str(best_template), str(multi_threading))
            os.system(cmd)
        print ("# Illumina mapping completed succesfully", file=logfile)


def nanoporeMapping(complete_path_nanopore_input, nanopore_input, best_template, target_dir, kma_database_path, logfile, multi_threading, bc, reference):
    # Nanopore input
    if reference != "":
        kma_database_path = target_dir + "temdb.ATG"

    if nanopore_input != "":
        for i in range(0, len(nanopore_input)):
            cmd = "kma -i " + complete_path_nanopore_input[i] + " -o " + target_dir + nanopore_input[
                i] + "_mapping_results" + " -t_db " + kma_database_path + " -mp 20 -1t1 -dense -vcf -ref_fsa -ca -bcNano -Mt1 " + str(
                best_template) + " -t " + str(multi_threading) + " -bc " + str(bc)
            os.system(cmd)
        print ("# Nanopore mapping completed succesfully", file=logfile)

def cleanUp( target_dir, illumina_input, nanopore_input, paired_end, reference):
    save_files_bool = True
    if save_files_bool == False:
        if illumina_input != "" and paired_end == False:
            for i in range(len(illumina_input)):
                cmd = "rm " + target_dir + illumina_input[i] + "_mapping_results.aln"
                os.system(cmd)
                cmd = "rm " + target_dir + illumina_input[i] + "_mapping_results.frag.gz"
                os.system(cmd)
                cmd = "rm " + target_dir + illumina_input[i] + "_mapping_results.res"
                os.system(cmd)
                cmd = "rm " + target_dir + illumina_input[i] + "_mapping_results.fsa"
                os.system(cmd)
                cmd = "mv " + target_dir + illumina_input[i] + "_mapping_results.vcf.gz " + target_dir + "DataFiles"
                os.system(cmd)
        elif illumina_input != "" and paired_end == True:
            for i in range(0, len(illumina_input), 2):
                cmd = "rm " + target_dir + illumina_input[i] + "_mapping_results.aln"
                os.system(cmd)
                cmd = "rm " + target_dir + illumina_input[i] + "_mapping_results.frag.gz"
                os.system(cmd)
                cmd = "rm " + target_dir + illumina_input[i] + "_mapping_results.res"
                os.system(cmd)
                cmd = "rm " + target_dir + illumina_input[i] + "_mapping_results.fsa"
                os.system(cmd)
                cmd = "mv " + target_dir + illumina_input[i] + "_mapping_results.vcf.gz " + target_dir + "DataFiles"
                os.system(cmd)
        if nanopore_input != "":
            for i in range(len(nanopore_input)):
                cmd = "rm " + target_dir + nanopore_input[i] + "_mapping_results.aln"
                os.system(cmd)
                cmd = "rm " + target_dir + nanopore_input[i] + "_mapping_results.frag.gz"
                os.system(cmd)
                cmd = "rm " + target_dir + nanopore_input[i] + "_mapping_results.res"
                os.system(cmd)
                cmd = "rm " + target_dir + nanopore_input[i] + "_mapping_results.fsa"
                os.system(cmd)
                cmd = "mv " + target_dir + nanopore_input[i] + "_mapping_results.vcf.gz " + target_dir + "DataFiles"
                os.system(cmd)
    elif save_files_bool == True:
        if illumina_input != "" and paired_end == False:
            for i in range(len(illumina_input)):
                cmd = "mv " + target_dir + illumina_input[i] + "_mapping_results.aln" + " " + target_dir +"DataFiles"
                os.system(cmd)
                cmd = "mv " + target_dir + illumina_input[i] + "_mapping_results.frag.gz" + " " + target_dir +"DataFiles"
                os.system(cmd)
                cmd = "mv " + target_dir + illumina_input[i] + "_mapping_results.res" + " " + target_dir +"DataFiles"
                os.system(cmd)
                cmd = "mv " + target_dir + illumina_input[i] + "_mapping_results.fsa" + " " + target_dir +"DataFiles"
                os.system(cmd)
                cmd = "mv " + target_dir + illumina_input[i] + "_mapping_results.vcf.gz " + target_dir + "DataFiles"
                os.system(cmd)
        elif illumina_input != "" and paired_end == True:
            for i in range(0, len(illumina_input), 2):
                cmd = "mv " + target_dir + illumina_input[i] + "_mapping_results.aln" + " " + target_dir +"DataFiles"
                os.system(cmd)
                cmd = "mv " + target_dir + illumina_input[i] + "_mapping_results.frag.gz" + " " + target_dir +"DataFiles"
                os.system(cmd)
                cmd = "mv " + target_dir + illumina_input[i] + "_mapping_results.res" + " " + target_dir +"DataFiles"
                os.system(cmd)
                cmd = "mv " + target_dir + illumina_input[i] + "_mapping_results.fsa" + " " + target_dir +"DataFiles"
                os.system(cmd)
                cmd = "mv " + target_dir + illumina_input[i] + "_mapping_results.vcf.gz " + target_dir + "DataFiles"
                os.system(cmd)
        if nanopore_input != "":
            for i in range(len(nanopore_input)):
                cmd = "mv " + target_dir + nanopore_input[i] + "_mapping_results.aln" + " " + target_dir +"DataFiles"
                os.system(cmd)
                cmd = "mv " + target_dir + nanopore_input[i] + "_mapping_results.frag.gz" + " " + target_dir +"DataFiles"
                os.system(cmd)
                cmd = "mv " + target_dir + nanopore_input[i] + "_mapping_results.res" + " " + target_dir +"DataFiles"
                os.system(cmd)
                cmd = "mv " + target_dir + nanopore_input[i] + "_mapping_results.fsa" + " " + target_dir +"DataFiles"
                os.system(cmd)
                cmd = "mv " + target_dir + nanopore_input[i] + "_mapping_results.vcf.gz " + target_dir + "DataFiles"
                os.system(cmd)

    if reference != "":
        cmd = "rm " + target_dir + "temdb.ATG.comp.b"
        os.system(cmd)
        cmd = "rm " + target_dir + "temdb.ATG.index.b"
        os.system(cmd)
        cmd = "rm " + target_dir + "temdb.ATG.length.b"
        os.system(cmd)
        cmd = "rm " + target_dir + "temdb.ATG.name"
        os.system(cmd)
        cmd = "rm " + target_dir + "temdb.ATG.seq.b"
        os.system(cmd)

def load_illumina(illumina_path_input):
    if illumina_path_input != "":
        path = illumina_path_input
        illumina_files = os.listdir(path)
        illumina_files.sort()
    else:
        illumina_files = ""
    return illumina_files

def load_nanopore(nanopore_path_input):
    if nanopore_path_input != "":
        path = nanopore_path_input
        nanopore_files = os.listdir(path)
        nanopore_files.sort()
    else:
        nanopore_files = ""
    return nanopore_files

def generate_complete_path_illumina_files(illumina_files, illumina_path_input):
    path = illumina_path_input
    complete_path_illumina_files = []
    for i in range(len(illumina_files)):
        complete_path_illumina_files.append(path + illumina_files[i])
    return complete_path_illumina_files

def generate_complete_path_nanopore_files(nanopore_files, nanopore_path_input):
    path = nanopore_path_input
    complete_path_nanopore_files = []
    for i in range(len(nanopore_files)):
        complete_path_nanopore_files.append(path + nanopore_files[i])
    return complete_path_nanopore_files

def combine_input_files(illumina_files, nanopore_files):
    if illumina_files == "":
        total_input_files = nanopore_files
    elif nanopore_files == "":
        total_input_files = illumina_files
    else:
        total_input_files = illumina_files + nanopore_files
    total_input_files = " ".join(total_input_files)
    return total_input_files

def logfileConditionsResearch(logfile, dcmMethylation, prune, prune_distance, bc, ref_kma_database, multi_threading, reference, output_name, paired_end):
    logdict = {}
    if dcmMethylation != "":
        logdict['dcmMethylation'] = dcmMethylation
    else:
        logdict['dcmMethylation'] = ""
    if prune != False:
        logdict['prune'] = str(prune)
    else:
        logdict['prune'] = str(False)
    if prune_distance != 10:
        logdict['prune_distance'] = prune_distance
    else:
        logdict['prune_distance'] = 10
    if bc != 0.7:
        logdict['bc'] = bc
    else:
        logdict['bc'] = 0.7
    if ref_kma_database != "":
        logdict['ref_kma_database'] = ref_kma_database
    else:
        logdict['ref_kma_database'] = ""
    if multi_threading != 1:
        logdict['multi_threading'] = multi_threading
    else:
        logdict['multi_threading'] = 1
    if reference != "":
        logdict['reference'] = reference
    else:
        logdict['reference'] = ""
    if output_name != "":
        logdict['output_name'] = output_name
    else:
        logdict['output_name'] = ""
    if paired_end != False:
        logdict['paired_end'] = str(paired_end)
    else:
        logdict['paired_end'] = str(False)
    print (logdict, file=logfile)

def varriansfileRenamer(target_dir):
    varfile =  open("{}nucleotideVarriance_file".format(target_dir),'w')
    readfile = open("{}nucleotideVarriance".format(target_dir), 'r')
    distfile = open("{}distmatrix".format(target_dir),'r')
    indexnames = {}
    for line in distfile:
        line = line.rstrip()
    return True



def mutationSpotter(target_dir):
    templatefile =  open("{}template_sequence".format(target_dir),'r')
    sequence = ""

    for line in templatefile:
        line = line.rstrip()
        if line[0] != ">":
            sequnce += line
    print ("Template lenght is {}".format(str(len(sequence))))
    return True
