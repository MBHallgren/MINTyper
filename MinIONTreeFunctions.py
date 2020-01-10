# Copyright (c) 2019, Malte Bjørn Hallgren Technical University of Denmark
# All rights reserved.
#

#Import Libraries
import sys
import os
import argparse
import operator
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
        cmd = "kma seq2fasta -t_db " + target_dir + "temdb -seqs " + str(best_template) + " > " + target_dir + "template_sequence"
        os.system(cmd)
        print("# Mapping reads to template", file=logfile)
        return best_template, templatename
    else:
        print("# Finding best template", file=logfile)
        cmd = "kma -i " + total_filenames + " -o " + target_dir + "template_kma_results" + " -ID 50 -t_db " + kma_database_path + " -Sparse -mp 20"
        #print ("DID KMA NO SPARSE TEST")
        os.system(cmd)
        "Tried finding best template"
        try:
            infile_template = open(target_dir + "template_kma_results.spa", 'r')
            line = infile_template.readlines()[1]
            best_template = line.split("\t")[1]
            templatename = line.split("\t")[0]
            #scoring of matches, loop through and fine all maatches of 97%> BTD
            best_template_score = float(line.split("\t")[6])
            """
            if best_template_score < 90.00:
                print (
                    "None of the given templates matches any of the entries in given ref_kma_database. The input reads will now be assembled and added to the reference ref_kma_database as a new reference. After this the program will be stopped, and thus no distance matrix based analysis will be carried out.")
                print (
                "None of the given templates matches any of the entries in given ref_kma_database. The input reads will now be assembled and added to the reference ref_kma_database as a new reference. After this the program will be stopped, and thus no distance matrix based analysis will be carried out.",
                file=logfile)
                # Perform assembly based on input
                no_template_found = True
                return f
            """
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
        kma_database_path = target_dir + "temdb"

    # Illumina input
    if illumina_input != "":
        for i in range(len(illumina_input)):
            cmd = "kma -i {} -o {}{}_mapping_results -t_db {} -ref_fsa -ca -dense -cge -vcf -bc90 -Mt1 {} -t {}".format(complete_path_illumina_input[i], target_dir, illumina_input[i], kma_database_path, str(best_template), str(multi_threading))
            #cmd = "kma -i " + complete_path_illumina_input[i] + " -o " + target_dir + illumina_input[i] + "_mapping_results"  + " -t_db " + kma_database_path + " -ref_fsa -ca -dense -cge -vcf -bc90 -Mt1 " + str(best_template) + " -t " + str(multi_threading)
            os.system(cmd)
        print ("# Illumina mapping completed succesfully", file=logfile)

def illuminaMappingPE(complete_path_illumina_input, illumina_input, best_template, target_dir, kma_database_path, logfile, multi_threading, reference):
    if reference != "":
        kma_database_path = target_dir + "temdb"

    # Illumina input
    if illumina_input != "":
        for i in range(0, len(illumina_input), 2):
            cmd = "kma -ipe {} {} -o {}{}_mapping_results -t_db {} -ref_fsa -ca -dense -cge -vcf -bc90 -Mt1 {} -t {}".format(complete_path_illumina_input[i], complete_path_illumina_input[i+1], target_dir, illumina_input[i], kma_database_path, str(best_template), str(multi_threading))
            #cmd = "kma -ipe " + complete_path_illumina_input[i] + " " + complete_path_illumina_input[i + 1] + " -o " + target_dir + illumina_input[i] + "_mapping_results" + " -t_db " + kma_database_path + " -ref_fsa -ca -dense -cge -vcf -bc90 -Mt1 " + str(best_template) + " -t " + str(multi_threading)
            os.system(cmd)
        print ("# Illumina mapping completed succesfully", file=logfile)


def nanoporeMapping(complete_path_nanopore_input, nanopore_input, best_template, target_dir, kma_database_path, logfile, multi_threading, bc, reference):
    # Nanopore input
    if reference != "":
        kma_database_path = target_dir + "temdb"

    if nanopore_input != "":
        for i in range(0, len(nanopore_input)):
            cmd = "kma -i " + complete_path_nanopore_input[i] + " -o " + target_dir + nanopore_input[
                i] + "_mapping_results" + " -t_db " + kma_database_path + " -mp 20 -1t1 -dense -vcf -ref_fsa -ca -bcNano -Mt1 " + str(
                best_template) + " -t " + str(multi_threading) + " -bc " + str(bc)
            #cmd = "kma -i {} -o {}{}_mapping_results -t_db {} -mp 20 -1t1 -dense -vcf -ref_fsa -ca -bcNano -Mt1 {} -t {} -bc {}".format(complete_path_nanopore_input[i], target_dir, nanopore_input[i], kma_database_path, str(best_template), str(multi_threading), str(bc))
            #cmd = "kma -i " + complete_path_nanopore_input[i] + " -o " + target_dir + nanopore_input[i] + "_mapping_results" + " -t_db " + kma_database_path + " -mp 20 -1t1 -dense -vcf -ref_fsa -ca -bcNano -Mt1 " + str(best_template) + " -t " + str(multi_threading) + " -bc " + str(bc)
            os.system(cmd)
        print ("# Nanopore mapping completed succesfully", file=logfile)

def dcmMethylationFunction(dcmMethylation, nanopore_input, illumina_input, target_dir, logfile, dcmfilename, paired_end):
    #Non-cleaned
    dcmfile = open(dcmfilename, 'w')

    if dcmMethylation == "":
        dcmMethylationBases = "ATGC"
    elif dcmMethylation.upper() == "A" or dcmMethylation.upper() == "C" or dcmMethylation.upper() == "G" or dcmMethylation.upper() == "T":
        dcmMethylationBases = dcmMethylation.upper()
    elif dcmMethylation.upper() == "R":
        dcmMethylationBases = "AG"
    elif dcmMethylation.upper() == "Y":
        dcmMethylationBases = "CT"
    elif dcmMethylation.upper() == "M":
        dcmMethylationBases = "AC"
    elif dcmMethylation.upper() == "K":
        dcmMethylationBases = "GT"
    elif dcmMethylation.upper() == "S":
        dcmMethylationBases = "CG"
    elif dcmMethylation.upper() == "W":
        dcmMethylationBases = "AT"
    elif dcmMethylation.upper() == "H":
        dcmMethylationBases = "ACT"
    elif dcmMethylation.upper() == "B":
        dcmMethylationBases = "CGT"
    elif dcmMethylation.upper() == "V":
        dcmMethylationBases = "ACG"
    elif dcmMethylation.upper() == "D":
        dcmMethylationBases = "AGT"
    elif dcmMethylation.upper() == "N":
        dcmMethylationBases = "ATGC"


    templateseqfile = open(target_dir + "template_sequence", 'r')
    sequence = ""
    for line in templateseqfile:
        if line[0] != ">":
            line = line.rstrip()
            sequence += line
    startFrame = "CC"
    endFrame = "GG"
    motifNumber = len(dcmMethylationBases)
    possibleFrames = []
    for i in range(motifNumber):
        frame = startFrame + dcmMethylationBases[i] + endFrame
        possibleFrames.append(frame)
    motif_positions = []
    motifcount = 0
    for i in range(len(sequence)-5):
        if sequence[i:i+5] in possibleFrames:
            motif_positions.append(i)
            motifcount += 1
    print ("{} occurrences of dcmMethylation motifs found in the template, but they have NOT yet been masked".format(str(motifcount)), file = logfile)
    print (">start{}end".format("\t"), file = dcmfile)
    for i in range(len(motif_positions)):
        print (str(motif_positions[i]) + "\t" + str(motif_positions[i]+5), file = dcmfile)
    if dcmMethylation != "":
        print (str(motifcount) + " occurrences of dcmMethylation motifs found in the template has been masked!", file = logfile)
            #removing dcm-positions in mapped reads
        if nanopore_input != "":
                for i in range(len(nanopore_input)):
                    openname = target_dir + nanopore_input[i] + "_mapping_results.fsa"
                    openfile = open(openname, 'r')
                    writename = target_dir + nanopore_input[i] + "_mapping_results.fsa" + "write"
                    writefile = open(writename, 'w')
                    sequence = ""
                    for line in openfile:
                        if line[0] != '>':
                            line = line.rstrip()
                            sequence += line
                        else:
                            line = line.rstrip()
                            header = line
                            print (header, file = writefile)
                    for i in range(len(motif_positions)):
                        sequence = sequence[:int(motif_positions[i])] + "NNNNN" + sequence[int(motif_positions[i]) + 5:]
                    for i in range(len(sequence)-5):
                        if sequence[i:i+5] in possibleFrames:
                            sequence = sequence[:i] + "NNNNN" + sequence[i+5:]
                    print (sequence, file = writefile)
                    openfile.close()
                    writefile.close()
                    cmd = "rm " + openname
                    os.system(cmd)
                    cmd = "mv " + writename + " " + openname
                    os.system(cmd)
        if illumina_input != "":
            if paired_end == False:
                for i in range(len(illumina_input)):
                    openname = target_dir + illumina_input[i] + "_mapping_results.fsa"
                    openfile = open(openname, 'r')
                    writename = target_dir + illumina_input[i] + "_mapping_results.fsa" + "write"
                    writefile = open(writename, 'w')
                    sequence = ""
                    for line in openfile:
                        if line[0] != '>':
                            line = line.rstrip()
                            sequence += line
                        else:
                            line = line.rstrip()
                            header = line
                            print (header, file = writefile)
                    for i in range(len(motif_positions)):
                        sequence = sequence[:int(motif_positions[i])] + "NNNNN" + sequence[int(motif_positions[i]) + 5:]
                    for i in range(len(sequence)-5):
                        if sequence[i:i+5] in possibleFrames:
                            sequence = sequence[:i] + "NNNNN" + sequence[i+5:]
                    print (sequence, file = writefile)
                    openfile.close()
                    writefile.close()
                    cmd = "rm " + openname
                    os.system(cmd)
                    cmd = "mv " + writename + " " + openname
                    os.system(cmd)
            if paired_end == True:
                for i in range(0, len(illumina_input), 2):
                    openname = target_dir + illumina_input[i] + "_mapping_results.fsa"
                    openfile = open(openname, 'r')
                    writename = target_dir + illumina_input[i] + "_mapping_results.fsa" + "write"
                    writefile = open(writename, 'w')
                    sequence = ""
                    for line in openfile:
                        if line[0] != '>':
                            line = line.rstrip()
                            sequence += line
                        else:
                            line = line.rstrip()
                            header = line
                            print (header, file = writefile)
                    for i in range(len(motif_positions)):
                        sequence = sequence[:int(motif_positions[i])] + "NNNNN" + sequence[int(motif_positions[i]) + 5:]
                    for i in range(len(sequence)-5):
                        if sequence[i:i+5] in possibleFrames:
                            sequence = sequence[:i] + "NNNNN" + sequence[i+5:]
                    print (sequence, file = writefile)
                    openfile.close()
                    writefile.close()
                    cmd = "rm " + openname
                    os.system(cmd)
                    cmd = "mv " + writename + " " + openname
                    os.system(cmd)
    dcmfile.close()

def pruneNanopore(nanopore_input, target_dir, prune_distance):
    for i in range(len(nanopore_input)):
        openname = "{}{}_mapping_results.aln".format(target_dir, nanopore_input[i])
        template = []
        templatefile = open(openname, 'r')
        for line in templatefile:
            if line[0:8] == 'template':
                line = line.rstrip()
                line = line.split('\t')
                for t in range(len(line[1])):
                    template.append(line[1][t])
        openname = target_dir + nanopore_input[i] + "_mapping_results.fsa"
        readfile = open(openname, 'r')
        query = []
        for line in readfile:
            if line[0] != '>':
                line = line.rstrip()
                for i in range(len(line)):
                    query.append(line[i])
            else:
                line = line.rstrip()
                header = line
        # locate snps
        bad_snps = []
        snp_list = []
        for i in range(len(template)):
            if query[i].upper() != 'N' and query[i].upper() != template[i]:
                snp_list.append(i)
        # remove too close snp areas
        for i in range(len(snp_list) - 1):
            if prune_distance >= (snp_list[i + 1] - snp_list[i]):
                bad_snps.append(snp_list[i])
                bad_snps.append(snp_list[i + 1])
        cleaned_bad_snips = []
        bad_snps.sort()
        if bad_snps != []:
            for i in range(len(bad_snps) - 1):
                if bad_snps[i] != bad_snps[i + 1]:
                    cleaned_bad_snips.append(bad_snps[i])
            if bad_snps[-1] != bad_snps[-2]:
                cleaned_bad_snips.append(bad_snps[-1])
            for i in range(len(cleaned_bad_snips)):
                if prune_distance >= cleaned_bad_snips[i]:
                    for t in range(0 - cleaned_bad_snips[i], cleaned_bad_snips[i] + (prune_distance + 1),
                                   1):
                        query[cleaned_bad_snips[i] + t] = 'n'
                elif cleaned_bad_snips[i] >= (len(query) - prune_distance):
                    for t in range(-prune_distance, len(query) - cleaned_bad_snips[i], 1):
                        query[cleaned_bad_snips[i] + t] = 'n'
                else:
                    for t in range(-prune_distance, (prune_distance + 1), 1):
                        query[cleaned_bad_snips[i] + t] = 'n'
            readfile.close()
            cmd = "rm " + openname
            os.system(cmd)
            newfile = open(openname, 'w')
            print(header, file=newfile)
            sequence = ("").join(query)
            print(sequence, file=newfile)
            newfile.close()

def pruneIllumina(paired_end, illumina_input, target_dir, prune_distance):
    if illumina_input != "":
        if paired_end == False:
            for i in range(len(illumina_input)):
                openname = target_dir + illumina_input[i] + "_mapping_results.aln"
                template = []
                templatefile = open(openname, 'r')
                for line in templatefile:
                    if line[0:8] == 'template':
                        line = line.rstrip()
                        line = line.split('\t')
                        for t in range(len(line[1])):
                            template.append(line[1][t])
                openname = target_dir + illumina_input[i] + "_mapping_results.fsa"
                readfile = open(openname, 'r')
                query = []
                for line in readfile:
                    if line[0] != '>':
                        line = line.rstrip()
                        for i in range(len(line)):
                            query.append(line[i])
                    else:
                        line = line.rstrip()
                        header = line
                # locate snps
                bad_snps = []
                snp_list = []
                for i in range(len(template)):
                    if query[i].upper() != 'N' and query[i].upper() != template[i]:
                        snp_list.append(i)
                # remove too close snp areas
                for i in range(len(snp_list) - 1):
                    if prune_distance >= (snp_list[i + 1] - snp_list[i]):
                        bad_snps.append(snp_list[i])
                        bad_snps.append(snp_list[i + 1])
                cleaned_bad_snips = []
                bad_snps.sort()
                if bad_snps != []:
                    for i in range(len(bad_snps) - 1):
                        if bad_snps[i] != bad_snps[i + 1]:
                            cleaned_bad_snips.append(bad_snps[i])
                    if bad_snps[-1] != bad_snps[-2]:
                        cleaned_bad_snips.append(bad_snps[-1])
                    for i in range(len(cleaned_bad_snips)):
                        if prune_distance >= cleaned_bad_snips[i]:
                            for t in range(0 - cleaned_bad_snips[i], cleaned_bad_snips[i] + (prune_distance + 1),
                                           1):
                                query[cleaned_bad_snips[i] + t] = 'n'
                        elif cleaned_bad_snips[i] >= (len(query) - prune_distance):
                            for t in range(-prune_distance, len(query) - cleaned_bad_snips[i], 1):
                                query[cleaned_bad_snips[i] + t] = 'n'
                        else:
                            for t in range(-prune_distance, (prune_distance + 1), 1):
                                query[cleaned_bad_snips[i] + t] = 'n'
                    readfile.close()
                    cmd = "rm " + openname
                    os.system(cmd)
                    newfile = open(openname, 'w')
                    print(header, file=newfile)
                    sequence = ("").join(query)
                    print(sequence, file=newfile)
                    newfile.close()
        if paired_end == True:
            for i in range(0, len(illumina_input), 2):
                openname = target_dir + illumina_input[i] + "_mapping_results.aln"
                template = []
                templatefile = open(openname, 'r')
                for line in templatefile:
                    if line[0:8] == 'template':
                        line = line.rstrip()
                        line = line.split('\t')
                        for t in range(len(line[1])):
                            template.append(line[1][t])
                openname = target_dir + illumina_input[i] + "_mapping_results.fsa"
                readfile = open(openname, 'r')
                query = []
                for line in readfile:
                    if line[0] != '>':
                        line = line.rstrip()
                        for i in range(len(line)):
                            query.append(line[i])
                    else:
                        line = line.rstrip()
                        header = line
                # locate snps
                bad_snps = []
                snp_list = []
                for i in range(len(template)):
                    if query[i].upper() != 'N' and query[i].upper() != template[i]:
                        snp_list.append(i)
                # remove too close snp areas
                for i in range(len(snp_list) - 1):
                    if prune_distance >= (snp_list[i + 1] - snp_list[i]):
                        bad_snps.append(snp_list[i])
                        bad_snps.append(snp_list[i + 1])
                cleaned_bad_snips = []
                bad_snps.sort()
                if bad_snps != []:
                    for i in range(len(bad_snps) - 1):
                        if bad_snps[i] != bad_snps[i + 1]:
                            cleaned_bad_snips.append(bad_snps[i])
                    if bad_snps[-1] != bad_snps[-2]:
                        cleaned_bad_snips.append(bad_snps[-1])
                    for i in range(len(cleaned_bad_snips)):
                        if prune_distance >= cleaned_bad_snips[i]:
                            for t in range(0 - cleaned_bad_snips[i], cleaned_bad_snips[i] + (prune_distance + 1), 1):
                                query[cleaned_bad_snips[i] + t] = 'n'
                        elif cleaned_bad_snips[i] >= (len(query) - prune_distance):
                            for t in range(-prune_distance, len(query) - cleaned_bad_snips[i], 1):
                                query[cleaned_bad_snips[i] + t] = 'n'
                        else:
                            for t in range(-prune_distance, (prune_distance + 1), 1):
                                query[cleaned_bad_snips[i] + t] = 'n'
                    readfile.close()
                    cmd = "rm " + openname
                    os.system(cmd)
                    newfile = open(openname, 'w')
                    print(header, file=newfile)
                    sequence = ("").join(query)
                    print(sequence, file=newfile)
                    newfile.close()

def ccphylo(illumina_input, nanopore_input, paired_end, target_dir):
    nctreefile = open(target_dir + "nctree_input", 'w')
    ###Prep for nctree
    # only illumina input and no PE
    if illumina_input != "" and paired_end == False:
        for i in range(len(illumina_input)):
            openname = target_dir + illumina_input[i] + "_mapping_results.fsa"
            readfile = open(openname, 'r')
            for line in readfile:
                line = line.rstrip()
                if line[0] == ">":
                    print(">" + illumina_input[i], file=nctreefile)
                else:
                    print(line, file=nctreefile)
    # only illumina input and PE
    elif illumina_input != "" and paired_end == True:
        for i in range(0, len(illumina_input), 2):
            openname = target_dir + illumina_input[i] + "_mapping_results.fsa"
            readfile = open(openname, 'r')
            for line in readfile:
                line = line.rstrip()
                if line[0] == ">":
                    print(">" + illumina_input[i], file=nctreefile)
                else:
                    print(line, file=nctreefile)
    # only nanopore
    if nanopore_input != "":
        for i in range(len(nanopore_input)):
            openname = target_dir + nanopore_input[i] + "_mapping_results.fsa"
            readfile = open(openname, 'r')
            for line in readfile:
                line = line.rstrip()
                if line[0] == ">":
                    print(">" + nanopore_input[i], file=nctreefile)
                else:
                    print(line, file=nctreefile)
    """
    # Both illumina and nanopore input and no PE
    elif illumina_input != "" and nanopore_input != "" and paired_end == False:
        for i in range(len(illumina_input)):
            openname = target_dir + illumina_input[i] + "_mapping_results.fsa"
            readfile = open(openname, 'r')
            for line in readfile:
                line = line.rstrip()
                if line[0] == ">":
                    print(">" + illumina_input[i], file=nctreefile)
                else:
                    print(line, file=nctreefile)
        for i in range(len(nanopore_input)):
            openname = target_dir + nanopore_input[i] + "_mapping_results.fsa"
            readfile = open(openname, 'r')
            for line in readfile:
                line = line.rstrip()
                if line[0] == ">":
                    print(">" + nanopore_input[i], file=nctreefile)
                else:
                    print(line, file=nctreefile)
    # Bpth illumina and nanopore input and PE
    elif illumina_input != "" and nanopore_input != "" and paired_end == True:
        for i in range(0, len(illumina_input), 2):
            openname = target_dir + illumina_input[i] + "_mapping_results.fsa"
            readfile = open(openname, 'r')
            for line in readfile:
                line = line.rstrip()
                if line[0] == ">":
                    print(">" + illumina_input[i], file=nctreefile)
                else:
                    print(line, file=nctreefile)
        for i in range(len(nanopore_input)):
            openname = target_dir + nanopore_input[i] + "_mapping_results.fsa"
            readfile = open(openname, 'r')
            for line in readfile:
                line = line.rstrip()
                if line[0] == ">":
                    print(">" + nanopore_input[i], file=nctreefile)
                else:
                    print(line, file=nctreefile)
    """
    nctreefile.close()

def nctreeFunction(ncfile, target_dir, logfile):
    t0 = time.time()
    #
    #
    #
    etta = 0.001
    #
    # Parse command line options

    inputfile = target_dir + ncfile
    outputfile = target_dir + ncfile + "_matrix"
    difffilename = target_dir + "difffile"
    allcalled = True
    #
    # Open files
    #
    #
    # File with reads
    #
    if inputfile != None:
        inputfile = open(inputfile, "r")
    else:
        inputfile = sys.stdin
    #
    # File for general output
    #
    if outputfile != None:
        outputfile = open(outputfile, "w")
    else:
        outputfile = sys.stdout
    #
    # File for differences between samples
    #
    if difffilename != None:
        difffile = open(difffilename, "w")
    #
    # Read Input fasta file
    #
    inputseq = []
    inputseqsegments = []
    consensusseq = []
    inputname = []
    inputdesc = []
    Ninputs = 0
    i = 0
    # if inputfile != None:
    if 1 != 0:
        # allways true
        t1 = time.time()
        # sys.stdout.write("%s %d %s\n" % ("# Time used: ", int(t1-t0)," seconds"))
        sys.stdout.write("%s\n" % ("# Reading inputfile"))
        print ("%s" % ("# Nctree Data: "), file=logfile)
        for line in inputfile:
            fields = line.split()
            if len(line) > 1:
                if fields[0][0] == ">":
                    # print "# input ",line
                    if (i > 0):
                        inputseq[-1] = ''.join(inputseqsegments)
                        # print len(inputseq),len(inputseq[-1]),len(inputseq[len(inputseq)-1])
                    del inputseqsegments
                    inputseqsegments = []
                    i = 0
                    inputseq.append("")
                    consensusseq.append("")
                    inputname.append(fields[0][1:])
                    inputdesc.append(re.sub(r"^[^\s]+\s", "", line.strip()))
                else:
                    # inputseq[-1] = inputseq[-1]  + fields[0]
                    # print i
                    # print fields[0]
                    inputseqsegments.append("")
                    inputseqsegments[i] = fields[0]
                    i += 1
        inputseq[-1] = ''.join(inputseqsegments)
        # print len(inputseq),len(inputseq[-1]),len(inputseq[len(inputseq)-1])
    del inputseqsegments

    #
    # Cast sequences in numpy datastructure
    #
    nseq = len(inputseq)
    lseq0 = len(inputseq[0])
    nchrom = 1
    lseqmax = lseq0
    while (len(inputseq[nchrom]) != lseq0):
        sys.stdout.write("# Length of chromosome %d: %d\n" % (nchrom, len(inputseq[nchrom - 1])))
        print ("# Length of chromosome %d: %d" % (nchrom, len(inputseq[nchrom - 1])), file=logfile)
        if (len(inputseq[nchrom]) > lseqmax):
            lseqmax = len(inputseq[nchrom])
        nchrom += 1
    sys.stdout.write("# Length of chromosome %d: %d\n" % (nchrom, len(inputseq[nchrom - 1])))
    print ("# Length of chromosome %d: %d" % (nchrom, len(inputseq[nchrom - 1])), file=logfile)
    nstrain = int(nseq / nchrom)
    sys.stdout.write("# Number of strains: %d\n" % (nstrain))
    print ("# Number of strains: %d" % (nstrain), file=logfile)
    sys.stdout.write("# Number of chromosomes: %d\n" % (nchrom))
    print ("# Number of chromosomes: %d" % (nchrom), file=logfile)
    # for h in range(0, nchrom):
    inputseqmat = np.zeros((nseq, lseqmax), dtype=np.int8)
    non_nucleotidemat = np.zeros((nseq, lseqmax), dtype=np.int8)

    #
    # Define nucleotides as numbers
    #
    nuc2num = {
        #
        # A  adenosine          C  cytidine             G  guanine
        # T  thymidine          N  A/G/C/T (any)        U  uridine
        # K  G/T (keto)         S  G/C (strong)         Y  T/C (pyrimidine)
        # M  A/C (amino)        W  A/T (weak)           R  G/A (purine)
        # B  G/T/C              D  G/A/T                H  A/C/T
        # V  G/C/A              -  gap of indeterminate length
        #   A C G T
        # A A M R W
        # C M C S Y
        # G R S G K
        # T W Y K T
        #
        'A': 1,
        'T': 2,
        'C': 3,
        'G': 4,
        'M': 5,
        'R': 6,
        'W': 7,
        'S': 8,
        'Y': 9,
        'K': 10
    }
    #
    # Cast sequences as numpy vectors
    #
    t1 = time.time()
    sys.stdout.write("%s %d %s\n" % ("# Time used: ", int(t1 - t0), " seconds"))
    sys.stdout.write("%s\n" % ("# Cast input sequences in numpy"))
    #
    # Find non nucleotide positions
    #
    if allcalled != None:
        for i in range(0, nseq):
            for j in range(0, len(inputseq[i])):
                try:
                    # Set A, T, C, G to 1, 2, 3, 4, respectively (vector is initialized to 0)
                    inputseqmat[i][j] = nuc2num[inputseq[i][j]]
                except:
                    # Set to 1 if position do not contain an A, T, C or G (vector is initialized to 0)
                    non_nucleotidemat[i][j] = 1
    #
    # Keep only positions which are called in all sequences - this code is commented out
    #
    for n in range(0, nchrom):
        bad = 0
        sys.stdout.write("# Chromosome: %d\n" % (n + 1))
        print ("# Chromosome: %d" % (n + 1), file=logfile)
        for j in range(0, len(inputseq[n])):
            # print "pos ",j
            for l in range(0, nstrain):
                # print "strain ",l
                i = n + l * nchrom
                try:
                    # Set A, T, C, G to 1, 2, 3, 4, respectively (vector is initialized to 0)
                    inputseqmat[i][j] = nuc2num[inputseq[i][j]]
                except:
                    # Set to 1 if position do not contain an A, T, C or G (vector is initialized to 0)
                    non_nucleotidemat[i][j] = 1
                    #
                    # if options.extendtemplate != None:
                    #
                    if allcalled != None:
                        # print n,j,l, inputseq[i][j]
                        bad += 1
                        for l in range(0, nstrain):
                            i = n + l * nchrom
                            inputseqmat[i][j] = 0
                            non_nucleotidemat[i][j] = 1
                        break
        if allcalled != None:
            sys.stdout.write(
                "# Number of positions used for phylogeny in chromosome %d: %d\n" % (n + 1, len(inputseq[n]) - bad))
            print ("# Number of positions used for phylogeny in chromosome %d: %d" % (n + 1, len(inputseq[n]) - bad), file=logfile)

    # print bad

    #
    # Calculate pairwise distances
    #
    t1 = time.time()
    sys.stdout.write("%s %d %s\n" % ("# Time used: ", int(t1 - t0), " seconds"))
    print ("%s %d %s" % ("# Time used: ", int(t1 - t0), " seconds"), file=logfile)
    sys.stdout.write("%s\n" % ("# Calculating pairwise distances"))
    mat = np.zeros((nstrain, nstrain))
    prog = 0

    if difffile == None:
        #
        # fast calculation not saving different positions
        #
        for l in range(0, nstrain):
            for m in range(0, l):
                for n in range(0, nchrom):
                    i = l * nchrom + n
                    j = m * nchrom + n
                    # for j in range(0, i):
                    prog += 1
                    dist = np.sum(inputseqmat[i] != inputseqmat[j]) - np.sum(
                        non_nucleotidemat[i] != non_nucleotidemat[j])
                    # calledi = np.sum(inputseqmat[i]!=0)
                    # calledj = np.sum(inputseqmat[j]!=0)
                    # chromlength = len(inputseq[i])
                    #
                    # Estimate genetic distance if all bases had been called in both sequences
                    #
                    # scalefactor = float((chromlength*chromlength)/(calledi*calledj+0.001))
                    # print "l: ",l,"m: ",m,"n: ",n,"dist: ",dist,"mat: ",mat[l][m]
                    mat[l][m] += dist
                    # mat[l][m] += dist*scalefactor
                    mat[m][l] = mat[l][m]
                # print "l: ",l,"m: ",m,"i: ",i,"n: ",n,"j: ",j, "dist: ", dist, "sumdist: ", mat[l][m], "Calledi: ", calledi, "calledj: ",calledj, "len: ", len(inputseq[i]), "dist: ", dist, "scalefactor: ", scalefactor,"mat: ", mat[l][m]
            # sys.stdout.write("\r# %s%s done" % (int(2*100*prog/(nseq*nseq+0.001)),"%"))
            # sys.stdout.flush()
        # sys.stdout.write("\n")


    else:
        #
        # Slow version keeping information about different positions
        # NB this has to be reworked to handle multi chromosome genomes
        #
        for i in range(0, nseq):
            for j in range(0, i):
                prog += 1
                sum = 0
                for k in range(0, lseq0):
                    if (inputseqmat[i][k] != inputseqmat[j][k] and inputseqmat[i][k] > 0 and inputseqmat[j][k] > 0):
                        sum += 1
                        difffile.write("%-8s %-8s %8d %1s %1s\n" % (
                            inputname[i], inputname[j], k + 1, inputseq[i][k], inputseq[j][k]))
                mat[i][j] = sum
                mat[j][i] = mat[i][j]
                sys.stdout.write("\r# %s%s done" % (int(2 * 100 * prog / (nseq * nseq + 0.001)), "%"))
                sys.stdout.flush()
        sys.stdout.write("\n")

    #
    # Write in neighbor format
    #
    outputfile.write("%s\n" % (nstrain))
    for i in range(0, nstrain):
        mynamei = inputname[i * nchrom]
        outputfile.write("%-11s " % (mynamei))
        for j in range(0, nstrain):
            mynamej = inputname[i * nchrom]
            outputfile.write("%0.8f" % (mat[i][j]))
            if (j == nstrain - 1):
                outputfile.write("\n")
            else:
                outputfile.write(" ")
            #
    # Close files
    #
    t1 = time.time()
    sys.stdout.write("%s %d %s\n" % ("# Finishing. Time used: ", int(t1 - t0), " seconds"))
    # inputfile.close()
    # templatefile.close()
    # outputfile.close()

    #convert to new phylip format

def nctreeFormatConverter(ncfile, target_dir, logfile):
    # Nctree
    nctreeFunction(ncfile, target_dir, logfile)
    openfile = open(target_dir + ncfile + "_matrix", 'r')
    writefile = open(target_dir + ncfile + "_matrix" + "2", 'w')
    for line in openfile:
        line = line.rstrip()
        line = line.replace("\t", "!?!")
        line = line.replace(" ", "!?!")
        line = line.replace("!?!!?!!?!!?!", "\t")
        line = line.replace("!?!", "\t")
        line = line.split("\t")
        if len(line) > 1:
            for i in range(len(line)-1):
                number = line[1+i]
                number = float(number)
                number = int(round(number))
                number = str(number)
                line[1+i] = number
        line = "\t".join(line)
        print (line, file=writefile)

    cmd = "rm " + target_dir + ncfile + "_matrix"
    os.system(cmd)
    cmd = "mv " + target_dir + ncfile + "_matrix" + "2 " + target_dir + ncfile + "_matrix"
    os.system(cmd)
    openfile.close()
    writefile.close()

def ccphylo(ncfile, output_name, target_dir):
    cmd = "ccphylo tree -i " + target_dir + ncfile + "_matrix" + " -o "+ target_dir + ncfile +"_outtree"
    os.system(cmd)

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
        cmd = "rm " + target_dir + "temdb.comp.b"
        os.system(cmd)
        cmd = "rm " + target_dir + "temdb.index.b"
        os.system(cmd)
        cmd = "rm " + target_dir + "temdb.length.b"
        os.system(cmd)
        cmd = "rm " + target_dir + "temdb.name"
        os.system(cmd)
        cmd = "rm " + target_dir + "temdb.seq.b"
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

