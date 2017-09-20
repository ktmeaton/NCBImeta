# -*- coding: utf-8 -*-
"""
Created on Sat Jan 30 23:09:50 2016

@author: Poinar
"""

import os

from genomeutilities import os_check

#-----------------------------------------------------------------------#
#                            Directory Setup                            #
#-----------------------------------------------------------------------#

OS_SEP = os_check()                                                        # Retrieve the directory separator by OS

cwd = os.getcwd()
data_path = cwd + OS_SEP + "data" + OS_SEP
genome_path = cwd + OS_SEP + "genomes" + OS_SEP
genbank_path = cwd + OS_SEP + "genbank" + OS_SEP

#-----------------------------------------------------------------------#
#                             Processing                                #
#-----------------------------------------------------------------------#

#-----Iterate through sample directories-----#
for directory in os.listdir(data_path):
    org_dir = data_path + directory

    directory = directory.replace(".","_")
    directory = directory.replace("(","_")
    directory = directory.replace(")","_")
    directory = directory.replace("+","plus")
    directory = directory.replace("__","_")

    org_genome_path = genome_path + directory + ".fasta"                    # Path to the final genome fasta file
    org_genbank_path = genbank_path + directory + ".gb"


    #-----Only create file if it doesn't already exist-----#
    if os.path.exists(org_genome_path): continue

    print("Creating genome file: " + directory + ".fasta")
    org_genome_file = open(org_genome_path, "a")

    #-----Iterate through sample files-----#
    for file in os.listdir(org_dir):
        if not file.endswith(".fasta"): continue
        org_contig_file = open(org_dir + OS_SEP + file, "r")
        org_contig_seq = org_contig_file.read().replace("| ","|")               # Get rid of the pesky space after the final bar
        org_contig_seq = org_contig_seq.replace(" ","_")                         # Read in sequence and replace all spaces with underscores
        org_genome_file.write(org_contig_seq)
        org_genome_file.write("\n\n")

     #-----Only create file if it doesn't already exist-----#
    if os.path.exists(org_genbank_path): continue

    print("Creating genbank file: " + directory + ".gb")
    org_genbank_file = open(org_genbank_path, "a")

    #-----Iterate through sample files-----#
    for file in os.listdir(org_dir):
        if not file.endswith(".gb"): continue
        org_genbank_contig_file = open(org_dir + OS_SEP + file, "r")
        org_genbank_file.write(org_genbank_contig_file.read())
        org_genbank_file.write("\n\n")

    


    org_contig_file.close()
    org_genome_file.close()
    org_genbank_file.close()


#-----------------------------------------------------------------------#
#                               Cleanup                                 #
#-----------------------------------------------------------------------#
