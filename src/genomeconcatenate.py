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

#-----------------------------------------------------------------------#
#                             Processing                                #
#-----------------------------------------------------------------------#

#-----Iterate through sample directories-----#
for directory in os.listdir(data_path):
    org_dir = data_path + directory
    org_genome_path = genome_path + directory + ".fasta"                    # Path to the final genome fasta file

    #-----Only create file if it doesn't already exist-----#
    if os.path.exists(org_genome_path): continue

    print("Creating genome file: " + directory + ".fasta")
    org_genome_file = open(org_genome_path, "w")

    #-----Iterate through sample files-----#
    for file in os.listdir(org_dir):
        org_contig_file = open(org_dir + OS_SEP + file, "r")
        org_contig_seq = org_contig_file.read().replace("| ","|")               # Get rid of the pesky space after the final bar
        org_contig_seq = org_contig_seq.replace(" ","_")                         # Read in sequence and replace all spaces with underscores
        org_genome_file.write(org_contig_seq)


    org_contig_file.close()
    org_genome_file.close()


#-----------------------------------------------------------------------#
#                               Cleanup                                 #
#-----------------------------------------------------------------------#
