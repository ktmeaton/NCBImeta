# -*- coding: utf-8 -*-
"""
Created on Sat Jan 30 23:09:50 2016

@author: Poinar
"""

import os  


cwd = os.getcwd()
data_path = cwd + "\\data\\"
genome_path = cwd + "\\genomes\\"
dir_list = []
title_list = ['Yersinia_pestis_CO92_(enterobacteria)_GCF_000009065.1']

for roots, dirs, filenames in os.walk(data_path):
    if dirs:
        organism_id = dirs[0]
        organism_directory = data_path + organism_id
        organism_genome = open(genome_path + organism_id + ".fasta", "w")
        for file in os.listdir(organism_directory):
            print file
            organism_contig = open(organism_directory + "\\" + file, "r")
            organism_sequence = organism_contig.read()
            organism_genome.write(organism_sequence)
            organism_contig.close()
