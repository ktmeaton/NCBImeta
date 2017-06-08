# -*- coding: utf-8 -*-
"""
Created on Sat Jan 30 14:09:59 2016
Edited: June 6, 2017 - post-GI number removal check

@author: Katherine Eaton

"Genome Collector Module"
"""

import argparse                                                                # Command-line argument parsing
import sqlite3
import os

from genomeerrors import *
import genomeassembly
import genomeutilities
import genomebioproject
import genomebiosample
import genomenucleotide


#-----------------------------------------------------------------------#
#                            Argument Parsing                           #
#-----------------------------------------------------------------------#

genomeutilities.check_accessory_dir()                                                           # Create accessory directories

parser = argparse.ArgumentParser(description='Description of the Genome Collector.',
                                 add_help=True,
                                 version = 'GenomeCollector 2.0')

mandatory = parser.add_argument_group('mandatory')
create_update = parser.add_argument_group('create/update')
bonus = parser.add_argument_group('bonus')


mandatory.add_argument('--mode',
                    help = 'Create, Update, or Delete Genome Database',
                    action = 'store',
                    dest = 'mode',
                    required = True)


mandatory.add_argument('--database',
                    help='Genome Database Name',
                    type = str,
                    action = 'store',
                    dest = 'dbname',
                    required=True)

bonus.add_argument('--assembly-status',
                   help='Filter by Assembly Status: Complete, Chromosome, Scaffold, or Contig',
                   type = str,
                   action = 'store',
                   dest = 'assembly_status',
                   required = False)

args = vars(parser.parse_args())


MODE = args['mode']
dbName = args['dbname'] + ".sqlite"
ASSEMBLY_STATUS = args['assembly_status']


#-----------------------------------------------------------------------#
#                           Argument Checking                           #
#-----------------------------------------------------------------------#


#---------------------------Delete Database-----------------------------#
if MODE.lower() == 'delete':                                                   # Delete mode check
    if os.path.exists(dbName):                                                 # Make sure database exists
        os.remove(dbName)                                                      # Delete database
        print('\nDeleting database: ' + dbName)
    else:
        raise ErrorDBNotExists(dbName)                                         # If database doesn't exists, raise error


#---------------------------Create Database-----------------------------#
elif MODE.lower() == 'create':                                                 # Create mode check
    if not os.path.exists(dbName):
        conn = sqlite3.connect(dbName)
        conn.commit()
        print('\nCreating database: ' + dbName)
    else:
        raise ErrorDBExists(dbName)


#---------------------------Update Database-----------------------------#
elif MODE.lower() == 'update':
    if os.path.exists(dbName):
        conn = sqlite3.connect(dbName)
        print('\nOpening database: ' + dbName)
    else:
        raise ErrorDBNotExists(dbName)

#----------------------------Invalid Mode------------------------------#
else:
    raise ErrorInvalidMode(MODE)


#------------------------Create Mode Processing------------------------#

if MODE.lower() == 'create' or MODE.lower() == 'update':

   #-----------------------------NCBI Info-----------------------------#
    EMAIL =  raw_input('\n\
    Please enter a valid email address for NCBI queries: ')
    ORGANISM = raw_input('\n\
    Please enter the name of the organism of interest: ')

    #------------------------Assembly User Input-----------------------#
    modify_assembly = raw_input('\n\
    Do you want to create/update the Assembly table? (Y/N)')
    while modify_assembly.lower() != 'y' and modify_assembly.lower() != 'n':
        print ("\tInvalid input, must be either 'Y','y','N', or 'n'.")
        modify_assembly = raw_input('\n\
        Do you want to create/update the Assembly table? (Y/N)')

    #---------------------------Assembly Table--------------------------#

    if modify_assembly.lower() == 'y':
        genomeassembly.AssemblyTable(dbName, ORGANISM, EMAIL)

    #------------------------Bioproject User Input----------------------#
    modify_bioproject = raw_input('\n\
    Do you want to create/update the BioProject table? (Y/N)')
    while modify_bioproject.lower() != 'y' and modify_bioproject.lower() != 'n':
        print ("\tInvalid input, must be either 'Y','y','N', or 'n'.")
        modify_bioproject = raw_input('\n\
        Do you want to create/update the Assembly table? (Y/N)')


    #---------------------------Bioproject Table------------------------#
    if modify_bioproject.lower() == 'y':
        genomebioproject.BioProjectTable(dbName, ORGANISM, EMAIL)

    #------------------------Biosample User Input----------------------#
    modify_biosample = raw_input('\n\
    Do you want to create/update the BioSample table? (Y/N)')
    while modify_biosample.lower() != 'y' and modify_biosample.lower() != 'n':
        print ("\tInvalid input, must be either 'Y','y','N', or 'n'.")
        modify_biosample = raw_input('\n\
        Do you want to create/update the BioSample table? (Y/N)')


    #---------------------------Biosample Table------------------------#
    if modify_biosample.lower() == 'y':
        genomebiosample.BioSampleTable(dbName, ORGANISM, EMAIL)


    #-----------------------Nucleotide User Input----------------------#
    modify_nucleotide = raw_input('\n\
    Do you want to create/update the Nucleotide table and download genomes? (Y/N)')
    while modify_nucleotide.lower() != 'y' and modify_nucleotide.lower() != 'n':
        print ("\tInvalid input, must be either 'Y','y','N', or 'n'.")
        modify_nucleotide = raw_input('\n\
        Do you want to create/update the nucleotide table and download genomes? (Y/N)')


    #--------------------------Nucleotide Table------------------------#
    if modify_nucleotide.lower() == 'y':
        genomenucleotide.NucleotideTable(dbName, ORGANISM, EMAIL)
