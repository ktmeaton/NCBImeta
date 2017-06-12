# -*- coding: utf-8 -*-
"""
Created on Sat Jan 30 14:09:59 2016
Edited: June 6, 2017 - post-GI number removal check

@author: Katherine Eaton

"Genome Collector Module"
"""

import argparse                                                         # Command-line argument parsing
import sqlite3
import os

from genomeerrors import *
import genomeassembly
import genomeutilities
import genomebioproject
import genomebiosample
import genomenucleotide
import genomesra


#-----------------------------------------------------------------------#
#                            Argument Parsing                           #
#-----------------------------------------------------------------------#

OS_SEP = genomeutilities.os_check()                                             # Retrieve the directory separator by OS

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
                    dest = 'db_name',
                    required=True)

mandatory.add_argument('--outputdir',
                    help='Output Directory',
                    type = str,
                    action = 'store',
                    dest = 'output_dir',
                    required=True)


bonus.add_argument('--assembly-status',
                   help='Filter by Assembly Status: Complete, Chromosome, Scaffold, or Contig',
                   type = str,
                   action = 'store',
                   dest = 'assembly_status',
                   required = False)

args = vars(parser.parse_args())


mode = args['mode']
db_name = args['db_name'] + ".sqlite"
output_dir = args['output_dir']
assembly_status = args['assembly_status']


#-------------------------------------------------------------------------------#
#                           Argument Checking                                   #
#-------------------------------------------------------------------------------#

genomeutilities.check_accessory_dir(output_dir)                                 # Create accessory directories

db_path = output_dir + OS_SEP + "database" + OS_SEP + db_name

#---------------------------Delete Database-------------------------------------#
if mode.lower() == 'delete':                                                    # Delete mode check
    if os.path.exists(db_path):                                                 # Make sure database exists in output directory
        os.remove(db_path)                                                      # Delete database
        print('\nDeleting database: ' + db_path)
    else:
        raise ErrorDBNotExists(db_path)                                        # If database doesn't exists, raise error


#---------------------------Create Database-------------------------------------#
elif mode.lower() == 'create':                                                  # Create mode check
    if not os.path.exists(db_path):                                             # Make sure database exists in output directory
        conn = sqlite3.connect(db_path)                                         # Create database
        conn.commit()
        print('\nCreating database: ' + db_path)
    else:
        raise ErrorDBExists(db_path)                                             # If database already exists, raise error


#---------------------------Update Database-----------------------------#
elif mode.lower() == 'update':
    if os.path.exists(db_path):
        conn = sqlite3.connect(db_path)
        print('\nOpening database: ' + db_path)
    else:
        raise ErrorDBNotExists(db_path)

#----------------------------Invalid Mode------------------------------#
else:
    raise ErrorInvalidMode(mode)



#------------------------Create Mode Processing------------------------#

if mode.lower() == 'create' or mode.lower() == 'update':

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
        genomeassembly.AssemblyTable(db_path, ORGANISM, EMAIL, output_dir)





    #------------------------Bioproject User Input----------------------#
    modify_bioproject = raw_input('\n\
    Do you want to create/update the BioProject table? (Y/N)')
    while modify_bioproject.lower() != 'y' and modify_bioproject.lower() != 'n':
        print ("\tInvalid input, must be either 'Y','y','N', or 'n'.")
        modify_bioproject = raw_input('\n\
        Do you want to create/update the Bioproject table? (Y/N)')


    #---------------------------Bioproject Table------------------------#
    if modify_bioproject.lower() == 'y':
        genomebioproject.BioProjectTable(db_path, ORGANISM, EMAIL, output_dir)






    #----------------------------SRA User Input-------------------------#
    modify_sra = raw_input('\n\
    Do you want to create/update the SRA table? (Y/N)')
    while modify_sra.lower() != 'y' and modify_sra.lower() != 'n':
        print ("\tInvalid input, must be either 'Y','y','N', or 'n'.")
        modify_sra = raw_input('\n\
        Do you want to create/update the SRA table? (Y/N)')


    #-----------------------------SRA Table-----------------------------#
    if modify_sra.lower() == 'y':
        genomesra.SRATable(db_path, ORGANISM, EMAIL, output_dir)





    #------------------------Biosample User Input----------------------#
    modify_biosample = raw_input('\n\
    Do you want to create/update the BioSample table? (Y/N)')
    while modify_biosample.lower() != 'y' and modify_biosample.lower() != 'n':
        print ("\tInvalid input, must be either 'Y','y','N', or 'n'.")
        modify_biosample = raw_input('\n\
        Do you want to create/update the BioSample table? (Y/N)')


    #---------------------------Biosample Table------------------------#
    if modify_biosample.lower() == 'y':
        genomebiosample.BioSampleTable(db_path, ORGANISM, EMAIL, output_dir)






    #-----------------------Nucleotide User Input----------------------#
    modify_nucleotide = raw_input('\n\
    Do you want to create/update the Nucleotide table and download genomes? (Y/N)')
    while modify_nucleotide.lower() != 'y' and modify_nucleotide.lower() != 'n':
        print ("\tInvalid input, must be either 'Y','y','N', or 'n'.")
        modify_nucleotide = raw_input('\n\
        Do you want to create/update the Nucleotide table and download genomes? (Y/N)')


    #--------------------------Nucleotide Table------------------------#
    if modify_nucleotide.lower() == 'y':
        genomenucleotide.NucleotideTable(db_path, ORGANISM, EMAIL, output_dir)

print("Genome Collector module has finished.")
