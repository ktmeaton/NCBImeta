
"""
Created on June 12, 2017
NCBI SRA Table Generator

@author: Katherine Eaton
"""

import sqlite3                                                                 # SQL database functionality
import datetime                                                                # Date and time for log files
import os

from Bio import Entrez                                                         # Entrez NCBI API from biopython
from genomeutilities import metadata_count, os_check

def SRATable(dbName, ORGANISM, EMAIL):
    ''' '''
    print("\nCreating/Updating the BioProject table using the following parameters: ")
    print("\t" + dbName)
    print("\t" + ORGANISM)
    print("\t" + EMAIL +"\n\n")

    Entrez.email = EMAIL

    OS_SEP = os_check()                                                        # Retrieve the directory separator by OS

    #-----------------------------------------------------------------------#
    #                                File Setup                             #
    #-----------------------------------------------------------------------#
    if not os.path.exists("log"):                                              # Check if log directory exists
        os.makedirs("log")

    str_sra_log_file = "log" + OS_SEP + ORGANISM.replace(" ", "_") + "_db_sra.log"

    if os.path.exists(str_bioproject_log_file):
        sra_log_file = open(str_sra_log_file, "a")               # Open logfile for appending
    else:
        sra_log_file = open(str_sra_log_file, "w")               # Open logfile for writing


    #-----------------------------------------------------------------------#
    #                                SQL Setup                              #
    #-----------------------------------------------------------------------#

    conn = sqlite3.connect(dbName)                                             # Connect to the DB
    cur = conn.cursor()                                                        # Create cursor for commands
    

    #---------------------------SRA Table----------------------------#
    cur.execute('''
    Create TABLE IF NOT EXISTS BioProject (bioproject_id TEXT,
                                           accession TEXT,
                                           status TEXT,
                                           data_type TEXT,
                                           objectives TEXT,
                                           project_scope TEXT,
                                           description TEXT,
                                           organization TEXT,
                                           registration_date TEXT)''')



 
    # ---------------------Get Bioproject ID-----------------------#
    search_term = ORGANISM + "[Orgn] AND " + "PRJNA269675"  + "[Bioproject]"
    handle = Entrez.esearch(db="sra",
                term=search_term,
                retmax = 1)
    record = Entrez.read(handle)
    ID = record['IdList'][0]
    
    #-------------------------Bioproject Record--------------------#
    ID_handle = Entrez.esummary(db="sra",id=ID)                 # Search for biorproject entry using ID
    ID_record = Entrez.read(ID_handle)                                 # Read in the search results

    key_run = ID_record[0]["Runs"]
    print(key_run)
    '''
    record_dict = ID_record['DocumentSummarySet']['DocumentSummary'][0]                                        # Store metadata as dictionary

    for item in record_dict:
        if isinstance(item, basestring):
            print(item + "\t" + record_dict[item])
    '''
SRATable("abc", "Yersinia pestis", "ktmeaton@gmail.com")
