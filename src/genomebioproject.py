# -*- coding: utf-8 -*-
"""
Created on Sat Jan 30 17:20:43 2016

NCBI BioProject Table Generator

@author: Katherine Eaton
"""

import sqlite3                                                                 # SQL database functionality
import datetime                                                                # Date and time for log files
import os


from Bio import Entrez                                                         # Entrez NCBI API from biopython
from genomeutilities import metadata_count, os_check

def BioProjectTable(dbName, ORGANISM, EMAIL):
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

    str_bioproject_log_file = "log" + OS_SEP + ORGANISM.replace(" ", "_") + "_db_bioproject.log"

    if os.path.exists(str_bioproject_log_file):
        bioproject_log_file = open(str_bioproject_log_file, "a")               # Open logfile for appending
    else:
        bioproject_log_file = open(str_bioproject_log_file, "w")               # Open logfile for writing


    #-----------------------------------------------------------------------#
    #                                SQL Setup                              #
    #-----------------------------------------------------------------------#

    conn = sqlite3.connect(dbName)                                             # Connect to the DB
    cur = conn.cursor()                                                        # Create cursor for commands


    #---------------------------BioProjects Table----------------------------#
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




    #-----------------------------------------------------------------------#
    #                         Create Assembly List                          #
    #-----------------------------------------------------------------------#

    cur.execute('''SELECT accession,
                            genbank_bioproject,
                            refseq_bioproject
                            FROM Assembly''')                                  # Select all assembly accession numbers

    assembly_list = cur.fetchall()                                             # Store assembly accessions as list
    num_assembly = len(assembly_list)                                          # Number of assemblies to process
    num_assembly_processed = 0                                                 # Counter for the number of assemblies processed



    for assembly in assembly_list:
        #-------------------Progress Log and Entry Counter-------------------#
        num_assembly_processed += 1                                            # Increment assembly entry counter

        print("Processing BioProject Accession: " +
                    str(num_assembly_processed) +
                    "/" +
                    str(num_assembly))                                         # Print record progress to screen

        asm_accession = str(assembly[0])                                       # Assembly accession string
        asm_bioproj = str(assembly[1])                                         # Assembly bioproject (for entrez search term)

        if asm_bioproj == "":
            continue

        # ---------------------Check if record exists-----------------------#
        cur.execute('''
        SELECT EXISTS(SELECT accession
                            FROM BioProject
                            WHERE accession=?)''',
                            (asm_bioproj,))                                 # Check if bioproject record is already in BioProject Table

        record_exists = cur.fetchone()[0]                                      # 0 if not found, 1 if found


        if not record_exists:
            '''
            IMPORTANT:
            The bioproject accession should not exists in the BioProject table
            UNLESS the BioProject record was fully parsed.
            ie. The database does not get updated until the end of each record.
            '''

            # ---------------------Get Bioproject ID-----------------------#
            search_term = ORGANISM + "[Orgn] AND " + asm_bioproj + "[Bioproject]"
            handle = Entrez.esearch(db="bioproject",
                        term=search_term,
                        retmax = 1)
            record = Entrez.read(handle)
            ID = record['IdList'][0]

            #-------------------------Bioproject Record--------------------#
            ID_handle = Entrez.esummary(db="bioproject",id=ID)                 # Search for biorproject entry using ID
            ID_record = Entrez.read(ID_handle)                                 # Read in the search results
            record_dict = ID_record['DocumentSummarySet']['DocumentSummary'][0]                                        # Store metadata as dictionary


            # -----------------------Bioproject attributes-----------------#
            accession = record_dict['Project_Acc']
            status = record_dict['Sequencing_Status']
            data_type = record_dict['Project_Data_Type']


            # Objectives
            obj_list = record_dict['Project_Objectives_List']
            obj_str = ""
            obj_counter = 1
            for item in obj_list:
               if obj_counter < len(obj_list):
                  obj_str += item['Project_ObjectivesType'] + ", "
               else:
                  obj_str += item['Project_ObjectivesType']
               obj_counter += 1


            proj_scope = record_dict['Project_Target_Scope']
            description = record_dict['Project_Description']
            organization = record_dict['Submitter_Organization']
            reg_date = record_dict['Registration_Date']


            # --------------------------Update Database--------------------------#
            print ("Writing: " + accession + " to the database.\n")
            cur.execute('''
            INSERT INTO BioProject (bioproject_id,
                                  accession,
                                  status,
                                  data_type,
                                  objectives,
                                  project_scope,
                                  description,
                                  organization,
                                  registration_date)
                                  VALUES
                                  (?, ?, ?, ?, ?, ?, ?, ?, ?)''',
                             (ID,
                             accession,
                              status,
                              data_type,
                              obj_str,
                              proj_scope,
                              description,
                              organization,
                              reg_date))

            # Write to logfile
            now = datetime.datetime.now()
            bioproject_log_file.write("[" + str(now) + "]" +
                         "\t" + "New accession number added:" +
                         "\t" + accession + "." + "\n")
            conn.commit()



    #----------------------------------------------------------------------#
    #                                    Cleanup                           #
    #----------------------------------------------------------------------#
    conn.commit()                                                              # Make sure all changes are committed to database
    cur.close()                                                                # Close the database
    bioproject_log_file.close()                                                # Close the logfile
