
"""
Created on June 12, 2017
NCBI SRA Table Generator

@author: Katherine Eaton
"""

import sqlite3                                                              # SQL database functionality
import datetime                                                             # Date and time for log files
import os
import xml.etree.ElementTree                                                # XML Parsing
from bs4 import BeautifulSoup

from Bio import Entrez                                                      # Entrez NCBI API from biopython
from genomeutilities import metadata_count, os_check, parseSRARunInfo, sra_xml_parse


def SRATable(dbName, ORGANISM, EMAIL, output_dir):
    ''' '''
    print("\nCreating/Updating the SRA table using the following parameters: ")
    print("Database: " + "\t" + dbName)
    print("Organism: " + "\t" + ORGANISM)
    print("Email: " + "\t" + EMAIL)
    print("Output Directory: " + "\t" + output_dir + "\n\n")

    Entrez.email = EMAIL

    OS_SEP = os_check()                                                     # Retrieve the directory separator by OS


    #-----------------------------------------------------------------------#
    #                                File Setup                             #
    #-----------------------------------------------------------------------#

    log_path = output_dir + OS_SEP + "log"

    str_sra_log_file = log_path + OS_SEP + ORGANISM.replace(" ", "_") + "_db_sra.log"

    if os.path.exists(str_sra_log_file):
        sra_log_file = open(str_sra_log_file, "a")                          # Open logfile for appending
    else:
        sra_log_file = open(str_sra_log_file, "w")                          # Open logfile for writing

    #-----------------------------------------------------------------------#
    #                                SQL Setup                              #
    #-----------------------------------------------------------------------#

    conn = sqlite3.connect(dbName)                                   # Connect to the DB
    cur = conn.cursor()                                                 # Create cursor for commands


    #---------------------------SRA Table-----------------------------------#
    cur.execute('''
    Create TABLE IF NOT EXISTS SRA (bioproject_id TEXT,
                                           bioproject_accession TEXT,
                                           run_accession TEXT,
                                           create_date TEXT,
                                           update_date TEXT,
                                           ext_links TEXT,
                                           sra_ID TEXT,
                                           total_bases int,
                                           total_spots int,
                                           scientific_name TEXT,
                                           library_layout TEXT,
                                           library_source TEXT,
                                           library_selection TEXT,
                                           name TEXT,
                                           title TEXT,
                                           sra_bioproject TEXT,
                                           sra_biosample TEXT,
                                           total_runs TEXT,
                                           experiment_acc TEXT,
                                           organism_taxid TEXT,
                                           instrument TEXT)''')

    cur.execute('''
    SELECT accession FROM BioProject''')
    bioproject_list = cur.fetchall()

    num_bioproject = len(bioproject_list)                                   # Number of bioprojects to process
    num_bioproject_processed = 0                                            # Counter for the number of bioprojects processed

    #-----------------------------------------------------------------------#
    #                                Processing                             #
    #-----------------------------------------------------------------------#

    for element in bioproject_list:
        bioproject_accession = element[0]

        #-------------------Progress Log and Entry Counter------------------#
        num_bioproject_processed += 1                                       # Increment bioproject entry counter

        print("Processing BioProject Accession: " +
                str(num_bioproject_processed) +
                "/" +
                str(num_bioproject))                                        # Print record progress to screen
        
        # ---------------------Check if record exists-----------------------#
        cur.execute('''
        SELECT EXISTS(SELECT bioproject_accession
                            FROM SRA
                            WHERE bioproject_accession=?)''',
                            (bioproject_accession,))                        # Check if sra record is already in SRA Table

        record_exists = cur.fetchone()[0]                               # 0 if not found, 1 if found

        if record_exists:
            continue                                                        # If bioproject has been fully processed, skip
            '''
            IMPORTANT:
            The SRA accession should not exists in the SRA table
            UNLESS the Bioproject record was fully parsed.
            ie. The database does not get updated until the end of each bioproject collection
            '''        

        # ---------------------Get SRA ID-----------------------------------#

        search_term = ORGANISM + "[Orgn] AND " + bioproject_accession  + "[Bioproject]"
        handle = Entrez.esearch(db="sra",
                    term=search_term,
                    retmax = 1000)
        record = Entrez.read(handle)
        ID_list = record['IdList']

        num_ID = len(ID_list)
        num_ID_processed = 0 
        
        #-------------------------SRA Record--------------------------------#
        for ID in ID_list:

            num_ID_processed += 1                                           # Increment assembly entry counter

            print("Processing SRA ID: " +
                str(num_ID_processed) +
                "/" +
                str(num_ID))                                                # Print record progress to screen
    
            ID_handle = Entrez.esummary(db="sra",id=ID)                     # Search for bioproject entry using ID
            ID_record = Entrez.read(ID_handle)                              # Read in the search results
            sra_metadata = ID_record[0]

            # Dates
            create_date = sra_metadata["CreateDate"]            
            update_date = sra_metadata["UpdateDate"]

            # Misc
            ext_links = sra_metadata["ExtLinks"]             
            exp_xml = sra_metadata["ExpXml"]
            item = sra_metadata["Item"]          
            sra_ID = sra_metadata["Id"]
            
            # Run Info
            run_info = sra_metadata["Runs"]        
            run_info_dict = parseSRARunInfo(run_info)

            run_accession = run_info_dict["Run_acc"]
            run_total_bases = run_info_dict["total_bases"]
            run_total_spots = run_info_dict["total_spots"]

            # Meta
            scientific_name = sra_xml_parse(exp_xml, "ScientificName")
            library_layout = sra_xml_parse(exp_xml, "LIBRARY_LAYOUT")
            library_source = sra_xml_parse(exp_xml, "LIBRARY_SOURCE")
            library_selection = sra_xml_parse(exp_xml, "LIBRARY_SELECTION")
            name = sra_xml_parse(exp_xml, " name")
            title = sra_xml_parse(exp_xml, "Title")
            sra_bioproject = sra_xml_parse(exp_xml, "Bioproject")
            sra_biosample = sra_xml_parse(exp_xml, "Biosample")
            total_runs = sra_xml_parse(exp_xml, "total_runs")
            experiment_acc = sra_xml_parse(exp_xml, "Experiment_acc")
            organism_taxid = sra_xml_parse(exp_xml, "Organism_taxid")
            instrument = sra_xml_parse(exp_xml, "instrument_model")

            # --------------------------Update Database--------------------------#
            print ("Writing: " + run_accession + " to the database.\n")
            cur.execute('''
            INSERT INTO SRA (bioproject_id,
                                       bioproject_accession,
                                       run_accession,
                                       create_date,
                                       update_date,
                                       ext_links,
                                       sra_ID,
                                       total_bases,
                                       total_spots,
                                       scientific_name,
                                       library_layout,
                                       library_source,
                                       library_selection,
                                       name,
                                       title,
                                       sra_bioproject,
                                       sra_biosample,
                                       total_runs,
                                       experiment_acc,
                                       organism_taxid,
                                       instrument)
                                  VALUES
                                  (?, ?, ?, ?, ?,
                                   ?, ?, ?, ?, ?,
                                   ?, ?, ?, ?, ?,
                                   ?, ?, ?, ?, ?,
                                   ?)''',
                             (ID,
                                       bioproject_accession,
                                       run_accession,
                                       create_date,
                                       update_date,
                                       ext_links,
                                       sra_ID,
                                       run_total_bases,
                                       run_total_spots,
                                       scientific_name,
                                       library_layout,
                                       library_source,
                                       library_selection,
                                       name,
                                       title,
                                       sra_bioproject,
                                       sra_biosample,
                                       total_runs,
                                       experiment_acc,
                                       organism_taxid,
                                       instrument))
            

        # Write to logfile
        now = datetime.datetime.now()
        sra_log_file.write("[" + str(now) + "]" +
                     "\t" + "New bioproject accession files added:" +
                     "\t" + bioproject_accession + "." + "\n")
        conn.commit()

                
    #----------------------------------------------------------------------#
    #                                    Cleanup                           #
    #----------------------------------------------------------------------#
    conn.commit()                                                      # Make sure all changes are committed to database
    cur.close()                                                        # Close the database
    sra_log_file.close()                                                   # Close the logfile


