# -*- coding: utf-8 -*-
"""
Created on Sat Jan 30 17:55:59 2016

NCBI BioSample Table Generator

@author: Katherine Eaton
"""

import sqlite3                                                                 # SQL database functionality
import datetime                                                                # Date and time for log files
import os


from Bio import Entrez                                                         # Entrez NCBI API from biopython
from genomeutilities import sampledata_retrieve,os_check

def BioSampleTable(dbName, ORGANISM, EMAIL, output_dir):
    ''' '''
    print("\nCreating/Updating the BioSample table using the following parameters: ")
    print("Database: " + "\t" + dbName)
    print("Organism: " + "\t" + ORGANISM)
    print("Email: " + "\t" + EMAIL)
    print("Output Directory: " + "\t" + output_dir + "\n\n")

    Entrez.email = EMAIL

    OS_SEP = os_check()                                                        # Retrieve the directory separator by OS

    #-----------------------------------------------------------------------#
    #                                File Setup                             #
    #-----------------------------------------------------------------------#

    log_path = output_dir + OS_SEP + "log"
    
    str_biosample_log_file = "log" + OS_SEP + ORGANISM.replace(" ", "_") + "_db_biosample.log"

    if os.path.exists(str_biosample_log_file):
        biosample_log_file = open(str_biosample_log_file, "a")                # Open logfile for appending
    else:
        biosample_log_file = open(str_biosample_log_file, "w")                 # Open logfile for writing


    #-----------------------------------------------------------------------#
    #                                SQL Setup                              #
    #-----------------------------------------------------------------------#

    conn = sqlite3.connect(dbName)                                             # Connect to the DB
    cur = conn.cursor()                                                        # Create cursor for commands


    #---------------------------BioSample Table----------------------------#
    cur.execute('''
    Create TABLE IF NOT EXISTS BioSample (accession TEXT,
                         organism TEXT,
                         strain TEXT,
                         sample_title TEXT,
                         organization TEXT,
                         collection_date TEXT,
                         geographic_location TEXT,
                         collected_by TEXT,
                         isolate_source TEXT,
                         habitat TEXT,
                         latitude TEXT,
                         longitude TEXT,
                         latitude_and_longitude TEXT,
                         host TEXT,
                         host_disease TEXT,
                         host_status TEXT,
                         biovar TEXT,
                         biotype TEXT,
                         biogroup TEXT,
                         description TEXT,
                         publication_date TEXT,
                         modification_date TEXT)''')


    #-----------------------------------------------------------------------#
    #                         Create Assembly List                          #
    #-----------------------------------------------------------------------#

    cur.execute('''SELECT accession,
                            biosample
                            FROM Assembly''')                                  # Select all assembly accession numbers

    assembly_list = cur.fetchall()                                             # Store assembly accessions as list
    num_assembly = len(assembly_list)                                          # Number of assemblies to process
    num_assembly_processed = 0                                                 # Counter for the number of assemblies processed



    for assembly in assembly_list:
        #-------------------Progress Log and Entry Counter-------------------#
        num_assembly_processed += 1                                            # Increment assembly entry counter

        print("Processing BioSample Accession: " +
                    str(num_assembly_processed) +
                    "/" +
                    str(num_assembly))                                         # Print record progress to screen

        asm_accession = str(assembly[0])                                       # Assembly accession string
        asm_biosample = str(assembly[1])                                       # Assembly biosample (for entrez search term)

        if asm_biosample == "":
            continue

        # ---------------------Check if record exists-----------------------#
        cur.execute('''
        SELECT EXISTS(SELECT accession
                            FROM BioSample
                            WHERE accession=?)''',
                            (asm_biosample,))                                  # Check if biosample record is already in BioSample Table

        record_exists = cur.fetchone()[0]                                      # 0 if not found, 1 if found


        if record_exists:
            continue
        '''
        IMPORTANT:
        The bioproject accession should not exists in the BioProject table
        UNLESS the BioProject record was fully parsed.
        ie. The database does not get updated until the end of each record.
        '''

        # ---------------------Get Biosample ID-----------------------#
        search_term = ORGANISM + "[Orgn] AND " + asm_biosample + "[Biosample]"
        handle = Entrez.esearch(db="biosample",
                    term=search_term,
                    retmax = 1)
        record = Entrez.read(handle)
        ID = record['IdList'][0]

        #-------------------------Bioproject Record--------------------#
        ID_handle = Entrez.esummary(db="biosample",id=ID)                                   # Search for biorproject entry using ID
        ID_record = Entrez.read(ID_handle, validate = False)                                 # Read in the search results
        record_dict = ID_record['DocumentSummarySet']['DocumentSummary'][0]                                        # Store metadata as dictionary

        # -----------------------BioSample attributes-----------------------#
        accession = record_dict['Accession']
        organism = record_dict['Organism']

        sampledata = record_dict['SampleData']


        # Try and get strain from the SampleData
        strain = sampledata_retrieve(sampledata, "strain")
        # If not possible, try and get from the identifiers
        if not strain:
           identifiers = record_dict['Identifiers'].split("; ")
           for element in identifiers:
              if element.startswith("Sample name"):
                 split_element = element.split()
                 strain = split_element[len(split_element) - 1]

        sample_title = record_dict['Title'].replace(";","")
        org = record_dict['Organization']
        col_date = sampledata_retrieve(sampledata, "collection_date")
        geo = sampledata_retrieve(sampledata, "geo")
        collected_by = sampledata_retrieve(sampledata, "collected_by")
        isolate = sampledata_retrieve(sampledata, "isolate_source")
        habitat = sampledata_retrieve(sampledata, "habitat")
        latitude = sampledata_retrieve(sampledata, "latitude")
        longitude = sampledata_retrieve(sampledata, "longitude")
        lat_long = sampledata_retrieve(sampledata, "lat_long")
        host = sampledata_retrieve(sampledata, "host")
        host_disease = sampledata_retrieve(sampledata, "host_disease")
        host_status = sampledata_retrieve(sampledata, "host_status")
        biovar = sampledata_retrieve(sampledata, "biovar")
        biotype = sampledata_retrieve(sampledata, "biotype")
        biogroup = sampledata_retrieve(sampledata, "group")
        description = sampledata_retrieve(sampledata, "description")
        mod_date = record_dict['ModificationDate']
        pub_date = record_dict['PublicationDate']

        # --------------------------Update Database--------------------------#
        cur.execute("SELECT EXISTS(SELECT accession FROM BioSample WHERE accession=?)", (accession,))
        record_exists = cur.fetchone()[0]

        if not record_exists:
          # Write to database
          print ("Writing: " + accession + " to the database.\n")
          cur.execute('''
          INSERT INTO BioSample (accession,
                             organism,
                             strain,
                             sample_title,
                             organization,
                             collection_date,
                             geographic_location,
                             collected_by,
                             isolate_source,
                             habitat,
                             latitude,
                             longitude,
                             latitude_and_longitude,
                             host,
                             host_disease,
                             host_status,
                             biovar,
                             biotype,
                             biogroup,
                             description,
                             publication_date,
                             modification_date)
                                  VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?,
                                  ?, ?, ?, ?, ?, ?, ?, ?, ?, ?,
                                  ?, ?)''',
                             (accession,
                             organism,
                             strain,
                             sample_title,
                             org,
                             col_date,
                             geo,
                             collected_by,
                             isolate,
                             habitat,
                             latitude,
                             longitude,
                             lat_long,
                             host,
                             host_disease,
                             host_status,
                             biovar,
                             biotype,
                             biogroup,
                             description,
                             pub_date,
                             mod_date))

          # Write to logfile
          now = datetime.datetime.now()
          biosample_log_file.write("[" + str(now) + "]" +
                         "\t" + "New accession number added:" +
                         "\t" + accession + "." + "\n")

          conn.commit()



    #----------------------------------------------------------------------#
    #                                    Cleanup                           #
    #----------------------------------------------------------------------#
    conn.commit()                                                              # Make sure all changes are committed to database
    cur.close()                                                                # Close the database
    biosample_log_file.close()                                                 # Close the logfile
