# -*- coding: utf-8 -*-
"""
NCBI Assembly Table Generator

"""

import sqlite3                                                                 # SQL database functionality
import datetime                                                                # Date and time for log files
import os


from Bio import Entrez                                                         # Entrez NCBI API from biopython
from genomeutilities import metadata_count, os_check


def AssemblyTable(dbName, ORGANISM, EMAIL):
    ''' '''
    print("\nCreating/Updating the Assembly table using the following parameters: ")
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

    str_assembly_log_file = "log" + OS_SEP + ORGANISM.replace(" ", "_") + "_db_assembly.log"

    if os.path.exists(str_assembly_log_file):
        assembly_log_file = open(str_assembly_log_file, "a")                   # Open logfile for appending
    else:
        assembly_log_file = open(str_assembly_log_file, "w")                   # Open logfile for writing
    #-----------------------------------------------------------------------#
    #                                SQL Setup                              #
    #-----------------------------------------------------------------------#

    conn = sqlite3.connect(dbName)                                             # Connect to the DB
    cur = conn.cursor()                                                        # Create cursor for commands

    #-----------------------------------------------------------------------#
    #                             Assembly Table                            #
    #-----------------------------------------------------------------------#
    cur.execute('''
    Create TABLE IF NOT EXISTS Assembly (id INTEGER NOT NULL PRIMARY KEY AUTOINCREMENT UNIQUE,
                                         assembly_id TEXT,
                                         accession TEXT UNIQUE,
                                         organism TEXT,
                                         infraspecies TEXT,
                                         assembly_status TEXT,
                                         organization TEXT,
                                         taxid TEXT,
                                         biosample INTEGER,
                                         coverage FLOAT,
                                         replicon_count int,
                                         chromosome_count int,
                                         non_chromosome_replicon_count int,
                                         contig_count int,
                                         contig_l50_count int,
                                         contig_n50_count int,
                                         scaffold int,
                                         scaffold_l50_count int,
                                         scaffold_n50_count int,
                                         total_length int,
                                         ungapped_length int,
                                         refseq_bioproject INTEGER,
                                         genbank_bioproject INTEGER,
                                         wgs_project TEXT,
                                         submission_date TEXT,
                                         sequence_release_date TEXT,
                                         assembly_release_date TEXT)''')



    #-----------------------------------------------------------------------#
    #                           Entrez Assemblies                           #
    #-----------------------------------------------------------------------#
    search_term = ORGANISM + "[Orgn]"
    handle = Entrez.esearch(db="assembly",
                            term=search_term,
                            retmax = 1000)

    record = Entrez.read(handle)                                                   # Read the record
    num_records = int(record['Count'])                                             # Count the total number of entries
    num_processed = 0                                                              # Counter for num entries processed


    #-----------------------------------------------------------------------#
    #                          Iterate Through ID List                      #
    #-----------------------------------------------------------------------#

    for ID in record['IdList']:

       #-------------------Progress Log and Entry Counter-------------------#
       num_processed += 1                                                          # Increment entry counter
       print("Processing record: " +
           str(num_processed) + \
           "/" + str(num_records))                                                 # Print record progress to screen


       #------------Check if Assembly Already Exists in Database------------#
       cur.execute('''
       SELECT EXISTS(SELECT assembly_id FROM Assembly WHERE assembly_id=?)''',\
       (ID,))                                                                      # Check if assembly exists in DB
       record_exists = cur.fetchone()[0]                                           # 0 if not found, 1 if found


       #---------------If Assembly Isn't in Database, Add it------------#
       if not record_exists:

           ID_handle = Entrez.esummary(db="assembly",id=ID)                        # Retrieve Assembly record using ID
           ID_record = Entrez.read(ID_handle, validate=False)                      # Read Assembly record
           record_dict = ID_record['DocumentSummarySet']['DocumentSummary'][0]     # Store Assembly information as dictionary


           tmp_file = open("assembly_record_dict.txt",'w')
           entrez_string_element = type(record_dict['Taxid'])                      # Get type of Bio String object (used in type checking)


           # -----------------------Asssembly attributes-----------------------#
           accession = record_dict['AssemblyAccession']                            # Assembly accession number
           organism = record_dict['Organism']                                      # Organism name (sometimes includes straing)
           asm_status = record_dict['AssemblyStatus']                              # Status of assembly (complete, chromosome, scaffold, contig)
           organization = record_dict['SubmitterOrganization']                     # Organization that submitted the record
           taxid = record_dict['Taxid']                                            # NCBI Taxonomic ID number (specific to strain!)
           biosample = record_dict['BioSampleAccn']                                # Biosample accession number
           coverage = float(record_dict['Coverage'])                               # Average coverage of assembly
           sub_date = record_dict['SubmissionDate']                                # Submission date
           seq_rel_date = record_dict['SeqReleaseDate']                            # Sequence release date
           asm_rel_date= record_dict['SubmissionDate']                             # Assembly release date
           wgs_proj = record_dict['WGS']                                           # Whole genome sequence project accession
           meta_string = record_dict['Meta']                                       # Metadata as entered by submitter



           # -----------------------Metadata attributes-----------------------#
           chromosome_count = metadata_count(meta_string, "chromosome")            # Number of chromosomes
           replicon_count = metadata_count(meta_string, "replicon")                # Total number of replicons
           non_chr_replicon_count = metadata_count(meta_string, "non_chr_rep")     # Total number of non-chromosome replicons
           contig_count = metadata_count(meta_string, "contig")                    # Total number of contigs
           contig_l50_count = metadata_count(meta_string, "contig_l50")            # Contig L50
           contig_n50_count = metadata_count(meta_string, "contig_n50")            # Contig N50
           scaffold_count = metadata_count(meta_string, "scaffold")                # Total number of scaffols
           scaffold_l50_count = metadata_count(meta_string, "scaffold_l50")        # Scaffold L50
           scaffold_n50_count = metadata_count(meta_string, "scaffold_n50")        # Scaffold N50
           total_length = metadata_count(meta_string, "length")                    # Total assembly length
           ungapped_length= metadata_count(meta_string, "ungap_len")               # Ungapped length



           # --------------------Inconsistent attributes---------------------#
           if record_dict['Biosource']['InfraspeciesList']:                        # Useful when 'Organism' does not include strain info
               infraspecies = record_dict['Biosource']['InfraspeciesList'][0]['Sub_value']
           else: infraspecies = ""                                                 # Sometimes the Infraspecies key is blank

           if record_dict['GB_BioProjects']:                                       # Genbank specific bioproject accession
               bioproj_gb_accession = record_dict['GB_BioProjects'][0]['BioprojectAccn']
           else: bioproj_gb_accession = ""                                         # Sometimes no Genabnk bioproject is associated

           if record_dict['RS_BioProjects']:                                       # Refseq specific bioproject accession
               bioproj_rs_accession = record_dict['RS_BioProjects'][0]['BioprojectAccn']
           else: bioproj_rs_accession = ""                                         # Sometimes no Refseq bioproject is associated





           # ------------------------Update Database-------------------------#
           print ("Writing: " + accession + " to the database.\n")                 # Print progress log to screen
           cur.execute('''
           INSERT INTO Assembly (assembly_id,
                                 accession,
                                organism,
                                infraspecies,
                                assembly_status,
                                organization,
                                taxid,
                                biosample,
                                coverage,
                                replicon_count,
                                chromosome_count,
                                non_chromosome_replicon_count,
                                contig_count,
                                contig_l50_count,
                                contig_n50_count,
                                scaffold,
                                scaffold_l50_count,
                                scaffold_n50_count,
                                total_length,
                                ungapped_length,
                                refseq_bioproject,
                                genbank_bioproject,
                                wgs_project,
                                submission_date,
                                sequence_release_date,
                                assembly_release_date)
                                VALUES (
                                ?, ?, ?, ?, ?, ?, ?, ?, ?, ?,
                                ?, ?, ?, ?, ?, ?, ?, ?, ?, ?,
                                ?, ?, ?, ?, ?, ?)''',
                             (ID,
                             accession,
                             organism,
                             infraspecies,
                             asm_status,
                             organization,
                             taxid,
                             biosample,
                             coverage,
                             replicon_count,
                             chromosome_count,
                             non_chr_replicon_count,
                             contig_count,
                             contig_l50_count,
                             contig_n50_count,
                             scaffold_count,
                             scaffold_l50_count,
                             scaffold_n50_count,
                             total_length,
                             ungapped_length,
                             bioproj_rs_accession,
                             bioproj_gb_accession,
                             wgs_proj,
                             sub_date,
                             seq_rel_date,
                             asm_rel_date))



           # ------------------------WWrite to Logfile-------------------------#
           now = datetime.datetime.now()                                           # Get timestamp
           assembly_log_file.write("[" + str(now) + "]" +
                         "\t" + "New accession number added:" +
                         "\t" + accession + "." + "\n")                            # Write to logfile
           conn.commit()                                                           # Commit assembly record changes to database





    #-----------------------------------------------------------------------#
    #                                    Cleanup                            #
    #-----------------------------------------------------------------------#
    conn.commit()                                                                  # Make sure all changes are committed to database
    cur.close()                                                                    # Close the database
    assembly_log_file.close()                                                      # Close the logfile
