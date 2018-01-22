# -*- coding: utf-8 -*-
# -*- coding: utf-8 -*-
"""
Yersinia pestis nucleotide database generator.

"""

import sqlite3                                                                 # SQL database functionality
import os                                                                      # Create directories
import time
from genomeutilities import *
import datetime                                                                # Date and time for log files
from Bio import Entrez                                                         # Entrez NCBI API from biopython
from Bio import SeqIO
from Bio.SeqIO import FastaIO

from genomeutilities import os_check


def NucleotideTable(dbName, ORGANISM, EMAIL, output_dir):
    ''' '''
    print("\nCreating/Updating the Nucleotide table using the following parameters: ")
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
    str_nucleotide_log_file = log_path + OS_SEP + ORGANISM.replace(" ", "_") + "_db_nucleotide.log"
    str_genbank_log_file = log_path + OS_SEP + ORGANISM.replace(" ", "_") + "_db_genbank.log"

    if os.path.exists(str_nucleotide_log_file):
        nucleotide_log_file = open(str_nucleotide_log_file, "a")               # Open logfile for appending
    else:
        nucleotide_log_file = open(str_nucleotide_log_file, "w")               # Open logfile for writing

    if os.path.exists(str_genbank_log_file):
        genbank_log_file = open(str_genbank_log_file,"a")
    else:
        genbank_log_file = open(str_genbank_log_file,"w")


    data_path = output_dir + OS_SEP + "data" + OS_SEP                              # Path to the data folder
    genbank_path = output_dir + OS_SEP + "genbank" + OS_SEP


    #-----------------------------------------------------------------------#
    #                                SQL Setup                              #
    #-----------------------------------------------------------------------#

    conn = sqlite3.connect(dbName)                                             # Connect to the DB
    cur = conn.cursor()                                                        # Create cursor for commands

    #cur.execute(''' Drop TABLE IF EXISTS Nucleotide''')
    #cur.execute(''' Drop TABLE IF EXISTS Fasta''')
    #cur.execute(''' Drop TABLE IF EXISTS Genbank''')   
    #return 0

    #-----------------------------------------------------------------------#
    #                            Nucleotide Table                           #
    #-----------------------------------------------------------------------#
    cur.execute('''
    Create TABLE IF NOT EXISTS Nucleotide (asm_accession TEXT,
                                           num_records TEXT)''')                   # Table that holds assembly accessions and number of nucleotide records

    cur.execute('''
    Create TABLE IF NOT EXISTS Fasta (accession TEXT,
                                      GI TEXT,
                                      title TEXT,
                                      taxid TEXT,
                                      length int,
                                      create_date TEXT,
                                      update_date TEXT
                                      )''')                                        # Table that holds nucleotide accessions  and contig information

    cur.execute('''
    Create TABLE IF NOT EXISTS Genbank (accession TEXT,
					GI TEXT
					)''')				# Table that holds genbank accessions and GI 



    #-----------------------------------------------------------------------#
    #                         Create Assembly List                          #
    #-----------------------------------------------------------------------#

    cur.execute('''SELECT accession,
                            organism,
                            genbank_bioproject,
                            biosample,
                            contig_count,
                            refseq_bioproject
                            FROM Assembly''')                                      # Select all assembly accession numbers

    assembly_list = cur.fetchall()                                                 # Store assembly accessions as list
    assembly_list = sorted(assembly_list, key=lambda tup: tup[4])

    num_assembly = len(assembly_list)                                              # Number of assemblies to process
    num_assembly_processed = 0                                                     # Counter for the number of assemblies processed
    


    for assembly in assembly_list:                                       # I choose to start from the bottom of the list

        #-------------------Progress Log and Entry Counter-------------------#
        num_assembly_processed += 1                                                # Increment assembly entry counter

        print("Processing Assembly Accession: " +
                    str(num_assembly_processed) +
                    "/" +
                    str(num_assembly))                                             # Print record progress to screen

        asm_accession = str(assembly[0])                                           # Assembly accession string
        organism = str(assembly[1]).replace(" ", "_")                              # Assembly organism name (for file output naming)
        asm_bioproj = str(assembly[2])                                             # Assembly bioproject (for entrez search term)
        asm_biosample = str(assembly[3])                                           # Assembly biosample (for entrez search term)
        asm_contig_count = assembly[4]                                             # Assembly contig count (for ret max)
        asm_refseq_bioproj = str(assembly[5])

        # check for newer version
        bool_found_newer = 0
        if asm_accession.endswith(".1"):
             asm_accession_newer = asm_accession.split(".")[0] + ".2"
             for asm_check in assembly_list:
                 if asm_accession_newer == asm_check[0]:
                     print("Skipping. Newer accession: " + asm_accession_newer + " found.")
                     bool_found_newer = 1             
        if bool_found_newer == 1 : continue    

        # Change bioproject listing to refseq if available
        # This is very important to get correct contig number for
        # assembly updates
        #if asm_refseq_bioproj:
        #    asm_bioproj = asm_refseq_bioproj
             
        #print(asm_bioproj)
        # ---------------------Check if record exists-----------------------#
        cur.execute('''
        SELECT EXISTS(SELECT asm_accession
                            FROM Nucleotide
                            WHERE asm_accession=?)''',
                            (asm_accession,))                                      # Check if assembly record is alread in Nucleotide DB

        record_exists = cur.fetchone()[0]                                          # 0 if not found, 1 if found


        if record_exists:
            continue
            '''
            IMPORTANT:
            The assembly accession should not exists in the Nucleotide table UNLESS
            all associated fasta contigs have already been download.
            ie. The database does not get updated until the end of each record.
            '''
            
        #----------------------------------------------------------------#
        #                           Entrez Nucleotide                    #
        #----------------------------------------------------------------#

        assembly_directory = data_path + OS_SEP + organism + "_" + asm_accession         # Check if assembly directory already exists
        if not os.path.exists(assembly_directory):
            os.makedirs(assembly_directory)                                    # If not, create assembly directory

        assembly_gb_directory = genbank_path + OS_SEP + organism + "_" + asm_accession 
        if not os.path.exists(assembly_gb_directory):
            os.makedirs(assembly_gb_directory)

        search_term = (ORGANISM + "[Orgn] AND "  +
                        asm_bioproj + "[Bioproject] AND " +
                        asm_biosample + "[Biosample] " +
                        "NOT project[Title]")                                  # Search term using organism, bioproject, biosample
                                                                               # I exclude "project" in title so as to not receive the WGS record entry
        print(search_term)
        assembly_nt_handle = Entrez.esearch(db="nucleotide",
                                term=search_term,
                                retmax = asm_contig_count)                     # Search for the assembly nt data, number returned is the number of contigs/scaffolds/replicons

        assembly_nt_record = Entrez.read(assembly_nt_handle)                   # Read in the search results


        #-------------------Progress Log and Entry Counter---------------#
        num_assembly_nt_records = len(assembly_nt_record['IdList'])             # Count the number of contigs/scaffols/replicons
        num_assembly_nt_processed = 0                                          # Counter for the number of nucleotide entries processed
        print("Number of records:\t", num_assembly_nt_records)                 # Print  number of nucleotide entries to screen



        #----------------------------------------------------------------#
        #                           Entrez Sequence                      #
        #----------------------------------------------------------------#
        for ID in assembly_nt_record["IdList"]:
            #-------------------Progress Log and Entry Counter-----------#
            num_assembly_nt_processed += 1                                     # Increment nucleotide entry counter
            print("Processing record: " +
                    str(num_assembly_nt_processed) +
                    "/" +
                    str(num_assembly_nt_records))                              # Print nucleotide progress to screen

            #------------------------Nucleotide Record-------------------#
            ID_handle = Entrez.esummary(db="nucleotide",id=ID)                 # Search for nucleotide entry using contig ID
            ID_record = Entrez.read(ID_handle)                                 # Read in the search results
            record_dict = ID_record[0]                                         # Store metadata as dictionary
            ID_accession = record_dict['Caption']                              # Retrieve the nucleotide accession number


            # ------------------Check if record exists-------------------#
            cur.execute('''SELECT EXISTS(
                            SELECT accession
                            FROM Fasta
                            WHERE accession=?)''',
                            (ID_accession,))                                   # Check if contig entry is already in Fasta Table

            record_exists = cur.fetchone()[0]                                  # 0 if not found, 1 if found

            if record_exists:
                continue
                '''
                IMPORTANT:
                The contig accession should not exists in the Fasta table UNLESS
                that fasta file has been downloaded
                ie. The database does not get updated until the end of each contig.
                '''

            # --------------------Get Attributes---------------------#
            ID_GI = record_dict['Gi']
            ID_title = record_dict['Title']
            ID_taxid = record_dict['TaxId']
            ID_length = int(record_dict['Length'])
            ID_create_date = record_dict['CreateDate']
            ID_update_date = record_dict['UpdateDate']



            # -------------------Create fasta filename---------------#
            fasta_file_name = record_dict['Title'].replace(" ", "_")       # Convert all spaces to underscore for file name
            fasta_file_name = fasta_file_name.strip(
                        ',_whole_genome_shotgun_sequence')                 # Strip off the end, it's uninformative and makes the filename too lon
            fasta_file_name = fasta_file_name.replace("/","_")
            fasta_file_name = fasta_file_name + ".fasta"                   # Add the .fasta extension
            fasta_file_path = assembly_directory + OS_SEP + fasta_file_name  # Create the fasta file path

            # ------------------Create genbank filename--------------#
            genbank_file_name = fasta_file_name.replace(".fasta",".gb")
            genbank_file_path = assembly_gb_directory + OS_SEP + genbank_file_name
            

   

            # ----------------Check if fasta file exists-------------------#
            if not os.path.exists(fasta_file_path):
                #!#!#!#!#!
                time.sleep(30)
                #!#!#!#!#!


                print("\tDownloading:\t" + fasta_file_name)                          # Print progress to screen

                fasta_file = open(fasta_file_path, "w")                    # Open the fasta file
                fasta_handle = Entrez.efetch(db="nucleotide",
                                             id=ID, rettype="fasta",
                                             retmode="text")               # Download the fasta file

                seq_record = SeqIO.parse(fasta_handle, "fasta")            # Parse out the fasta record
                list_seq_records = list(seq_record)[0]
                seq_id = list_seq_records.id
                seq_seq = list_seq_records.seq
                header = ">gi|" + str(ID_GI) + "||" + seq_id + "|" + fasta_file_name.replace(".fasta","")
                fasta_file.write(header + "\n")
                fasta_file.write(str(seq_seq))
                fasta_file.close()                                         # Close the fasta file
                

            #----------------Check if genbank file exists------------#
            if not os.path.exists(genbank_file_path):
                
                print("\tDownloading:\t" + genbank_file_name)
                genbank_file = open(genbank_file_path,"w")
                genbank_handle = Entrez.efetch(db="nucleotide",
						id=ID,
						rettype="gbwithparts",
						retmode="text")
                genbank_file.write(genbank_handle.read())
                genbank_file.close()


            #--------------------------------------------------------#
            #                   Update FASTA Table                   #
            #--------------------------------------------------------#
            cur.execute('''
            INSERT INTO Fasta (accession,
                              GI,
                              title,
                              taxid,
                              length,
                              create_date,
                              update_date)
                              VALUES
                              (?, ?, ?, ?, ?, ?, ?)''',
                                (ID_accession,
                                 ID_GI,
                                 ID_title,
                                 ID_taxid,
                                 ID_length,
                                 ID_create_date,
                                 ID_update_date))




            #--------------------------------------------------------#
            #                   Update Genbank Table                 #
            #--------------------------------------------------------#
            cur.execute('''
            INSERT INTO Genbank (accession,
                              GI)
                              VALUES
                              (?, ?)''',
                                (ID_accession,
                                 ID_GI))

            '''
            Remove the following line if you are confident the
            downloading
            will proceed 100% without failing, and want to speed up
            the running time.
            '''


            conn.commit()                                                  # Commit changes to the Fasta Table


        #----------------------------------------------------------------#
        #                    Update Nucleotide Table                     #
        #----------------------------------------------------------------#

        print ("Writing: " + asm_accession + " to the database.\n")            # Print progress log

        cur.execute('''
        INSERT INTO Nucleotide (asm_accession,
                               num_records)
                               VALUES
                               (?,?)''',
                               (asm_accession,
                                num_assembly_nt_records))


        # ----------------------Write to Log File------------------------#
        now = datetime.datetime.now()
        nucleotide_log_file.write("[" + str(now) + "]" +
                     "\t" + "New accession number added:" +
                     "\t" + asm_accession + "." + "\n")
        genbank_log_file.write("[" + str(now) + "]" + 
		     "\t" + "New accession number added:" + 
		     "\t" + asm_accession + "." + "\n")
        conn.commit()                                                          # Commit changes to the Nucleotide Table

    #-----------------------------------------------------------------------#
    #                               Cleanup                                 #
    #-----------------------------------------------------------------------#
    conn.commit()
    cur.close()
    nucleotide_log_file.close()
    genbank_log_file.close()

#NucleotideTable("database/Yersinia_pestis_db.sqlite", "Yersinia pestis", "ktmeaton@gmail.com", "./")
