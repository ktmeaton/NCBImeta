"""
Created on Thursday June 08 2017

NCBI Genome Database Annotator

@author: Katherine Eaton
"""

import argparse 
import sqlite3
import datetime 
import os

from genomeerrors import *
from genomeutilities import os_check

#-----------------------------------------------------------------------#
#                            Argument Parsing                           #
#-----------------------------------------------------------------------#

parser = argparse.ArgumentParser(description='Description of the Genome Collector Annotation Tool.',
                                 add_help=True,
                                 version = 'GenomeCollector 2.0')

mandatory = parser.add_argument_group('mandatory')
bonus = parser.add_argument_group('bonus')

mandatory.add_argument('--database',
                    help='Genome Database Name',
                    type = str,
                    action = 'store',
                    dest = 'dbname',
                    required=True)

mandatory.add_argument('--annotfile',
                    help='Annotation File',
                    type = str,
                    action = 'store',
                    dest = 'annotfile',
                    required=True)

mandatory.add_argument('--organism',
                    help='Organism Name',
                    type = str,
                    action = 'store',
                    dest = 'orgname',
                    required=True)

args = vars(parser.parse_args())


dbName = args['dbname'] + ".sqlite"
annotFileName = args['annotfile']
orgName = args['orgname']

#-----------------------------------------------------------------------#
#                           Argument Checking                           #
#-----------------------------------------------------------------------#

#---------------------------Check Database------------------------------#

if os.path.exists(dbName):
    conn = sqlite3.connect(dbName)                                      # Make sure database exists
    print('\nOpening database: ' + dbName)
else:
    raise ErrorDBNotExists(dbName)                                      # If database doesn't exists, raise error

if not os.path.exists(annotFileName):
    raise ErrorAnnotFileNotExists(annotFileName)                        # If annotation file doesn't exists, raise error

cur = conn.cursor()                                                     # Create cursor for commands



#-----------------------------------------------------------------------#
#                                File Setup                             #
#-----------------------------------------------------------------------#

OS_SEP = os_check()                                                     # Retrieve the directory separator by OS

if not os.path.exists("log"):                                           # Check if log directory exists
    os.makedirs("log")

str_annot_log_file = "log" + OS_SEP + orgName.replace(" ", "_") + "_db_annotate.log"

if os.path.exists(str_annot_log_file):
    annotate_log_file = open(str_annot_log_file, "a")                   # Open logfile for appending
else:
    annotate_log_file = open(str_annot_log_file, "w")                   # Open logfile for writing
        

annot_file = open(annotFileName, "r")                                    # Open annotation file for writing


#-------------------------Annotation Dictionary-------------------------#

annot_dict = {}

annot_line = annot_file.readline()                                        # Read in the headers
annot_line = annot_file.readline()                                        # Read in the first data line

while annot_line:
    split_line = annot_line.split("\t")
    
    accession = split_line[13].split(" ")[0]                            # Retrieve accession, get rid of spaces and accessory numbers
    latitude = split_line[11]
    longitude = split_line[10]
    source = split_line[9]
    host = split_line[8]
    year = split_line[7]

    annot_dict[accession] = [latitude, longitude, source, host, year]
    annot_line = annot_file.readline()                                    # Read in the next line  


print(annot_dict[accession])  

#-----------------------------------------------------------------------#
#                   Iterate Through Annotation Dict                     #
#-----------------------------------------------------------------------#

for element in annot_dict:
    if element.startswith("A"):
        if element.startswith("AAOS0"):
            element_short = element[0:5] + "2"   
        else:
            element_short = element[0:5] + "1"                                       # db only stores first 5 digits of WGS information plus a 1
        #------------Retrieve Biosample according to WGS Project------------#
        cur.execute('''
        SELECT biosample FROM Assembly WHERE wgs_project=?''',\
        (element_short,))
        
        element_biosample = cur.fetchone()[0]

        # enter element metadata into the biosample table
        
        
    elif element.startswith("G"):
        #------------Retrieve Biosample according to Accession--------------#
        cur.execute('''
        SELECT accession FROM Assembly WHERE accession=?''',\
        (element,))                                                                      # Check if assembly exists in DB
        element_biosample = cur.fetchone()[0]                                           # 0 if not found, 1 if found
        
    else:
        raise ErrorAccessionNotExistInDB(element)                        # If accession doesn't exists, raise error

    print(element_biosample)                   
    latitude = annot_dict[element][0]
    longitude = annot_dict[element][1]
    source = annot_dict[element][2]
    host = annot_dict[element][3]
    year = annot_dict[element][4]
    
    #------------Annotate Metadata for biosample--------------#

    cur.execute('''
    UPDATE BioSample SET latitude=?,longitude=?,geographic_location=?,host=?,collection_date=? WHERE accession=?''',
                (latitude,
                 longitude,
                 source,
                 host,
                 year,
                 element_biosample,))

    #------------------------WWrite to Logfile-------------------------#
    now = datetime.datetime.now()                                           # Get timestamp
    annotate_log_file.write("[" + str(now) + "]" +
                 "\t" + "Accession number meta data modified:" +
                 "\t" + element + "." + "\n")                               # Write to logfile
    #conn.commit()                                                           # Commit assembly record changes to database
                
 

#-----------------------------------------------------------------------#
#                                    Cleanup                            #
#-----------------------------------------------------------------------#
conn.commit()                                                                  # Make sure all changes are committed to database
cur.close()                                                                    # Close the database
annotate_log_file.close()                                                      # Close the annotation file
annot_file.close()                                                      # Close the logfile


