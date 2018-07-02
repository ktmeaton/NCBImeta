"""
NCBI Metadata Database Annotator

@author: Katherine Eaton
"""

import argparse
import sqlite3
import os
import sys

from NCBImeta_Errors import *

def flushprint(message):
    print(message)
    sys.stdout.flush()

#-----------------------------------------------------------------------#
#                            Argument Parsing                           #
#-----------------------------------------------------------------------#

parser = argparse.ArgumentParser(description=("NCBInfect Annotation Tool"),
                                 add_help=True)

mandatory = parser.add_argument_group('mandatory')

mandatory.add_argument('--database',
                    help='Path to the sqlite database generated by NCBImeta.',
                    type = str,
                    action = 'store',
                    dest = 'dbName',
                    required=True)

mandatory.add_argument('--outputdir',
                    help = 'Output directory.',
                    action = 'store',
                    dest = 'outputDir',
                    required = True)

args = vars(parser.parse_args())

db_name = args['dbName']
output_dir = args['outputDir']



#-----------------------------------------------------------------------#
#                           Argument Checking                           #
#-----------------------------------------------------------------------#

# Check if database exists
if os.path.exists(db_name):
    conn = sqlite3.connect(db_name)
    flushprint('\nOpening database: ' + db_name)
else:
    raise ErrorDBNotExists(db_name)

# Check if output dir exists
if not os.path.exists(output_dir):
    raise ErrorOutputDirNotExists(output_dir)

# no errors were raised, safe to connect to db
cur = conn.cursor()


#-----------------------------------------------------------------------#
#                         Process Database                              #
#-----------------------------------------------------------------------#

# Get a list of tables
cur.execute("SELECT name FROM sqlite_master WHERE type='table';")
table_list = cur.fetchall()

# Iterate through each table
for table in table_list:
    # Skip sqlite_sequence default table
    if table[0] == "sqlite_sequence": continue
    # Convert from tuple where second value is empty
    table_name=table[0]

    # Name of output table file
    db_basename = os.path.basename(db_name).split(".")[0]
    table_file_path = os.path.join(output_dir, "", db_basename + "_" + table_name + ".txt")
    table_file = open(table_file_path, "w")

    # get list of column names in Table
    query = "SELECT * FROM {0}".format(table_name)
    cur.execute(query)
    table_col_names = [description[0] for description in cur.description]

    # Retrieve and write the header
    header = ""
    for column in table_col_names:
        header += column + "\t"
    header = header.rstrip("\t") + "\n"
    table_file.write(header)

    # Retrieve and write the lines
    query = "SELECT * FROM {0}".format(table_name)
    for row in cur.execute(query):
        line = ""
        for cell in row:
            line += str(cell) + "\t"
        line = line.rstrip("\t")  + "\n"
        table_file.write(line)

    # Close table file
    table_file.close()

#-----------------------------------------------------------------------#
#                                    Cleanup                            #
#-----------------------------------------------------------------------#

flushprint("Closing database: " + db_name)
cur.close()