#!/usr/bin/env python3
"""
NCBImeta Join Tool - Joins accessory tables to an anchor table.

@author: Katherine Eaton
"""

# -----------------------------------------------------------------------------#
#                         Modules and Packages                                 #
# -----------------------------------------------------------------------------#

import argparse  # Command-line argument parsing
import sqlite3  # Database storage and queries
import os  # Filepath operations

from ncbimeta import NCBImetaUtilities  # NCBImeta helper functions
from ncbimeta import NCBImetaErrors  # NCBImeta Error classes

# -----------------------------------------------------------------------------#
#                            Argument Parsing                                  #
# -----------------------------------------------------------------------------#

parser = argparse.ArgumentParser(
    description=("NCBImeta Join Tool - Joins accessory tables to an anchor table."),
    add_help=True,
)

mandatory = parser.add_argument_group("mandatory")
bonus = parser.add_argument_group("bonus")

mandatory.add_argument(
    "--database",
    help="Path to the sqlite database generated by NCBImeta.",
    type=str,
    action="store",
    dest="dbName",
    required=True,
)

mandatory.add_argument(
    "--anchor",
    help="Table in NCBImeta database to use as an anchor for joining",
    type=str,
    action="store",
    dest="dbAnchor",
    required=True,
)

mandatory.add_argument(
    "--accessory",
    help=(
        "Accessory tables to join to the anchor as "
        + "space separated str list ex. 'BioProject Assembly' "
    ),
    type=str,
    action="store",
    dest="dbAccessory",
    required=True,
)

mandatory.add_argument(
    "--final",
    help="Name of the final table after the join",
    type=str,
    action="store",
    dest="dbFinal",
    required=True,
)

mandatory.add_argument(
    "--unique",
    help=(
        "Column names in anchor table, to search for "
        + "values in accessory tables ex. 'Accession BioProject' "
    ),
    type=str,
    action="store",
    dest="dbUnique",
    required=True,
)

parser.add_argument("--version", action="version", version="%(prog)s v0.6.6dev")


args = vars(parser.parse_args())

db_name = args["dbName"]
db_anchor = args["dbAnchor"]
db_accessory_str = args["dbAccessory"]
db_accessory_list = db_accessory_str.split(" ")
db_final = args["dbFinal"]
db_all_tables = [db_anchor] + db_accessory_list
unique_header_str = args["dbUnique"]
unique_header_list = unique_header_str.split(" ")
DB_VALUE_SEP = ";"

# -----------------------------------------------------------------------------#
#                           Argument Checking                                  #
# -----------------------------------------------------------------------------#

print("START")

# ---------------------------Check Database------------------------------------#

if os.path.exists(db_name):
    conn = sqlite3.connect(db_name)
    print("\nOpening database: " + db_name + "\n", flush=True)
else:
    raise NCBImetaErrors.ErrorDBNotExists(db_name)

# no errors were raised, safe to connect to db
cur = conn.cursor()

# ---------------------------Check Table Names---------------------------------#
# Check the table names for problematic char
table_name_list = [db_anchor] + [db_final] + db_accessory_list
for table_name in table_name_list:
    table_name_sanitize = NCBImetaUtilities.sql_sanitize(table_name)
    if table_name != table_name_sanitize:
        raise NCBImetaErrors.ErrorSQLNameSanitize(table_name, table_name_sanitize)

# Check the column names
for col_name in unique_header_list:
    col_name_sanitize = NCBImetaUtilities.sql_sanitize(col_name)
    if col_name != col_name_sanitize:
        raise NCBImetaErrors.ErrorSQLNameSanitize(col_name, col_name_sanitize)

# ---------------------------Check Tables Existence----------------------------#

if not NCBImetaUtilities.table_exists(cur, db_anchor):
    raise NCBImetaErrors.ErrorTableNotInDB(db_anchor)
for table in db_accessory_list:
    if not NCBImetaUtilities.table_exists(cur, table):
        raise NCBImetaErrors.ErrorTableNotInDB(table)

# -----------------------------------------------------------------------------#
#                                File Setup                                    #
# -----------------------------------------------------------------------------#

# get list of column names in anchor table
cur.execute("""SELECT * FROM {}""".format(db_anchor))
db_col_names = [
    description[0] for description in cur.description if description[0] != "id"
]
anchor_col_names = [description[0] for description in cur.description]

# Check to make sure the unique header is present in the anchor table
for unique_header in unique_header_list:
    if unique_header not in db_col_names:
        print("Column not in DB: " + unique_header + ".", flush=True)
        raise NCBImetaErrors.ErrorEntryNotInDB(unique_header)

# get list of column names in accessory tables
for table in db_accessory_list:
    cur.execute("""SELECT * FROM {}""".format(table))
    for header in [
        description[0] for description in cur.description if description[0] != "id"
    ]:
        db_col_names.append(header)

# Prepare the columns for the final join table
dupl_col_names = set([col for col in db_col_names if db_col_names.count(col) > 1])
if len(dupl_col_names) > 0:
    raise NCBImetaErrors.ErrorColumnsNotUnique(dupl_col_names)

# -----------------------------------------------------------------------------#
#                              Init Join Table                                 #
# -----------------------------------------------------------------------------#

# Create the database with dynamic variables to store joined info
sql_query = (
    "Create TABLE IF NOT EXISTS "
    + db_final
    + " (id INTEGER NOT NULL PRIMARY KEY AUTOINCREMENT UNIQUE "
)

for column_name in db_col_names:
    sql_query += ", " + column_name + " TEXT"
sql_query += ")"

cur.execute(sql_query)

# -----------------------------------------------------------------------------#
#                              Join Tables                                     #
# -----------------------------------------------------------------------------#

# Retrieve all the anchor table records
cur.execute("""SELECT * FROM {}""".format(db_anchor))
fetch_records = cur.fetchall()

num_processed = 0
num_records = len(fetch_records)


# Iterate through each record
for record in fetch_records:
    # Initialize/Reinitialize the master table value dictionary
    master_column_dict = {}
    for col in db_col_names:
        master_column_dict[col] = ""

    # Store the anchor table values
    for i in range(0, len(record)):
        # Skip id
        if anchor_col_names[i] == "id":
            continue
        if record[i] is None:
            record_val = ""
        else:
            # No quote wrapping needed after parameterizing input?
            record_val = str(record[i])
        master_column_dict[anchor_col_names[i]] = record_val

    # Find the first unique column and its headers
    unique_val = master_column_dict[unique_header_list[0]]

    # -------------------Progress Log and Entry Counter------------------------#
    # Increment entry counter and record progress to screen

    num_processed += 1
    print("Unique Value: " + unique_val, flush=True)
    print(
        "Processing record: " + str(num_processed) + "/" + str(num_records), flush=True
    )

    # Check if this record already exists in the master join table
    sql_query = (
        "SELECT EXISTS(SELECT "
        + unique_header_list[0]
        + " FROM "
        + db_final
        + " WHERE "
        + unique_header_list[0]
        + "=?)"
    )

    cur.execute(sql_query, (unique_val,))

    # 0 if not found, 1 if found
    record_exists = cur.fetchone()[0]
    if record_exists:
        continue

    # Now grab the associated records in the accessory tables
    sql_query = """SELECT {0} FROM {1} WHERE {2}=?""".format(
        ",".join(unique_header_list), db_anchor, unique_header_list[0]
    )
    cur.execute(sql_query, (unique_val,))
    fetch_unique_vals_anchor = cur.fetchall()

    for unique_values in fetch_unique_vals_anchor:
        # Convert from tuple to list for index positions
        unique_values = list(unique_values)
        # Ignore any supposedly unique values that wound up being missing
        while None in unique_values:
            unique_values.remove(None)
        # Properly remove empty string values
        while "" in unique_values:
            unique_values.remove("")

        # search for this value in each accessory table
        for table in db_accessory_list:
            # Get list of each column
            cur.execute(""" SELECT * FROM {}""".format(table))
            table_col_names = [description[0] for description in cur.description]

            # Initialize matching variables
            match_found = False
            match_column = ""
            # Iterate through each possible unique values (or until match is found)
            for uniq_val in unique_values:
                val = str(uniq_val)
                # Iterate through each column
                for table_col in table_col_names:
                    cur.execute("""SELECT {0} FROM {1}""".format(table_col, table))
                    table_col_vals = cur.fetchall()
                    # Search through every value for a match
                    for val in table_col_vals:
                        val = str(val[0])

                        # If a match, store the value, and set the boolean flag
                        if val == uniq_val:
                            match_found = True
                            match_column = table_col
                            match_val = val
                            # Found the match, stop searching through this column
                            break

                    # Found a match, stop searching through columns
                    if match_found:
                        break
                # Found  a match, stop searching values
                if match_found:
                    break

            if match_found:
                query = """SELECT * FROM {0} WHERE {1}=?""".format(table, match_column)

                cur.execute(query, (match_val,))
                match_records = cur.fetchall()
                record_dict = {}

                # If multiple records are found, concat into a pseudo single record
                if len(match_records) > 1:
                    match_records_concat = [None] * len(match_records[0])
                    for i in range(0, len(match_records[0])):
                        tmp_record_list = []
                        for record in match_records:
                            if record[i] is None or record[i] == "":
                                continue
                            elif type(record[i]) == int:
                                tmp_record_list.append(str(record[i]))
                            else:
                                tmp_record_list.append(record[i])
                        if len(tmp_record_list) > 1:
                            # Check if all values are identical
                            dupl_values = set(
                                [
                                    val
                                    for val in tmp_record_list
                                    if tmp_record_list.count(val)
                                    == len(tmp_record_list)
                                ]
                            )
                            if len(dupl_values) == 1:
                                match_records_concat[i] = list(dupl_values)[0]
                            else:
                                match_records_concat[i] = DB_VALUE_SEP.join(
                                    tmp_record_list
                                )

                        else:
                            match_records_concat[i] = None
                    match_records = [match_records_concat]

                match_records = list(match_records[0])
                for i in range(0, len(table_col_names)):
                    # Check if this is a valid col name in the master join table
                    if table_col_names[i] not in master_column_dict.keys():
                        continue
                    record_val = match_records[i]
                    # Check for None and handle unicode
                    if record_val is None:
                        record_val = ""
                    record_val = str(record_val)
                    # Assign record to dictionary
                    master_column_dict[table_col_names[i]] = record_val
                    # Uncommented the following line when
                    # excessive semi-colons appeared in between every char
                    master_column_dict[table_col_names[i]] = master_column_dict[
                        table_col_names[i]
                    ]

        # We've now finished processing all accessory tables
        # Add values to new join table

        # Experiment in placeholders to avoid escaping issue of special char
        sql_q_marks = ",".join(["?"] * len(master_column_dict.keys()))
        sql_q_marks = "(" + sql_q_marks + ")"
        sql_dynamic_colnames = "(" + ",".join(master_column_dict.keys()) + ")"
        sql_values_placeholder = [
            master_column_dict[column] for column in master_column_dict.keys()
        ]

        sql_query = (
            "INSERT INTO "
            + db_final
            + " "
            + sql_dynamic_colnames
            + " VALUES "
            + sql_q_marks
        )

        cur.execute(sql_query, sql_values_placeholder)
        # Save Changes
        conn.commit()

# -----------------------------------------------------------------------------#
#                                    Cleanup                                   #
# -----------------------------------------------------------------------------#
# Commit changes
conn.commit()
print("Closing database: " + db_name, flush=True)
cur.close()
