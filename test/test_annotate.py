"""
NCBImeta Test - Annotate Replace

@author: Katherine Eaton
"""

# -----------------------------------------------------------------------------#
#                         Modules and Packages                                 #
# -----------------------------------------------------------------------------#

import os  # Filepath operations
import test.test_ncbimeta as test_ncbimeta  # Run main program to create test db
import subprocess  # Execute CLI/Shell

# -----------------------------------------------------------------------------#
#                           Test Function                                      #
# -----------------------------------------------------------------------------#


def test_annotate_concatenate_run():
    """Test the NCBImetaAnnotate application in concatenate mode for run completion"""
    # Use the test database
    test_db = os.path.join(os.path.dirname(os.path.abspath(__file__)), "test.sqlite")
    test_annotfile = os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "test_annot.txt"
    )
    # If the test_db doesn't alread exist, run the test cmd from test_ncbimeta
    if not os.path.exists(test_db):
        test_ncbimeta.test_ncbimeta_run()
    test_table = "BioSample"
    test_cmd = (
        "ncbimeta/NCBImetaAnnotate --database "
        + test_db
        + " --table  "
        + test_table
        + " --annotfile "
        + test_annotfile
        + " --concatenate"
    )
    # test NCBImetaAnnotate through a subprocess
    returned_value = subprocess.call(test_cmd, shell=True)
    # If it returns a non-zero value, it failed
    assert returned_value == 0


def test_annotate_replace_run():
    """Test the NCBImetaAnnotate application for run completion"""
    # Use the test database
    test_db = os.path.join(os.path.dirname(os.path.abspath(__file__)), "test.sqlite")
    test_annotfile = os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "test_annot.txt"
    )
    # If the test_db doesn't already exist, run the test cmd from test_ncbimeta
    if not os.path.exists(test_db):
        test_ncbimeta.test_ncbimeta_run()
    test_table = "BioSample"
    test_cmd = (
        "ncbimeta/NCBImetaAnnotate --database "
        + test_db
        + " --table  "
        + test_table
        + " --annotfile "
        + test_annotfile
    )
    # test NCBImetaAnnotate through a subprocess
    returned_value = subprocess.call(test_cmd, shell=True)
    # If it returns a non-zero value, it failed
    assert returned_value == 0


# def test_annotate_multipleEntry(tmpdir):
# To Be done

# def test_annotate_ErrorDBNotExists(tmpdir):
# To Be done

# def test_annotate_ErrorAnnotFileNotExists(tmpdir):
# To Be done

# def test_annotate_ErrorTableNotInDB(tmpdir):
# To Be done
