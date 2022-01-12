#!/usr/bin/env python3

import os

os.system(
    """
    python -m coverage run -m pytest --cov=ncbimeta --cov-report=xml \
        test/test_errors.py \
        test/test_utilities.py \
        test/test_ncbimeta.py \
        test/test_annotate.py \
        test/test_join.py \
        test/test_export.py
"""
)
