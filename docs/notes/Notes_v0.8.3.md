# v0.8.3

1. Drop support for `Python 3.7`.
1. Test support for `Python 3.10` (conda `numpy` isn't ready).
1. Implement `str` enconding changes that accomdate all versions of `biopython` (>=1.74).
1. Consolidate `pytest` and `pycov` command to `test/test.py`.
1. `Build` and `Example` workflow include more comprehensive tests for `conda`.
1. Create a fallback conda environment for installation (`environment.yaml`).
1. Applied a hotfix for problematic packaging of `lxml` in `conda`.
