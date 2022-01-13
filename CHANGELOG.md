# CHANGELOG

## Development

### Commits

## v0.8.3

### Notes

1. Drop support for `Python 3.7`.
1. Test support for `Python 3.10` (conda `numpy` isn't ready).
1. Implement `str` enconding changes that accomdate all versions of `biopython` (>=1.74).
1. Consolidate `pytest` and `pycov` command to `test/test.py`.
1. `Build` and `Example` workflow include more comprehensive tests for `conda`.
1. Create a fallback conda environment for installation (`environment.yaml`).
1. Applied a hotfix for problematic packaging of `lxml` in `conda`.

### Commits

* [```a3bb443c```](https://github.com/ktmeaton/NCBImeta/commit/a3bb443c) update pre-commit dependencies
* [```e5d6420b```](https://github.com/ktmeaton/NCBImeta/commit/e5d6420b) change v0.8.3dev to v0.8.3
* [```1d45f5bb```](https://github.com/ktmeaton/NCBImeta/commit/1d45f5bb) update recent projects in README
* [```b0fab147```](https://github.com/ktmeaton/NCBImeta/commit/b0fab147) update paths to docs
* [```22d7f0f7```](https://github.com/ktmeaton/NCBImeta/commit/22d7f0f7) update notes for v0.8.3
* [```30bc84e3```](https://github.com/ktmeaton/NCBImeta/commit/30bc84e3) add lxml hotfix to build
* [```3d1e7040```](https://github.com/ktmeaton/NCBImeta/commit/3d1e7040) restore gcc to environment
* [```6f4b3f9c```](https://github.com/ktmeaton/NCBImeta/commit/6f4b3f9c) update install documentation
* [```560b8175```](https://github.com/ktmeaton/NCBImeta/commit/560b8175) add auto yes to pip uninstall
* [```9087056e```](https://github.com/ktmeaton/NCBImeta/commit/9087056e) remove python 3.10, try lxml hotfix
* [```f892e9b9```](https://github.com/ktmeaton/NCBImeta/commit/f892e9b9) fix missing activate
* [```c003f02a```](https://github.com/ktmeaton/NCBImeta/commit/c003f02a) update development docs
* [```9b5c819e```](https://github.com/ktmeaton/NCBImeta/commit/9b5c819e) create named conda environment
* [```de731fef```](https://github.com/ktmeaton/NCBImeta/commit/de731fef) fix missing conda activate
* [```6d3e91fc```](https://github.com/ktmeaton/NCBImeta/commit/6d3e91fc) try to run example from conda using different python versions
* [```636931dd```](https://github.com/ktmeaton/NCBImeta/commit/636931dd) try to install from conda using different python versions
* [```2c9e2c62```](https://github.com/ktmeaton/NCBImeta/commit/2c9e2c62) restrict workflows py 3.10 version to 3.10.1
* [```b970d73a```](https://github.com/ktmeaton/NCBImeta/commit/b970d73a) simplify calls to pytest and pycov
* [```ececd98c```](https://github.com/ktmeaton/NCBImeta/commit/ececd98c) implement functionality for different biopython versions
* [```88061b34```](https://github.com/ktmeaton/NCBImeta/commit/88061b34) specify python 3.10.1 is supported
* [```0acf7a5f```](https://github.com/ktmeaton/NCBImeta/commit/0acf7a5f) relax biopython dependencies to allow most recent package
* [```1abff891```](https://github.com/ktmeaton/NCBImeta/commit/1abff891) delete old annotate replace file
* [```42907aa9```](https://github.com/ktmeaton/NCBImeta/commit/42907aa9) lint test files
* [```801614df```](https://github.com/ktmeaton/NCBImeta/commit/801614df) try putting 3.10 in quotes
* [```709a8817```](https://github.com/ktmeaton/NCBImeta/commit/709a8817) add a full environment for conda install
* [```11024c4c```](https://github.com/ktmeaton/NCBImeta/commit/11024c4c) drop 3.6 add 3.10 in CI
* [```47dd5849```](https://github.com/ktmeaton/NCBImeta/commit/47dd5849) add a stock environment for dev
* [```0828d35e```](https://github.com/ktmeaton/NCBImeta/commit/0828d35e) relax dependencies
* [```d93afb2b```](https://github.com/ktmeaton/NCBImeta/commit/d93afb2b) enable debugging mode
* [```b9f0afdb```](https://github.com/ktmeaton/NCBImeta/commit/b9f0afdb) update sars-cov-2 search term
* [```8ff08d70```](https://github.com/ktmeaton/NCBImeta/commit/8ff08d70) update README with v0.8.2 project as released
* [```f81e8184```](https://github.com/ktmeaton/NCBImeta/commit/f81e8184) note conda install issue

## v0.8.2

### Notes

1. Add quiet mode with `--quiet`.
1. Create `NCBImetaAnnotate` that consolidates `NCBImetaAnnotateReplace` and `NCBImetaAnnotateConcatenate`.
1. Add configs `yersinia_pestis`, `mammuthus`, `hmp_oral`, `ncov_canada`.
1. Remove `sphinx` related dev dependencies.

### Commits

* [```089858ee```](https://github.com/ktmeaton/NCBImeta/commit/089858ee) include notes in CHANGELOG
* [```c1a4b2a1```](https://github.com/ktmeaton/NCBImeta/commit/c1a4b2a1) update notes
* [```aed461b8```](https://github.com/ktmeaton/NCBImeta/commit/aed461b8) update CHANGELOG
* [```d040d1ce```](https://github.com/ktmeaton/NCBImeta/commit/d040d1ce) update README and CHANGELOG
* [```79a35353```](https://github.com/ktmeaton/NCBImeta/commit/79a35353) re add autolots as submodule
* [```e1712b86```](https://github.com/ktmeaton/NCBImeta/commit/e1712b86) remove autologs
* [```372589b2```](https://github.com/ktmeaton/NCBImeta/commit/372589b2) change config outputs to database dir
* [```38869cbc```](https://github.com/ktmeaton/NCBImeta/commit/38869cbc) add bioproject and sra to ncov_canada
* [```c29cab00```](https://github.com/ktmeaton/NCBImeta/commit/c29cab00) update help documentation
* [```653fe168```](https://github.com/ktmeaton/NCBImeta/commit/653fe168) add quiet mode
* [```b11ed919```](https://github.com/ktmeaton/NCBImeta/commit/b11ed919) add execute permissions to new annotate script
* [```72d4b870```](https://github.com/ktmeaton/NCBImeta/commit/72d4b870) add ncov canada config
* [```f2b54b21```](https://github.com/ktmeaton/NCBImeta/commit/f2b54b21) adjust test workflow for renamed annotate
* [```f4a9bd65```](https://github.com/ktmeaton/NCBImeta/commit/f4a9bd65) fix code block in README
* [```741185b1```](https://github.com/ktmeaton/NCBImeta/commit/741185b1) new configs for hmp_oral and yersinia_pestis
* [```e4efa58f```](https://github.com/ktmeaton/NCBImeta/commit/e4efa58f) remove wheel and sphinx related from setup
* [```bdd193f3```](https://github.com/ktmeaton/NCBImeta/commit/bdd193f3) remove wheel and sphinx related dependencies
* [```58c162b5```](https://github.com/ktmeaton/NCBImeta/commit/58c162b5) adjust wheel pos in requirements
* [```14fedacd```](https://github.com/ktmeaton/NCBImeta/commit/14fedacd) try adding wheel to requirements
* [```07f4bb2c```](https://github.com/ktmeaton/NCBImeta/commit/07f4bb2c) adjust workflows for new example and annotation
* [```2e24ec28```](https://github.com/ktmeaton/NCBImeta/commit/2e24ec28) consolidate the annotation tool to a single file
* [```be0d6036```](https://github.com/ktmeaton/NCBImeta/commit/be0d6036) update install channels in build
* [```31ce01ac```](https://github.com/ktmeaton/NCBImeta/commit/31ce01ac) update install instructions
* [```ae039b34```](https://github.com/ktmeaton/NCBImeta/commit/ae039b34) add png versions of paper figures
* [```2706b1a8```](https://github.com/ktmeaton/NCBImeta/commit/2706b1a8) relax dependencies
* [```07bde416```](https://github.com/ktmeaton/NCBImeta/commit/07bde416) try to limit codecov to os and python
* [```55ec78d0```](https://github.com/ktmeaton/NCBImeta/commit/55ec78d0) remove the executable suffix in pypi testing

## v0.8.1

### Notes

1. Add Python 3.9 Support.
1. Bugfix the test export workflow.
1. Update dependncies: lxml, PyYAML, NumPy
1. Disable example workflow.

### Pull Requests

* [```pull/22```](https://github.com/ktmeaton/NCBImeta/pull/22) Add Python 3.9 Support

### Commits

* [```f0ebe0a0```](https://github.com/ktmeaton/NCBImeta/commit/f0ebe0a0) remove dev suffix from version
* [```8a7e17f7```](https://github.com/ktmeaton/NCBImeta/commit/8a7e17f7) update docs
* [```45baf436```](https://github.com/ktmeaton/NCBImeta/commit/45baf436) update autologs with standardized hash
* [```dab96e71```](https://github.com/ktmeaton/NCBImeta/commit/dab96e71) Merge pull request #22 from ktmeaton/dev
* [```dcaca613```](https://github.com/ktmeaton/NCBImeta/commit/dcaca613) add python 3.9 to the install job of build
* [```038d4f3b```](https://github.com/ktmeaton/NCBImeta/commit/038d4f3b) document example workflow disabled
* [```6d491ac5```](https://github.com/ktmeaton/NCBImeta/commit/6d491ac5) relocate codecov config
* [```2c798bca```](https://github.com/ktmeaton/NCBImeta/commit/2c798bca) add config file for codecov
* [```e0d7920f```](https://github.com/ktmeaton/NCBImeta/commit/e0d7920f) update docs before PR
* [```3d79271c```](https://github.com/ktmeaton/NCBImeta/commit/3d79271c) test numpy-1.19.5
* [```3bbbe596```](https://github.com/ktmeaton/NCBImeta/commit/3bbbe596) bugfix in example workflow trigger
* [```05d19602```](https://github.com/ktmeaton/NCBImeta/commit/05d19602) remove print statement
* [```d3028113```](https://github.com/ktmeaton/NCBImeta/commit/d3028113) update dependencies lxml pyaml and numpy
* [```d5eaf85b```](https://github.com/ktmeaton/NCBImeta/commit/d5eaf85b) fix bug in test_export biosample date
* [```4a1d4a7b```](https://github.com/ktmeaton/NCBImeta/commit/4a1d4a7b) add python 3.9 to the testing matrix
* [```31bfcafe```](https://github.com/ktmeaton/NCBImeta/commit/31bfcafe) update ver number to v0.8.1
* [```4272d9ed```](https://github.com/ktmeaton/NCBImeta/commit/4272d9ed) update changelog for v0.8.0

## v0.8.0

### Notes

1. Start new dev branch.
1. Update miniconda actions and use mamba.
1. Update lxml for security vulnerability.
1. Add autologs as a submodule.
1. Create CHANGELOG with autologs.
1. Simplify test config with less fields to check.
1. Remove .py extension from executable scripts.

### Pull Requests

* [```pull/19```](https://github.com/ktmeaton/NCBImeta/pull/19) Security fix for lxml, autologs, and workflow overhaul

### Commits

* [```3b3d630f```](https://github.com/ktmeaton/NCBImeta/commit/3b3d630f) update autologs for branch and tag link
* [```b5ee965a```](https://github.com/ktmeaton/NCBImeta/commit/b5ee965a) update submodules for release
* [```2290cfce```](https://github.com/ktmeaton/NCBImeta/commit/2290cfce) restrict testing workflows to master and dev
* [```3efed457```](https://github.com/ktmeaton/NCBImeta/commit/3efed457) Merge pull request #19 from ktmeaton/dev
* [```b325be9e```](https://github.com/ktmeaton/NCBImeta/commit/b325be9e) disable fail fast and restrict codecov upload
* [```5cb358eb```](https://github.com/ktmeaton/NCBImeta/commit/5cb358eb) workflow overhaul
* [```cc99998e```](https://github.com/ktmeaton/NCBImeta/commit/cc99998e) restrict python versions to >=3.6,<3.9
* [```4eaa9bc6```](https://github.com/ktmeaton/NCBImeta/commit/4eaa9bc6) overhaul test workflow
* [```c174036f```](https://github.com/ktmeaton/NCBImeta/commit/c174036f) overhaul build workflow
* [```cb0ba1f5```](https://github.com/ktmeaton/NCBImeta/commit/cb0ba1f5) add execute permissions to the newly renamed files
* [```0a519e0b```](https://github.com/ktmeaton/NCBImeta/commit/0a519e0b) replace .py extensions for Utilities script
* [```bf4f608b```](https://github.com/ktmeaton/NCBImeta/commit/bf4f608b) update autologs
* [```db06fd2f```](https://github.com/ktmeaton/NCBImeta/commit/db06fd2f) update to v0.8.0
* [```6b858776```](https://github.com/ktmeaton/NCBImeta/commit/6b858776) Remove .py extension from executable scripts
* [```c73cb41b```](https://github.com/ktmeaton/NCBImeta/commit/c73cb41b) Update dev changelog
* [```7c221988```](https://github.com/ktmeaton/NCBImeta/commit/7c221988) simplify database columns for export
* [```9aacd232```](https://github.com/ktmeaton/NCBImeta/commit/9aacd232) add autolog notes for v0.7.0
* [```34ba5495```](https://github.com/ktmeaton/NCBImeta/commit/34ba5495) update lxml for security fixes
* [```2e203c12```](https://github.com/ktmeaton/NCBImeta/commit/2e203c12) update miniconda action and use mamba for macos
* [```78b42297```](https://github.com/ktmeaton/NCBImeta/commit/78b42297) lint new mammuthus config
* [```6dd43d81```](https://github.com/ktmeaton/NCBImeta/commit/6dd43d81) update miniconda action and use mamba
* [```d275545b```](https://github.com/ktmeaton/NCBImeta/commit/d275545b) start new dev branch for v0.7.1
* [```c96c186f```](https://github.com/ktmeaton/NCBImeta/commit/c96c186f) resolve merges before v0.7.1
* [```42239bc3```](https://github.com/ktmeaton/NCBImeta/commit/42239bc3) add mammuthus config
* [```97fcd956```](https://github.com/ktmeaton/NCBImeta/commit/97fcd956) add new logo to README
* [```e466196f```](https://github.com/ktmeaton/NCBImeta/commit/e466196f) add new logo

## v0.7.0

### Notes

1. Allow config parameters to be specified at run-time.
1. Restrict biopython to >=1.74,<1.77 because of Issue #13
1. Specify versions for all user and dev dependencies.
1. Updated and moved PR template.
1. Updated Contributor's Guideline
1. Remove Python 3.5 support because of incompatibility with black.
1. Remove user email and API key from stdout.

### Pull Requests

* [```pull/16```](https://github.com/ktmeaton/NCBImeta/pull/16) Cli Params
* [```pull/14```](https://github.com/ktmeaton/NCBImeta/pull/14) Typeerror

### Commits

* [```18bf0c8d```](https://github.com/ktmeaton/NCBImeta/commit/18bf0c8d) remove dev suffix from version number
* [```fc0a25d4```](https://github.com/ktmeaton/NCBImeta/commit/fc0a25d4) Merge pull request #16 from ktmeaton/cli-params
* [```1592e45e```](https://github.com/ktmeaton/NCBImeta/commit/1592e45e) update help output
* [```28233181```](https://github.com/ktmeaton/NCBImeta/commit/28233181) add tests for cli param and bad api key
* [```2710af1c```](https://github.com/ktmeaton/NCBImeta/commit/2710af1c) cli parameters for NCBI api
* [```b5d48112```](https://github.com/ktmeaton/NCBImeta/commit/b5d48112) reorganize community and additional contributors
* [```34cf2259```](https://github.com/ktmeaton/NCBImeta/commit/34cf2259) update docs with final changes from typeerror
* [```adb51d7f```](https://github.com/ktmeaton/NCBImeta/commit/adb51d7f) Merge pull request #14 from ktmeaton/typeerror
* [```2baf5ba0```](https://github.com/ktmeaton/NCBImeta/commit/2baf5ba0) update the PR checklist
* [```1faf1bc6```](https://github.com/ktmeaton/NCBImeta/commit/1faf1bc6) change ver to 0.7.0dev and restrict all ver
* [```894f428f```](https://github.com/ktmeaton/NCBImeta/commit/894f428f) restrict pytest and cov ver
* [```8c490bcd```](https://github.com/ktmeaton/NCBImeta/commit/8c490bcd) fix lint workflow desc
* [```6dd00337```](https://github.com/ktmeaton/NCBImeta/commit/6dd00337) remove python 3.5 support
* [```d3a5a060```](https://github.com/ktmeaton/NCBImeta/commit/d3a5a060) update changelog for v0.6.7dev
* [```e9ce2d09```](https://github.com/ktmeaton/NCBImeta/commit/e9ce2d09) fix build ver broken link
* [```ffe033b4```](https://github.com/ktmeaton/NCBImeta/commit/ffe033b4) merge pr changes for lint workflow
* [```4146f54f```](https://github.com/ktmeaton/NCBImeta/commit/4146f54f) Merge branch 'dev' into typeerror
* [```5e69d662```](https://github.com/ktmeaton/NCBImeta/commit/5e69d662) change PR template
* [```bf2e903c```](https://github.com/ktmeaton/NCBImeta/commit/bf2e903c) restrict pre-commit to <= 2.6.0
* [```6da67ac9```](https://github.com/ktmeaton/NCBImeta/commit/6da67ac9) reformat workflows and use pip dev install
* [```fa35c1e7```](https://github.com/ktmeaton/NCBImeta/commit/fa35c1e7) restrict biopython ver to <1.77
* [```b6379306```](https://github.com/ktmeaton/NCBImeta/commit/b6379306) update to 0.6.7dev
* [```dca6bb2d```](https://github.com/ktmeaton/NCBImeta/commit/dca6bb2d) remove outdated project ver number
* [```f7d748fc```](https://github.com/ktmeaton/NCBImeta/commit/f7d748fc) add python support and black linting
* [```aa1f84c1```](https://github.com/ktmeaton/NCBImeta/commit/aa1f84c1) update test table assembly values
* [```b7094501```](https://github.com/ktmeaton/NCBImeta/commit/b7094501) add black and flake8 to dev dependencies
* [```41290068```](https://github.com/ktmeaton/NCBImeta/commit/41290068) fix list numbering indentation
* [```f302cce8```](https://github.com/ktmeaton/NCBImeta/commit/f302cce8) add logo asset attribution
* [```5acd4b78```](https://github.com/ktmeaton/NCBImeta/commit/5acd4b78) add draft logo
* [```d9ce0451```](https://github.com/ktmeaton/NCBImeta/commit/d9ce0451) add logo draft
* [```2723e189```](https://github.com/ktmeaton/NCBImeta/commit/2723e189) add url link for ncbi ref
* [```86776275```](https://github.com/ktmeaton/NCBImeta/commit/86776275) note add contributors
* [```d9b47239```](https://github.com/ktmeaton/NCBImeta/commit/d9b47239) conda link

## v0.6.6.post1

### Commits

* [```e931c8be```](https://github.com/ktmeaton/NCBImeta/commit/e931c8be) test post suffix again
* [```c4b682b5```](https://github.com/ktmeaton/NCBImeta/commit/c4b682b5) update to a1 ver
* [```5aafdbc0```](https://github.com/ktmeaton/NCBImeta/commit/5aafdbc0) run example workflow when src py is changed
* [```ffe86dd9```](https://github.com/ktmeaton/NCBImeta/commit/ffe86dd9) just increment the minor ver for bug fix
* [```d4e906b2```](https://github.com/ktmeaton/NCBImeta/commit/d4e906b2) commit to post suffix after test
* [```5be6eb63```](https://github.com/ktmeaton/NCBImeta/commit/5be6eb63) update ver use the post suffix
* [```7103b41b```](https://github.com/ktmeaton/NCBImeta/commit/7103b41b) remove tag hook causes problems for pypi dup
* [```32dcb8ad```](https://github.com/ktmeaton/NCBImeta/commit/32dcb8ad) update ver
* [```cb2997e4```](https://github.com/ktmeaton/NCBImeta/commit/cb2997e4) delete bad position tags
* [```047e96b6```](https://github.com/ktmeaton/NCBImeta/commit/047e96b6) move tag hook under push action
* [```93d6b819```](https://github.com/ktmeaton/NCBImeta/commit/93d6b819) also run workflows on tags
* [```72c05853```](https://github.com/ktmeaton/NCBImeta/commit/72c05853) update to ver v0.6.7a
* [```43747be1```](https://github.com/ktmeaton/NCBImeta/commit/43747be1) update pip and setuptools
* [```6116cdb1```](https://github.com/ktmeaton/NCBImeta/commit/6116cdb1) remove mention of setuptools install
* [```792e28c5```](https://github.com/ktmeaton/NCBImeta/commit/792e28c5) mark project v0.6.6 as RELEASED
* [```9582d78c```](https://github.com/ktmeaton/NCBImeta/commit/9582d78c) make pypi install plain no ver

## v0.6.6

### Pull Requests

* [```pull/12```](https://github.com/ktmeaton/NCBImeta/pull/12) Contributors Guide Test

### Commits

* [```739cad3d```](https://github.com/ktmeaton/NCBImeta/commit/739cad3d) switch release trigger to published
* [```c68e5404```](https://github.com/ktmeaton/NCBImeta/commit/c68e5404) run all workflows on release
* [```ceb5b369```](https://github.com/ktmeaton/NCBImeta/commit/ceb5b369) fix github links
* [```b5681abd```](https://github.com/ktmeaton/NCBImeta/commit/b5681abd) Merge pull request #12 from ktmeaton/dev
* [```1d40f9a9```](https://github.com/ktmeaton/NCBImeta/commit/1d40f9a9) change ver back to v0.6.6
* [```6d99f145```](https://github.com/ktmeaton/NCBImeta/commit/6d99f145) update links
* [```1d939dec```](https://github.com/ktmeaton/NCBImeta/commit/1d939dec) update pypi workflows
* [```6f383e32```](https://github.com/ktmeaton/NCBImeta/commit/6f383e32) update version numbers with dev suffix
* [```1c30f04c```](https://github.com/ktmeaton/NCBImeta/commit/1c30f04c) change version to v0.6.6dev
* [```74dffef0```](https://github.com/ktmeaton/NCBImeta/commit/74dffef0) remove test dir from package
* [```89875fa0```](https://github.com/ktmeaton/NCBImeta/commit/89875fa0) tidy up pull request template
* [```ebceb166```](https://github.com/ktmeaton/NCBImeta/commit/ebceb166) tidy up pull request template
* [```c0cc2e2e```](https://github.com/ktmeaton/NCBImeta/commit/c0cc2e2e) force linters to wait for pre-commit
* [```dcf8bf8e```](https://github.com/ktmeaton/NCBImeta/commit/dcf8bf8e) change flake8 py find file
* [```94762f74```](https://github.com/ktmeaton/NCBImeta/commit/94762f74) trial yaml lint
* [```90e10878```](https://github.com/ktmeaton/NCBImeta/commit/90e10878) contributors guide file update
* [```6add9ae9```](https://github.com/ktmeaton/NCBImeta/commit/6add9ae9) contributors guide file update
* [```f11e68f4```](https://github.com/ktmeaton/NCBImeta/commit/f11e68f4) install notes and project update
* [```5a515e4e```](https://github.com/ktmeaton/NCBImeta/commit/5a515e4e) separate install and run for flake8
* [```764cbac0```](https://github.com/ktmeaton/NCBImeta/commit/764cbac0) add back in flake8
* [```69d8a7cd```](https://github.com/ktmeaton/NCBImeta/commit/69d8a7cd) add back in python linting
* [```d96c0d05```](https://github.com/ktmeaton/NCBImeta/commit/d96c0d05) remove python linting
* [```4ad65ae9```](https://github.com/ktmeaton/NCBImeta/commit/4ad65ae9) remove python linting
* [```445c69aa```](https://github.com/ktmeaton/NCBImeta/commit/445c69aa) remove python linting
* [```8f51169c```](https://github.com/ktmeaton/NCBImeta/commit/8f51169c) just test flake8 install
* [```66376349```](https://github.com/ktmeaton/NCBImeta/commit/66376349) temp remove pr hook and test flake8
* [```522e1243```](https://github.com/ktmeaton/NCBImeta/commit/522e1243) add macos conda install test
* [```56b3ec1d```](https://github.com/ktmeaton/NCBImeta/commit/56b3ec1d) add macos build workflow
* [```2750e066```](https://github.com/ktmeaton/NCBImeta/commit/2750e066) use conda setup action
* [```36b88dad```](https://github.com/ktmeaton/NCBImeta/commit/36b88dad) comment out PyPI repo install for now
* [```48b15d51```](https://github.com/ktmeaton/NCBImeta/commit/48b15d51) remove setuptools install
* [```d44db942```](https://github.com/ktmeaton/NCBImeta/commit/d44db942) MAJOR bug fix for build workflow
* [```86c90653```](https://github.com/ktmeaton/NCBImeta/commit/86c90653) try to force build workflow with python3
* [```d7c1f9f2```](https://github.com/ktmeaton/NCBImeta/commit/d7c1f9f2) test lint workflow no req
* [```9870d239```](https://github.com/ktmeaton/NCBImeta/commit/9870d239) update markdownlint-cli to ignore paper dir
* [```062928c4```](https://github.com/ktmeaton/NCBImeta/commit/062928c4) lint setup
* [```c684f97c```](https://github.com/ktmeaton/NCBImeta/commit/c684f97c) lint bug_report
* [```f7a62b57```](https://github.com/ktmeaton/NCBImeta/commit/f7a62b57) new pull request template
* [```9cdc4ce6```](https://github.com/ktmeaton/NCBImeta/commit/9cdc4ce6) more detailed contibute guide
* [```0c1fcd16```](https://github.com/ktmeaton/NCBImeta/commit/0c1fcd16) more paths to lint
* [```12078ef7```](https://github.com/ktmeaton/NCBImeta/commit/12078ef7) lint conftest
* [```4976309a```](https://github.com/ktmeaton/NCBImeta/commit/4976309a) lint schema README
* [```cd724f9e```](https://github.com/ktmeaton/NCBImeta/commit/cd724f9e) lint config README
* [```93035ff7```](https://github.com/ktmeaton/NCBImeta/commit/93035ff7) update linting workflow
* [```aab8b906```](https://github.com/ktmeaton/NCBImeta/commit/aab8b906) lint conftest max 80
* [```3c4f9a8b```](https://github.com/ktmeaton/NCBImeta/commit/3c4f9a8b) lint conftest max 80
* [```7d6217db```](https://github.com/ktmeaton/NCBImeta/commit/7d6217db) lint test_ncbimeta max 80
* [```f8952f59```](https://github.com/ktmeaton/NCBImeta/commit/f8952f59) lint test_ncbimeta max 80
* [```267b401d```](https://github.com/ktmeaton/NCBImeta/commit/267b401d) lint test_xml max 80
* [```534f3971```](https://github.com/ktmeaton/NCBImeta/commit/534f3971) lint test_utilities max 80
* [```5cc8065f```](https://github.com/ktmeaton/NCBImeta/commit/5cc8065f) lint test_utilities max 80
* [```499c75f3```](https://github.com/ktmeaton/NCBImeta/commit/499c75f3) lint test_join max 80
* [```1dd3536f```](https://github.com/ktmeaton/NCBImeta/commit/1dd3536f) lint test_annotate max 80
* [```06938061```](https://github.com/ktmeaton/NCBImeta/commit/06938061) lint test_export max 80
* [```c8c67cca```](https://github.com/ktmeaton/NCBImeta/commit/c8c67cca) lint test_errors
* [```db12ca1f```](https://github.com/ktmeaton/NCBImeta/commit/db12ca1f) lint NCBImetaAnnotate max 80
* [```b8030efb```](https://github.com/ktmeaton/NCBImeta/commit/b8030efb) lint NCBImetaJoin max 80
* [```6b78c19c```](https://github.com/ktmeaton/NCBImeta/commit/6b78c19c) lint NCBImetaJoin max 80
* [```49a4e918```](https://github.com/ktmeaton/NCBImeta/commit/49a4e918) lint NCBImetaErrors max 80
* [```8e511183```](https://github.com/ktmeaton/NCBImeta/commit/8e511183) lint NCBImetaErrors max 80
* [```1c909129```](https://github.com/ktmeaton/NCBImeta/commit/1c909129) redo py linting max 80
* [```3a912e59```](https://github.com/ktmeaton/NCBImeta/commit/3a912e59) lint NCBImetaUtilities max 80
* [```69abba67```](https://github.com/ktmeaton/NCBImeta/commit/69abba67) lint NCBImetaUtilities max 80
* [```c635f83f```](https://github.com/ktmeaton/NCBImeta/commit/c635f83f) lint NCBImetaUtilities max 80
* [```f0818fd5```](https://github.com/ktmeaton/NCBImeta/commit/f0818fd5) lint NCBImetaUtilities max 80
* [```4997e76f```](https://github.com/ktmeaton/NCBImeta/commit/4997e76f) lint NCBImetaUtilities max 80
* [```d2c2f6d6```](https://github.com/ktmeaton/NCBImeta/commit/d2c2f6d6) lint NCBImetaExport max 80
* [```6cc3f06d```](https://github.com/ktmeaton/NCBImeta/commit/6cc3f06d) lint NCBImetaExport max 80
* [```8f3f898d```](https://github.com/ktmeaton/NCBImeta/commit/8f3f898d) update badges and requirements
* [```54327273```](https://github.com/ktmeaton/NCBImeta/commit/54327273) more informative names for build jobs
* [```229fc561```](https://github.com/ktmeaton/NCBImeta/commit/229fc561) lint CHANGELOG with line len
* [```35a6482a```](https://github.com/ktmeaton/NCBImeta/commit/35a6482a) update build badge to gh actions
* [```64f51457```](https://github.com/ktmeaton/NCBImeta/commit/64f51457) simplify workflow names
* [```96259597```](https://github.com/ktmeaton/NCBImeta/commit/96259597) simplify workflow names
* [```8bf6dbea```](https://github.com/ktmeaton/NCBImeta/commit/8bf6dbea) lint CHANGELOG
* [```3158f783```](https://github.com/ktmeaton/NCBImeta/commit/3158f783) add markdownlint-cli to pre-commit lint
* [```863aaf23```](https://github.com/ktmeaton/NCBImeta/commit/863aaf23) ncbimeta adhere to B950 max 80
* [```f6e58b2f```](https://github.com/ktmeaton/NCBImeta/commit/f6e58b2f) ignore B950 in conftest
* [```a609a1e7```](https://github.com/ktmeaton/NCBImeta/commit/a609a1e7) fix i var in loop
* [```212fd6e4```](https://github.com/ktmeaton/NCBImeta/commit/212fd6e4) revert to len 80
* [```eb1b9f0f```](https://github.com/ktmeaton/NCBImeta/commit/eb1b9f0f) revert to len 80
* [```a3c2b937```](https://github.com/ktmeaton/NCBImeta/commit/a3c2b937) lint ncbimeta with bugbear
* [```abdb2eb5```](https://github.com/ktmeaton/NCBImeta/commit/abdb2eb5) lint ncbimeta with bugbear
* [```b223578b```](https://github.com/ktmeaton/NCBImeta/commit/b223578b) lint ncbimeta with bugbear
* [```f60b4258```](https://github.com/ktmeaton/NCBImeta/commit/f60b4258) lint ncbimeta with bugbear
* [```ec29ac9d```](https://github.com/ktmeaton/NCBImeta/commit/ec29ac9d) lint ncbimeta with bugbear
* [```0e1b5198```](https://github.com/ktmeaton/NCBImeta/commit/0e1b5198) lint ncbimeta
* [```a790dd61```](https://github.com/ktmeaton/NCBImeta/commit/a790dd61) attempte bugbear for flake8
* [```3e93d0f4```](https://github.com/ktmeaton/NCBImeta/commit/3e93d0f4) more str format for black
* [```37beb315```](https://github.com/ktmeaton/NCBImeta/commit/37beb315) reformat str to make black happy
* [```caf0312c```](https://github.com/ktmeaton/NCBImeta/commit/caf0312c) lint with new rules
* [```15b5e41f```](https://github.com/ktmeaton/NCBImeta/commit/15b5e41f) lint with new rules
* [```a7a8e327```](https://github.com/ktmeaton/NCBImeta/commit/a7a8e327) lint with new rules
* [```e63ddc95```](https://github.com/ktmeaton/NCBImeta/commit/e63ddc95) lint with new rules
* [```048855c4```](https://github.com/ktmeaton/NCBImeta/commit/048855c4) lint with new rules
* [```f69c9bfd```](https://github.com/ktmeaton/NCBImeta/commit/f69c9bfd) lint with new rules
* [```9ccea9c9```](https://github.com/ktmeaton/NCBImeta/commit/9ccea9c9) fix test module import
* [```55da36c0```](https://github.com/ktmeaton/NCBImeta/commit/55da36c0) lint conftest
* [```e41452e7```](https://github.com/ktmeaton/NCBImeta/commit/e41452e7) test file specific rule exclude
* [```cad97a45```](https://github.com/ktmeaton/NCBImeta/commit/cad97a45) change E501 exclude to conftest
* [```294bc50b```](https://github.com/ktmeaton/NCBImeta/commit/294bc50b) flake8 config in setup.cfg
* [```ab237975```](https://github.com/ktmeaton/NCBImeta/commit/ab237975) test flake8 ignore conftest
* [```4563f59a```](https://github.com/ktmeaton/NCBImeta/commit/4563f59a) lint test_ncbimeta
* [```9fa14aa5```](https://github.com/ktmeaton/NCBImeta/commit/9fa14aa5) lint test_ncbimeta
* [```b22f4915```](https://github.com/ktmeaton/NCBImeta/commit/b22f4915) lint test_utilities
* [```2be38960```](https://github.com/ktmeaton/NCBImeta/commit/2be38960) lint test_utilities
* [```082f9156```](https://github.com/ktmeaton/NCBImeta/commit/082f9156) lint test_join
* [```b5742fd9```](https://github.com/ktmeaton/NCBImeta/commit/b5742fd9) lint test_export
* [```7d4f29c1```](https://github.com/ktmeaton/NCBImeta/commit/7d4f29c1) lint test_join
* [```efa7c750```](https://github.com/ktmeaton/NCBImeta/commit/efa7c750) lint test_export
* [```93cffa9f```](https://github.com/ktmeaton/NCBImeta/commit/93cffa9f) lint test_errors
* [```c25df97f```](https://github.com/ktmeaton/NCBImeta/commit/c25df97f) lint test_errors
* [```cde103fe```](https://github.com/ktmeaton/NCBImeta/commit/cde103fe) fix comment spacing
* [```5d574303```](https://github.com/ktmeaton/NCBImeta/commit/5d574303) fix comment spacing
* [```25ee239a```](https://github.com/ktmeaton/NCBImeta/commit/25ee239a) lint test_xml
* [```000d9060```](https://github.com/ktmeaton/NCBImeta/commit/000d9060) lint test_xml
* [```810a70cc```](https://github.com/ktmeaton/NCBImeta/commit/810a70cc) lint test_annotateconcatenate
* [```c76738a5```](https://github.com/ktmeaton/NCBImeta/commit/c76738a5) lint test_annotateconcatenate
* [```04dd0992```](https://github.com/ktmeaton/NCBImeta/commit/04dd0992) lint test_annotatereplace
* [```ce0f9e0d```](https://github.com/ktmeaton/NCBImeta/commit/ce0f9e0d) lint test_annotatereplace
* [```a6acbc0a```](https://github.com/ktmeaton/NCBImeta/commit/a6acbc0a) lint test_annotatereplace
* [```f039442c```](https://github.com/ktmeaton/NCBImeta/commit/f039442c) lint NCBImeta
* [```489fc1b5```](https://github.com/ktmeaton/NCBImeta/commit/489fc1b5) lint NCBImeta
* [```cf980d10```](https://github.com/ktmeaton/NCBImeta/commit/cf980d10) lint NCBImetaAnnotate
* [```2829ebfe```](https://github.com/ktmeaton/NCBImeta/commit/2829ebfe) lint NCBImetaAnnotate
* [```8388eb77```](https://github.com/ktmeaton/NCBImeta/commit/8388eb77) lint NCBImetaAnnotateReplace
* [```4ccfc225```](https://github.com/ktmeaton/NCBImeta/commit/4ccfc225) lint NCBImetaExport
* [```c2aab367```](https://github.com/ktmeaton/NCBImeta/commit/c2aab367) lint NCBImetaExport
* [```b00abcd0```](https://github.com/ktmeaton/NCBImeta/commit/b00abcd0) add lint checklist for py files
* [```cc8243df```](https://github.com/ktmeaton/NCBImeta/commit/cc8243df) lint NCBImetaJoin
* [```7614502d```](https://github.com/ktmeaton/NCBImeta/commit/7614502d) lint single quotes
* [```3fe54360```](https://github.com/ktmeaton/NCBImeta/commit/3fe54360) pull request template
* [```6d2c1b7b```](https://github.com/ktmeaton/NCBImeta/commit/6d2c1b7b) remove old dev require file
* [```58b0669c```](https://github.com/ktmeaton/NCBImeta/commit/58b0669c) draft contributors guide
* [```6fdf6f2e```](https://github.com/ktmeaton/NCBImeta/commit/6fdf6f2e) update build comments
* [```348dd6df```](https://github.com/ktmeaton/NCBImeta/commit/348dd6df) linting and dev depend update
* [```854a4371```](https://github.com/ktmeaton/NCBImeta/commit/854a4371) remove old contributor notes
* [```ce36f022```](https://github.com/ktmeaton/NCBImeta/commit/ce36f022) remove poetry toml
* [```db9889b6```](https://github.com/ktmeaton/NCBImeta/commit/db9889b6) update Lat and Lon for BioSample schema yaml
* [```e3259858```](https://github.com/ktmeaton/NCBImeta/commit/e3259858) notes on replacing old schema
* [```ebaf5eb2```](https://github.com/ktmeaton/NCBImeta/commit/ebaf5eb2) dev testing toml file
* [```e83d45c9```](https://github.com/ktmeaton/NCBImeta/commit/e83d45c9) formatting notes
* [```4c9de0a2```](https://github.com/ktmeaton/NCBImeta/commit/4c9de0a2) fix line wrap
* [```ff0bb82e```](https://github.com/ktmeaton/NCBImeta/commit/ff0bb82e) add python file trigger
* [```cb07360a```](https://github.com/ktmeaton/NCBImeta/commit/cb07360a) black and flake8 format
* [```1b1454f5```](https://github.com/ktmeaton/NCBImeta/commit/1b1454f5) poetry experiment
* [```3f1e9117```](https://github.com/ktmeaton/NCBImeta/commit/3f1e9117) force python 3.7
* [```dc9da600```](https://github.com/ktmeaton/NCBImeta/commit/dc9da600) rename dev requirements
* [```323627ea```](https://github.com/ktmeaton/NCBImeta/commit/323627ea) basic deploy workflow for pypi
* [```b9d387a0```](https://github.com/ktmeaton/NCBImeta/commit/b9d387a0) fix indendation
* [```506d2c30```](https://github.com/ktmeaton/NCBImeta/commit/506d2c30) create deploy.yaml
* [```22652e0a```](https://github.com/ktmeaton/NCBImeta/commit/22652e0a) line endings fix
* [```b0729285```](https://github.com/ktmeaton/NCBImeta/commit/b0729285) developer dependencies install guide
* [```a35875fe```](https://github.com/ktmeaton/NCBImeta/commit/a35875fe) fix mixed line endings
* [```5b32ea92```](https://github.com/ktmeaton/NCBImeta/commit/5b32ea92) pre-commit configuration
* [```52e9e799```](https://github.com/ktmeaton/NCBImeta/commit/52e9e799) remove trailing whitespace and file endings
* [```094ff71f```](https://github.com/ktmeaton/NCBImeta/commit/094ff71f) gh actions and linting
* [```1cc59ead```](https://github.com/ktmeaton/NCBImeta/commit/1cc59ead) fix indentation
* [```5bdb96ee```](https://github.com/ktmeaton/NCBImeta/commit/5bdb96ee) correct lint extension
* [```0dc356f9```](https://github.com/ktmeaton/NCBImeta/commit/0dc356f9) restrict paths and add more workflow
* [```6a5590c6```](https://github.com/ktmeaton/NCBImeta/commit/6a5590c6) remove codecov command
* [```58d10a9c```](https://github.com/ktmeaton/NCBImeta/commit/58d10a9c) update all names
* [```e6485616```](https://github.com/ktmeaton/NCBImeta/commit/e6485616) path in quotes
* [```8abaea04```](https://github.com/ktmeaton/NCBImeta/commit/8abaea04) linting yaml with markdown
* [```d1a446a2```](https://github.com/ktmeaton/NCBImeta/commit/d1a446a2) add workflow to path
* [```fce336df```](https://github.com/ktmeaton/NCBImeta/commit/fce336df) wildcard paths in quotes
* [```e5332c86```](https://github.com/ktmeaton/NCBImeta/commit/e5332c86) all branches restrict paths
* [```ff2e1c62```](https://github.com/ktmeaton/NCBImeta/commit/ff2e1c62) restrict paths remove example
* [```76d1a77b```](https://github.com/ktmeaton/NCBImeta/commit/76d1a77b) change repo name to action
* [```ba746033```](https://github.com/ktmeaton/NCBImeta/commit/ba746033) remove version number
* [```c2d46a19```](https://github.com/ktmeaton/NCBImeta/commit/c2d46a19) change repo tarball ver
* [```906fc27a```](https://github.com/ktmeaton/NCBImeta/commit/906fc27a) quick start example
* [```e5e516c1```](https://github.com/ktmeaton/NCBImeta/commit/e5e516c1) codecov upload
* [```01b938f2```](https://github.com/ktmeaton/NCBImeta/commit/01b938f2) fix indentation
* [```c65775fb```](https://github.com/ktmeaton/NCBImeta/commit/c65775fb) ensure unique names
* [```e49c4b11```](https://github.com/ktmeaton/NCBImeta/commit/e49c4b11) formatting and new test
* [```49a80089```](https://github.com/ktmeaton/NCBImeta/commit/49a80089) rename to build
* [```07f4bc79```](https://github.com/ktmeaton/NCBImeta/commit/07f4bc79) add main installation
* [```744a0d4c```](https://github.com/ktmeaton/NCBImeta/commit/744a0d4c) use the pip tutorial gh actions
* [```d6c38cd1```](https://github.com/ktmeaton/NCBImeta/commit/d6c38cd1) sudo for pip installs
* [```b31a524c```](https://github.com/ktmeaton/NCBImeta/commit/b31a524c) try sudo instead
* [```3bec2315```](https://github.com/ktmeaton/NCBImeta/commit/3bec2315) remember to activate
* [```d37982a8```](https://github.com/ktmeaton/NCBImeta/commit/d37982a8) try in conda env
* [```eb7ea002```](https://github.com/ktmeaton/NCBImeta/commit/eb7ea002) Try github actions
* [```16799dcb```](https://github.com/ktmeaton/NCBImeta/commit/16799dcb) test out github actions
* [```3cc15721```](https://github.com/ktmeaton/NCBImeta/commit/3cc15721) fix line ending
* [```0e901059```](https://github.com/ktmeaton/NCBImeta/commit/0e901059) MANIFEST for pypi packaging
* [```228c3698```](https://github.com/ktmeaton/NCBImeta/commit/228c3698) update v0.6.5 to v0.6.6
* [```7f92720f```](https://github.com/ktmeaton/NCBImeta/commit/7f92720f) minor update
* [```77403356```](https://github.com/ktmeaton/NCBImeta/commit/77403356) config for manifest
* [```5c4479bf```](https://github.com/ktmeaton/NCBImeta/commit/5c4479bf) document dev dependencies
* [```42ba5eac```](https://github.com/ktmeaton/NCBImeta/commit/42ba5eac) ignore fusion table
* [```3c9abf3b```](https://github.com/ktmeaton/NCBImeta/commit/3c9abf3b) Convert BioProjectTitle to XPath query
* [```a150cd64```](https://github.com/ktmeaton/NCBImeta/commit/a150cd64) SQL Update fix for char escape and security
* [```bb52506a```](https://github.com/ktmeaton/NCBImeta/commit/bb52506a) fix slack PR notification bool
* [```e55544e3```](https://github.com/ktmeaton/NCBImeta/commit/e55544e3) temp pip instructions
* [```3a88be83```](https://github.com/ktmeaton/NCBImeta/commit/3a88be83) pypi install note

## v0.6.5

### Pull Requests

* [```pull/9```](https://github.com/ktmeaton/NCBImeta/pull/9) Add ability to use XPath queries from within the .yaml config

### Commits

* [```ccfabe3a```](https://github.com/ktmeaton/NCBImeta/commit/ccfabe3a) Merge branch 'master' into dev
* [```4531fdf6```](https://github.com/ktmeaton/NCBImeta/commit/4531fdf6) update dates, ver, community guidelines
* [```a455008c```](https://github.com/ktmeaton/NCBImeta/commit/a455008c) Move XPath query explanation to README_schema doc
* [```b61f931d```](https://github.com/ktmeaton/NCBImeta/commit/b61f931d) update version number to v0.6.5
* [```138cc5f7```](https://github.com/ktmeaton/NCBImeta/commit/138cc5f7) Merge pull request #9 from hellothisisMatt/feature/add-xpath-support
* [```a7415066```](https://github.com/ktmeaton/NCBImeta/commit/a7415066) Update all other issues template [skip ci]
* [```7d375742```](https://github.com/ktmeaton/NCBImeta/commit/7d375742) Update feature request template
* [```86ccf013```](https://github.com/ktmeaton/NCBImeta/commit/86ccf013) Update bug report templates
* [```66dfb737```](https://github.com/ktmeaton/NCBImeta/commit/66dfb737) remove branch pr9 from travis testing
* [```fd166157```](https://github.com/ktmeaton/NCBImeta/commit/fd166157) update unittests to include new values retrieved with advanced Xpath
* [```564c18fb```](https://github.com/ktmeaton/NCBImeta/commit/564c18fb) test implementation and error class when Xpath query is empty and unspecified
* [```9b0b11fd```](https://github.com/ktmeaton/NCBImeta/commit/9b0b11fd) define an error class for when the Xpath query is empty and unspecified
* [```ade807c0```](https://github.com/ktmeaton/NCBImeta/commit/ade807c0) check if the Xpath query was left empty and unspecified
* [```2c5e9a20```](https://github.com/ktmeaton/NCBImeta/commit/2c5e9a20) add XPath queries to match the example config
* [```35b7b950```](https://github.com/ktmeaton/NCBImeta/commit/35b7b950) Merge pull request #1 from ktmeaton/pr9
* [```659238ae```](https://github.com/ktmeaton/NCBImeta/commit/659238ae) update changelog with error classes, run CI
* [```c8b6584f```](https://github.com/ktmeaton/NCBImeta/commit/c8b6584f) New error classes that can be raised in adv_xml_search
* [```c04e981d```](https://github.com/ktmeaton/NCBImeta/commit/c04e981d) revise adv_xml_search to check type of result
* [```1817705c```](https://github.com/ktmeaton/NCBImeta/commit/1817705c) remove outdated test xpath commands
* [```a5c2fff6```](https://github.com/ktmeaton/NCBImeta/commit/a5c2fff6) fix non-specific BioSampleBioProjectAccession
* [```c390ff66```](https://github.com/ktmeaton/NCBImeta/commit/c390ff66) refine xpath search for attr value
* [```84180028```](https://github.com/ktmeaton/NCBImeta/commit/84180028) Enable attribute returns and error checking
* [```a9b80573```](https://github.com/ktmeaton/NCBImeta/commit/a9b80573) enable dev and py9 for CI
* [```ddc6a065```](https://github.com/ktmeaton/NCBImeta/commit/ddc6a065) correct NucleotideBioSampleAccession xpath query
* [```8f0e2ead```](https://github.com/ktmeaton/NCBImeta/commit/8f0e2ead) simplified test example for tip search
* [```0a5fb699```](https://github.com/ktmeaton/NCBImeta/commit/0a5fb699) add contrib info
* [```bffe8294```](https://github.com/ktmeaton/NCBImeta/commit/bffe8294) XPATH use, remove nonspecific NucleotideAssemblyAccession
* [```16451e7d```](https://github.com/ktmeaton/NCBImeta/commit/16451e7d) clarify comma separated with space
* [```6e972c21```](https://github.com/ktmeaton/NCBImeta/commit/6e972c21) 2 tests of adv_xml_search starting from root or starting from tip
* [```51584860```](https://github.com/ktmeaton/NCBImeta/commit/51584860) xpath param of adv_xml_search as preformatted xml
* [```c22c8596```](https://github.com/ktmeaton/NCBImeta/commit/c22c8596) document the space and comma
* [```65353dca```](https://github.com/ktmeaton/NCBImeta/commit/65353dca) Merge upstream dev changes into pr9
* [```55859648```](https://github.com/ktmeaton/NCBImeta/commit/55859648) Add XPATH Information to config readme
* [```56e0b3c5```](https://github.com/ktmeaton/NCBImeta/commit/56e0b3c5) Add ability to use full XPath for XML searching
* [```3206959a```](https://github.com/ktmeaton/NCBImeta/commit/3206959a) new to do fixes
* [```62ab0e7c```](https://github.com/ktmeaton/NCBImeta/commit/62ab0e7c) SRABioProjectAccession edge cases
* [```0cadcec0```](https://github.com/ktmeaton/NCBImeta/commit/0cadcec0) SRABioSampleAccession typo
* [```5232f359```](https://github.com/ktmeaton/NCBImeta/commit/5232f359) update SRABioProjectAccession
* [```24bdfe01```](https://github.com/ktmeaton/NCBImeta/commit/24bdfe01) v0.6.5 init changes
* [```cd5efb4d```](https://github.com/ktmeaton/NCBImeta/commit/cd5efb4d) update BioSample schema and example for SRABioSampleAccession
* [```fb579507```](https://github.com/ktmeaton/NCBImeta/commit/fb579507) badge update and rearrange [skip ci]

## v0.6.4

### Commits

* [```33cdca3f```](https://github.com/ktmeaton/NCBImeta/commit/33cdca3f) switch travis back to master only
* [```03481fe5```](https://github.com/ktmeaton/NCBImeta/commit/03481fe5) expand issues section [skip ci]
* [```1b94dfeb```](https://github.com/ktmeaton/NCBImeta/commit/1b94dfeb) upcoming features section [skip ci]
* [```a55df18f```](https://github.com/ktmeaton/NCBImeta/commit/a55df18f) doc update [skip ci]
* [```80889436```](https://github.com/ktmeaton/NCBImeta/commit/80889436) test multi match
* [```7e5cabd9```](https://github.com/ktmeaton/NCBImeta/commit/7e5cabd9) update nucleotide data
* [```575b153c```](https://github.com/ktmeaton/NCBImeta/commit/575b153c) switch travis to dev branch
* [```9e56cec6```](https://github.com/ktmeaton/NCBImeta/commit/9e56cec6) small progress on nucleotide data
* [```a81357e4```](https://github.com/ktmeaton/NCBImeta/commit/a81357e4) v0.6.4 notes
* [```c7dd8faa```](https://github.com/ktmeaton/NCBImeta/commit/c7dd8faa) ignore some example and test files
* [```d3fa6ccd```](https://github.com/ktmeaton/NCBImeta/commit/d3fa6ccd) ver update
* [```a4a3d9ee```](https://github.com/ktmeaton/NCBImeta/commit/a4a3d9ee) commit to encode removal
* [```d3abba52```](https://github.com/ktmeaton/NCBImeta/commit/d3abba52) working on nucleotide export
* [```300055a7```](https://github.com/ktmeaton/NCBImeta/commit/300055a7) encode decode purge
* [```cdf506f8```](https://github.com/ktmeaton/NCBImeta/commit/cdf506f8) uncomment encode decode calls
* [```33192faf```](https://github.com/ktmeaton/NCBImeta/commit/33192faf) SQL parameter and table col check
* [```31f92037```](https://github.com/ktmeaton/NCBImeta/commit/31f92037) Extra param checking Issue #7
* [```585a200a```](https://github.com/ktmeaton/NCBImeta/commit/585a200a) table and column name catching
* [```5b46aaa0```](https://github.com/ktmeaton/NCBImeta/commit/5b46aaa0) SQL parameter complete
* [```b11f401b```](https://github.com/ktmeaton/NCBImeta/commit/b11f401b) sanitize output format again
* [```4eca45d9```](https://github.com/ktmeaton/NCBImeta/commit/4eca45d9) sanitize output format
* [```69b1f769```](https://github.com/ktmeaton/NCBImeta/commit/69b1f769) test for table and column name format
* [```8104b940```](https://github.com/ktmeaton/NCBImeta/commit/8104b940) SQL table name add
* [```73af61a4```](https://github.com/ktmeaton/NCBImeta/commit/73af61a4) before more sql select change
* [```e4c17c70```](https://github.com/ktmeaton/NCBImeta/commit/e4c17c70) slack travis-ci integration

## v0.6.3

### Commits

* [```77e6fe17```](https://github.com/ktmeaton/NCBImeta/commit/77e6fe17) enable master branch CI
* [```8850b76b```](https://github.com/ktmeaton/NCBImeta/commit/8850b76b) codecov spacing syntax error
* [```523687be```](https://github.com/ktmeaton/NCBImeta/commit/523687be) allow dev CI
* [```16a4f56f```](https://github.com/ktmeaton/NCBImeta/commit/16a4f56f) bioconda autobump bot rely
* [```4446f94f```](https://github.com/ktmeaton/NCBImeta/commit/4446f94f) remove bioconda scripts
* [```b92681f4```](https://github.com/ktmeaton/NCBImeta/commit/b92681f4) autobump script update
* [```49e4bab7```](https://github.com/ktmeaton/NCBImeta/commit/49e4bab7) hidden script names and autobump script [skip ci]
* [```502338b1```](https://github.com/ktmeaton/NCBImeta/commit/502338b1) resume full testing
* [```a76ff011```](https://github.com/ktmeaton/NCBImeta/commit/a76ff011) use the bioconda recommended install code [skip ci]
* [```052a82c5```](https://github.com/ktmeaton/NCBImeta/commit/052a82c5) solved mystery of branch bc-tbd [skip ci]
* [```324b668e```](https://github.com/ktmeaton/NCBImeta/commit/324b668e) automate ver number update [skip ci]
* [```41ed42ed```](https://github.com/ktmeaton/NCBImeta/commit/41ed42ed) switch zenodo citation to the all ver DOI [skip ci]
* [```c00f4784```](https://github.com/ktmeaton/NCBImeta/commit/c00f4784) remove unnecessary tag condition
* [```2d5dda85```](https://github.com/ktmeaton/NCBImeta/commit/2d5dda85) if statement spacing
* [```f139b226```](https://github.com/ktmeaton/NCBImeta/commit/f139b226) simplify commit message [skip ci]
* [```d6c3c027```](https://github.com/ktmeaton/NCBImeta/commit/d6c3c027) travis tag match v
* [```553124d1```](https://github.com/ktmeaton/NCBImeta/commit/553124d1) tag check
* [```fb349cad```](https://github.com/ktmeaton/NCBImeta/commit/fb349cad) file format linux
* [```086851a2```](https://github.com/ktmeaton/NCBImeta/commit/086851a2) executable mode
* [```4053814a```](https://github.com/ktmeaton/NCBImeta/commit/4053814a) move conda update commands to script
* [```23292884```](https://github.com/ktmeaton/NCBImeta/commit/23292884) set branch upstream origin
* [```6b3d224f```](https://github.com/ktmeaton/NCBImeta/commit/6b3d224f) add, commit, push recipe
* [```20b58f61```](https://github.com/ktmeaton/NCBImeta/commit/20b58f61) enforce semicolon
* [```98609ac1```](https://github.com/ktmeaton/NCBImeta/commit/98609ac1) split up branch creation and check out
* [```c9ac4a55```](https://github.com/ktmeaton/NCBImeta/commit/c9ac4a55) no tag enforcement yet
* [```03b675ee```](https://github.com/ktmeaton/NCBImeta/commit/03b675ee) forgotten symbols
* [```5d742c2e```](https://github.com/ktmeaton/NCBImeta/commit/5d742c2e) sed replacement check
* [```d6268fe0```](https://github.com/ktmeaton/NCBImeta/commit/d6268fe0) src sha256 refine
* [```4c11045c```](https://github.com/ktmeaton/NCBImeta/commit/4c11045c) path and url fix
* [```d8fc06b9```](https://github.com/ktmeaton/NCBImeta/commit/d8fc06b9) conda var first pass
* [```5f6d8ca5```](https://github.com/ktmeaton/NCBImeta/commit/5f6d8ca5) hide user ref [skip ci]
* [```e48c1e67```](https://github.com/ktmeaton/NCBImeta/commit/e48c1e67) remove dir checks
* [```87266f3c```](https://github.com/ktmeaton/NCBImeta/commit/87266f3c) remove comments
* [```4e1606b1```](https://github.com/ktmeaton/NCBImeta/commit/4e1606b1) travis-ci test post local verify
* [```4554da93```](https://github.com/ktmeaton/NCBImeta/commit/4554da93) fix incorrect slug usage
* [```a0fa57ee```](https://github.com/ktmeaton/NCBImeta/commit/a0fa57ee) show git config
* [```664a2aa6```](https://github.com/ktmeaton/NCBImeta/commit/664a2aa6) directory checking
* [```3d711f34```](https://github.com/ktmeaton/NCBImeta/commit/3d711f34) config param fix
* [```a1a9fa07```](https://github.com/ktmeaton/NCBImeta/commit/a1a9fa07) bioconda update attempt 2
* [```dfbd7e38```](https://github.com/ktmeaton/NCBImeta/commit/dfbd7e38) bioconda update attempt
* [```8f0c064e```](https://github.com/ktmeaton/NCBImeta/commit/8f0c064e) more pwd checks
* [```95801e52```](https://github.com/ktmeaton/NCBImeta/commit/95801e52) pwd test
* [```b6b0856d```](https://github.com/ktmeaton/NCBImeta/commit/b6b0856d) condition testing
* [```cb497125```](https://github.com/ktmeaton/NCBImeta/commit/cb497125) Document the bioconda repo master branch switch

## v0.6.2

### Pull Requests

* [```pull/5```](https://github.com/ktmeaton/NCBImeta/pull/5) Adding Bioconda as an installation option

### Commits

* [```3f8fac57```](https://github.com/ktmeaton/NCBImeta/commit/3f8fac57) PR and Issue linking
* [```4c942075```](https://github.com/ktmeaton/NCBImeta/commit/4c942075) branch match ver release name [skip ci]
* [```d421128d```](https://github.com/ktmeaton/NCBImeta/commit/d421128d) date update v0.6.2 [skip ci]
* [```7dc34dd5```](https://github.com/ktmeaton/NCBImeta/commit/7dc34dd5) automate pypi deploy [skip ci]
* [```e636c6a0```](https://github.com/ktmeaton/NCBImeta/commit/e636c6a0) CI branch on master
* [```37aedb52```](https://github.com/ktmeaton/NCBImeta/commit/37aedb52) Merge branch 'dev'
* [```71ab0e87```](https://github.com/ktmeaton/NCBImeta/commit/71ab0e87) merge dev into master
* [```3de43a83```](https://github.com/ktmeaton/NCBImeta/commit/3de43a83) Merge branch 'dev' of https://github.com/ktmeaton/NCBImeta into dev
* [```7a2a3f3a```](https://github.com/ktmeaton/NCBImeta/commit/7a2a3f3a) version update
* [```a9605de3```](https://github.com/ktmeaton/NCBImeta/commit/a9605de3) version update [skip ci]
* [```2dc59202```](https://github.com/ktmeaton/NCBImeta/commit/2dc59202) ver update to 0.6.2 [skip ci]
* [```92d810d1```](https://github.com/ktmeaton/NCBImeta/commit/92d810d1) Ouput dir is created instead of raising error [skip ci]
* [```bd10fbd7```](https://github.com/ktmeaton/NCBImeta/commit/bd10fbd7) remove slash [skip ci]
* [```f8200686```](https://github.com/ktmeaton/NCBImeta/commit/f8200686) bioconda installation and media as raw links [skip ci]
* [```f994d973```](https://github.com/ktmeaton/NCBImeta/commit/f994d973) Merge branch 'dev' of https://github.com/ktmeaton/NCBImeta into dev
* [```f89ccc53```](https://github.com/ktmeaton/NCBImeta/commit/f89ccc53) CI on dev
* [```5be403b9```](https://github.com/ktmeaton/NCBImeta/commit/5be403b9) v0.6.2 bioconda and output dir
* [```619e1dba```](https://github.com/ktmeaton/NCBImeta/commit/619e1dba) v0.6.2 bioconda and output dir
* [```c6409285```](https://github.com/ktmeaton/NCBImeta/commit/c6409285) Remove the output dir error class [skip ci]
* [```a68d17ad```](https://github.com/ktmeaton/NCBImeta/commit/a68d17ad) Make output dir instead of raising error
* [```a99a4bf4```](https://github.com/ktmeaton/NCBImeta/commit/a99a4bf4) Test if output dir is created [skip ci]
* [```c0c1019d```](https://github.com/ktmeaton/NCBImeta/commit/c0c1019d) Merge pull request #5 from druvus/patch-1
* [```c3fc9a43```](https://github.com/ktmeaton/NCBImeta/commit/c3fc9a43) Paper remove unnecessary reference ncbi website [skip ci]
* [```55e8bf7d```](https://github.com/ktmeaton/NCBImeta/commit/55e8bf7d) Adding Bioconda as an installation option
* [```17d6ba13```](https://github.com/ktmeaton/NCBImeta/commit/17d6ba13) Zenodo citation and badge [skip ci]

## v0.6.1

### Commits

* [```0ab131a5```](https://github.com/ktmeaton/NCBImeta/commit/0ab131a5) update schema doc for new xml parsing [skip ci]
* [```680ddb26```](https://github.com/ktmeaton/NCBImeta/commit/680ddb26) ref and capitalization [skip ci]
* [```edaa3b00```](https://github.com/ktmeaton/NCBImeta/commit/edaa3b00) python 3 spacing [skip ci]
* [```49132f98```](https://github.com/ktmeaton/NCBImeta/commit/49132f98) re-enable master CI
* [```fb3a01cf```](https://github.com/ktmeaton/NCBImeta/commit/fb3a01cf) aeruginosa text db
* [```7dd98044```](https://github.com/ktmeaton/NCBImeta/commit/7dd98044) aeruginosa paper db
* [```c50677f2```](https://github.com/ktmeaton/NCBImeta/commit/c50677f2) add release names
* [```a93d1a11```](https://github.com/ktmeaton/NCBImeta/commit/a93d1a11) Name for upcoming v0.6.1
* [```a24a82f6```](https://github.com/ktmeaton/NCBImeta/commit/a24a82f6) release to development
* [```32d162ba```](https://github.com/ktmeaton/NCBImeta/commit/32d162ba) prematurely add v0.6.1 links
* [```88d6d0c1```](https://github.com/ktmeaton/NCBImeta/commit/88d6d0c1) format and typos
* [```5ccd379a```](https://github.com/ktmeaton/NCBImeta/commit/5ccd379a) conceptual rearrange
* [```be940044```](https://github.com/ktmeaton/NCBImeta/commit/be940044) Initial commit desc
* [```37068fc7```](https://github.com/ktmeaton/NCBImeta/commit/37068fc7) complicated compare
* [```8637310f```](https://github.com/ktmeaton/NCBImeta/commit/8637310f) link fix
* [```900c6a55```](https://github.com/ktmeaton/NCBImeta/commit/900c6a55) test commit compare
* [```8efd0815```](https://github.com/ktmeaton/NCBImeta/commit/8efd0815) hyperlink exp
* [```3b17a514```](https://github.com/ktmeaton/NCBImeta/commit/3b17a514) aeruginosa config file
* [```fdc4672b```](https://github.com/ktmeaton/NCBImeta/commit/fdc4672b) list format tables
* [```4bbeb109```](https://github.com/ktmeaton/NCBImeta/commit/4bbeb109) rename
* [```f7f53111```](https://github.com/ktmeaton/NCBImeta/commit/f7f53111) remove old gif
* [```27d1a2de```](https://github.com/ktmeaton/NCBImeta/commit/27d1a2de) biosample gif
* [```e5716aa2```](https://github.com/ktmeaton/NCBImeta/commit/e5716aa2) gif larger
* [```79b82f1e```](https://github.com/ktmeaton/NCBImeta/commit/79b82f1e) gif rename
* [```867775af```](https://github.com/ktmeaton/NCBImeta/commit/867775af) image swap
* [```37dd9c90```](https://github.com/ktmeaton/NCBImeta/commit/37dd9c90) gif full path, reformat print
* [```107fa144```](https://github.com/ktmeaton/NCBImeta/commit/107fa144) State CLI Application
* [```69f17655```](https://github.com/ktmeaton/NCBImeta/commit/69f17655) specify CLI application
* [```8a72fc83```](https://github.com/ktmeaton/NCBImeta/commit/8a72fc83) asciicast
* [```7d0a4ccc```](https://github.com/ktmeaton/NCBImeta/commit/7d0a4ccc) dependency reference to file
* [```e2674d74```](https://github.com/ktmeaton/NCBImeta/commit/e2674d74) JOSS file additions
* [```bec9abe9```](https://github.com/ktmeaton/NCBImeta/commit/bec9abe9) update record number
* [```08b5e589```](https://github.com/ktmeaton/NCBImeta/commit/08b5e589) uncomment xml print
* [```084f2006```](https://github.com/ktmeaton/NCBImeta/commit/084f2006) jpg rename
* [```d1bd7ae9```](https://github.com/ktmeaton/NCBImeta/commit/d1bd7ae9) doc xml printing, disable master ci
* [```bd92205b```](https://github.com/ktmeaton/NCBImeta/commit/bd92205b) doc update API optional [skip ci]
* [```2c72df89```](https://github.com/ktmeaton/NCBImeta/commit/2c72df89) doc update node parse [skip ci]
* [```c7a77e46```](https://github.com/ktmeaton/NCBImeta/commit/c7a77e46) new prettyprint func [skip ci]
* [```c78b8c25```](https://github.com/ktmeaton/NCBImeta/commit/c78b8c25) update version [skip ci]
* [```ec6ebf8f```](https://github.com/ktmeaton/NCBImeta/commit/ec6ebf8f) edit description [skip ci]
* [```daa228b4```](https://github.com/ktmeaton/NCBImeta/commit/daa228b4) JOSS Paper
* [```82b7b4f3```](https://github.com/ktmeaton/NCBImeta/commit/82b7b4f3) Reformat headings [skip ci]
* [```7cbabc9a```](https://github.com/ktmeaton/NCBImeta/commit/7cbabc9a) config file name update [skip ci]
* [```9808a6bd```](https://github.com/ktmeaton/NCBImeta/commit/9808a6bd) Clarify and formatting [skip ci]

## v0.6.0

### Commits

* [```44805c70```](https://github.com/ktmeaton/NCBImeta/commit/44805c70) remove exp code
* [```717431a8```](https://github.com/ktmeaton/NCBImeta/commit/717431a8) license badge update [skip ci]
* [```05b3aeaa```](https://github.com/ktmeaton/NCBImeta/commit/05b3aeaa) Merge branch 'dev'
* [```cb575aea```](https://github.com/ktmeaton/NCBImeta/commit/cb575aea) conflict resolve
* [```b489ca62```](https://github.com/ktmeaton/NCBImeta/commit/b489ca62) v0.6.0 updates [skip ci]
* [```778ab45f```](https://github.com/ktmeaton/NCBImeta/commit/778ab45f) ignore test database and log files [skip ci]
* [```b3a310c5```](https://github.com/ktmeaton/NCBImeta/commit/b3a310c5) Single quote sql query fix
* [```0e0be9cb```](https://github.com/ktmeaton/NCBImeta/commit/0e0be9cb) XML overhaul test update
* [```b0dc3bcc```](https://github.com/ktmeaton/NCBImeta/commit/b0dc3bcc) successful run
* [```848cb532```](https://github.com/ktmeaton/NCBImeta/commit/848cb532) quotation fix
* [```e03597af```](https://github.com/ktmeaton/NCBImeta/commit/e03597af) bugfixes for xml_search
* [```e0ff551a```](https://github.com/ktmeaton/NCBImeta/commit/e0ff551a) remove old search and flatten functions
* [```546c368e```](https://github.com/ktmeaton/NCBImeta/commit/546c368e) Before function move
* [```9d97bbe0```](https://github.com/ktmeaton/NCBImeta/commit/9d97bbe0) efetch part functional for sra
* [```ab88a2c4```](https://github.com/ktmeaton/NCBImeta/commit/ab88a2c4) before nucleotide switch to efetch
* [```393466d3```](https://github.com/ktmeaton/NCBImeta/commit/393466d3) bioproject with efetch
* [```be00272b```](https://github.com/ktmeaton/NCBImeta/commit/be00272b) efetch switch
* [```ee6cf173```](https://github.com/ktmeaton/NCBImeta/commit/ee6cf173) lxml overhaul
* [```31de27e9```](https://github.com/ktmeaton/NCBImeta/commit/31de27e9) before efetch switch
* [```1b720a82```](https://github.com/ktmeaton/NCBImeta/commit/1b720a82) xml in xml parsing
* [```d33d58e5```](https://github.com/ktmeaton/NCBImeta/commit/d33d58e5) cdata testing
* [```940ba360```](https://github.com/ktmeaton/NCBImeta/commit/940ba360) lxml experiment
* [```88b7c0de```](https://github.com/ktmeaton/NCBImeta/commit/88b7c0de) Config and Schema README as full path [skip ci]
* [```0b1030dd```](https://github.com/ktmeaton/NCBImeta/commit/0b1030dd) ignore build dir [skip ci]
* [```e4be83d4```](https://github.com/ktmeaton/NCBImeta/commit/e4be83d4) change codecov to branch master [skip ci]
* [```f3b30de3```](https://github.com/ktmeaton/NCBImeta/commit/f3b30de3) v0.5.0 version update [skip ci]

## v0.5.0

### Commits

* [```046f68be```](https://github.com/ktmeaton/NCBImeta/commit/046f68be) travis-ci only master
* [```cecb1aae```](https://github.com/ktmeaton/NCBImeta/commit/cecb1aae) Merge branch 'dev'
* [```0d844574```](https://github.com/ktmeaton/NCBImeta/commit/0d844574) v0.5.0 updates [skip ci]
* [```849f7335```](https://github.com/ktmeaton/NCBImeta/commit/849f7335) Explain default slow download [skip ci]
* [```b810c2ad```](https://github.com/ktmeaton/NCBImeta/commit/b810c2ad) Fix codecov badge link [skip ci]
* [```5d4c83bc```](https://github.com/ktmeaton/NCBImeta/commit/5d4c83bc) Switch release badge to PyPI [skip ci]
* [```4112027b```](https://github.com/ktmeaton/NCBImeta/commit/4112027b) codecov url fix [skip ci]
* [```29a80f52```](https://github.com/ktmeaton/NCBImeta/commit/29a80f52) try codecov badge with dev [skip ci]
* [```1488da73```](https://github.com/ktmeaton/NCBImeta/commit/1488da73) Clear debugging output
* [```199559b5```](https://github.com/ktmeaton/NCBImeta/commit/199559b5) formatting [skip ci]
* [```e6bd3b70```](https://github.com/ktmeaton/NCBImeta/commit/e6bd3b70) document column_index bugfix [skip ci]
* [```f560dfe5```](https://github.com/ktmeaton/NCBImeta/commit/f560dfe5) Correct bioproject accession
* [```10e1d2ae```](https://github.com/ktmeaton/NCBImeta/commit/10e1d2ae) column_index reposition
* [```75edf10e```](https://github.com/ktmeaton/NCBImeta/commit/75edf10e) http error catching for esearch
* [```1bc53175```](https://github.com/ktmeaton/NCBImeta/commit/1bc53175) strpath for tmpdir [skip ci]
* [```fd509f42```](https://github.com/ktmeaton/NCBImeta/commit/fd509f42) Fix sys path and list comparison
* [```3eb83673```](https://github.com/ktmeaton/NCBImeta/commit/3eb83673) Remove Python 3.4 linux build [skip ci]
* [```6bd4d8d8```](https://github.com/ktmeaton/NCBImeta/commit/6bd4d8d8) require time module
* [```8df42459```](https://github.com/ktmeaton/NCBImeta/commit/8df42459) test db path fix [skip ci]
* [```1488745a```](https://github.com/ktmeaton/NCBImeta/commit/1488745a) Python 3.5+ required (no more 3.4) [skip ci]
* [```d5624e19```](https://github.com/ktmeaton/NCBImeta/commit/d5624e19) travis troubleshooting
* [```763c54f8```](https://github.com/ktmeaton/NCBImeta/commit/763c54f8) bash uploader for codecov
* [```b807bb02```](https://github.com/ktmeaton/NCBImeta/commit/b807bb02) example config reset
* [```595959b9```](https://github.com/ktmeaton/NCBImeta/commit/595959b9) Merge branch 'dev' of https://github.com/ktmeaton/NCBImeta into dev
* [```4c9b15bc```](https://github.com/ktmeaton/NCBImeta/commit/4c9b15bc) full travis-ci test build
* [```8f5f38b2```](https://github.com/ktmeaton/NCBImeta/commit/8f5f38b2) full travis-ci test build
* [```e2ede2c0```](https://github.com/ktmeaton/NCBImeta/commit/e2ede2c0) pytest for non-flat mode
* [```7d4ec415```](https://github.com/ktmeaton/NCBImeta/commit/7d4ec415) prototype master table check
* [```ea642b8c```](https://github.com/ktmeaton/NCBImeta/commit/ea642b8c) Move the HTTPErrorCatch method to utilities [skip ci]
* [```dd55698e```](https://github.com/ktmeaton/NCBImeta/commit/dd55698e) pytest db verify SRA
* [```8de4febd```](https://github.com/ktmeaton/NCBImeta/commit/8de4febd) pytesting and db verification [skip ci]
* [```bc7f6b29```](https://github.com/ktmeaton/NCBImeta/commit/bc7f6b29) Update test biosample metadata
* [```6771d87a```](https://github.com/ktmeaton/NCBImeta/commit/6771d87a) Exclude temporary test files [skip ci]
* [```b2df525d```](https://github.com/ktmeaton/NCBImeta/commit/b2df525d) pytest pubmed [skip ci]
* [```a88ba92d```](https://github.com/ktmeaton/NCBImeta/commit/a88ba92d) pytest export nucleotide table values
* [```a21df9eb```](https://github.com/ktmeaton/NCBImeta/commit/a21df9eb) pytest conftest assembly and bioproject
* [```fe28e8b3```](https://github.com/ktmeaton/NCBImeta/commit/fe28e8b3) pytest export assembly and bioproject
* [```918f1e3c```](https://github.com/ktmeaton/NCBImeta/commit/918f1e3c) remove annot file [skip ci]
* [```26b8d15f```](https://github.com/ktmeaton/NCBImeta/commit/26b8d15f) pytest join
* [```87ca364a```](https://github.com/ktmeaton/NCBImeta/commit/87ca364a) Comment clarification [skip ci]
* [```e84af09b```](https://github.com/ktmeaton/NCBImeta/commit/e84af09b) pytest annotatereplace
* [```90da7506```](https://github.com/ktmeaton/NCBImeta/commit/90da7506) annotateconcat 90% cov
* [```34725c9e```](https://github.com/ktmeaton/NCBImeta/commit/34725c9e) fix indentation
* [```4eea6b32```](https://github.com/ktmeaton/NCBImeta/commit/4eea6b32) annotation file for pytest [skip ci]
* [```1a4f0682```](https://github.com/ktmeaton/NCBImeta/commit/1a4f0682) specify pytest files in order
* [```98a4221f```](https://github.com/ktmeaton/NCBImeta/commit/98a4221f) module import in test dir [skip ci]
* [```88279923```](https://github.com/ktmeaton/NCBImeta/commit/88279923) sra table metadata re-fix [skip ci]
* [```debd22ab```](https://github.com/ktmeaton/NCBImeta/commit/debd22ab) debug help points [skip ci]
* [```7150dc32```](https://github.com/ktmeaton/NCBImeta/commit/7150dc32) generic test db name [skip ci]
* [```cac8a129```](https://github.com/ktmeaton/NCBImeta/commit/cac8a129) pytest annotateconcatenate [skip ci]
* [```e16d6e8c```](https://github.com/ktmeaton/NCBImeta/commit/e16d6e8c) formatting [skip ci]
* [```469bb6a6```](https://github.com/ktmeaton/NCBImeta/commit/469bb6a6) proper support for SRA [skip ci]
* [```b5e57539```](https://github.com/ktmeaton/NCBImeta/commit/b5e57539) codecov badge [skip ci]
* [```26221240```](https://github.com/ktmeaton/NCBImeta/commit/26221240) remove extra script
* [```081a76e3```](https://github.com/ktmeaton/NCBImeta/commit/081a76e3) execute permissions [skip ci]
* [```ef5b084c```](https://github.com/ktmeaton/NCBImeta/commit/ef5b084c) remove the test1 ref
* [```14133ba0```](https://github.com/ktmeaton/NCBImeta/commit/14133ba0) update execute permissions again
* [```90d9bfbc```](https://github.com/ktmeaton/NCBImeta/commit/90d9bfbc) test for output dir
* [```9a1533e3```](https://github.com/ktmeaton/NCBImeta/commit/9a1533e3) config data yaml test [skip ci]
* [```12b67361```](https://github.com/ktmeaton/NCBImeta/commit/12b67361) old cov [skip ci]
* [```14ce1bf7```](https://github.com/ktmeaton/NCBImeta/commit/14ce1bf7) remove extraneous [skip ci]
* [```7d8037b6```](https://github.com/ktmeaton/NCBImeta/commit/7d8037b6) simplify name [skip ci]
* [```f7b4502e```](https://github.com/ktmeaton/NCBImeta/commit/f7b4502e) yaml file error test [skip ci]
* [```727e547b```](https://github.com/ktmeaton/NCBImeta/commit/727e547b) main test [skip ci]
* [```1efd2793```](https://github.com/ktmeaton/NCBImeta/commit/1efd2793) give unique names [skip ci]
* [```dac651b8```](https://github.com/ktmeaton/NCBImeta/commit/dac651b8) print config file path [skip ci]
* [```e2de1968```](https://github.com/ktmeaton/NCBImeta/commit/e2de1968) yaml file test [skip ci]
* [```0eb16810```](https://github.com/ktmeaton/NCBImeta/commit/0eb16810) yaml files for testing [skip ci]
* [```805efa69```](https://github.com/ktmeaton/NCBImeta/commit/805efa69) super simple config [skip ci]
* [```de49481a```](https://github.com/ktmeaton/NCBImeta/commit/de49481a) Proper os join
* [```7a9be44c```](https://github.com/ktmeaton/NCBImeta/commit/7a9be44c) pytest Errors complete
* [```8afd1ce9```](https://github.com/ktmeaton/NCBImeta/commit/8afd1ce9) pytest errors
* [```ebfc8af1```](https://github.com/ktmeaton/NCBImeta/commit/ebfc8af1) duplicate column debugging [skip ci]
* [```b663a974```](https://github.com/ktmeaton/NCBImeta/commit/b663a974) str repr of errors proper return [skip ci]
* [```c51cd55b```](https://github.com/ktmeaton/NCBImeta/commit/c51cd55b) check for duplicate column names
* [```1725bc6b```](https://github.com/ktmeaton/NCBImeta/commit/1725bc6b) pytest and coverage [skip ci]
* [```35ffd83b```](https://github.com/ktmeaton/NCBImeta/commit/35ffd83b) cleanup docstring [skip ci]
* [```640ed49e```](https://github.com/ktmeaton/NCBImeta/commit/640ed49e) use tmp dir and files
* [```13b0edb6```](https://github.com/ktmeaton/NCBImeta/commit/13b0edb6) remove unnecessary code [skip ci]
* [```fb3bff1d```](https://github.com/ktmeaton/NCBImeta/commit/fb3bff1d) Re-enable full travis-ci [skip ci]
* [```ac82f6b3```](https://github.com/ktmeaton/NCBImeta/commit/ac82f6b3) remove unnecessary type checking
* [```1a801bc2```](https://github.com/ktmeaton/NCBImeta/commit/1a801bc2) codecov take 2
* [```e6b8977f```](https://github.com/ktmeaton/NCBImeta/commit/e6b8977f) codecov coverage test upload [skip ci]
* [```40a601f1```](https://github.com/ktmeaton/NCBImeta/commit/40a601f1) remove testing modules to travis requirements
* [```60213802```](https://github.com/ktmeaton/NCBImeta/commit/60213802) codecov require and test
* [```a35e2513```](https://github.com/ktmeaton/NCBImeta/commit/a35e2513) pytest-cov requirement
* [```ef914477```](https://github.com/ktmeaton/NCBImeta/commit/ef914477) pytest troubleshooting
* [```34850960```](https://github.com/ktmeaton/NCBImeta/commit/34850960) --cov-report= proper parameter
* [```2bbd4258```](https://github.com/ktmeaton/NCBImeta/commit/2bbd4258) pytest and codecov
* [```ba90d3f2```](https://github.com/ktmeaton/NCBImeta/commit/ba90d3f2) Longer description [skip ci]
* [```971a6335```](https://github.com/ktmeaton/NCBImeta/commit/971a6335) Comment out debugging [skip ci]
* [```eb75e654```](https://github.com/ktmeaton/NCBImeta/commit/eb75e654) v0.4.3 changelog updates [skip ci]
* [```3fd7cec2```](https://github.com/ktmeaton/NCBImeta/commit/3fd7cec2) Remove unicode ref
* [```e9086255```](https://github.com/ktmeaton/NCBImeta/commit/e9086255) Remove XmlXXXConfig Functions
* [```120cfa5e```](https://github.com/ktmeaton/NCBImeta/commit/120cfa5e) Extended description, remove troubleshooting printout [skip ci]
* [```a462b23d```](https://github.com/ktmeaton/NCBImeta/commit/a462b23d) Namespace clarify
* [```6bde5f20```](https://github.com/ktmeaton/NCBImeta/commit/6bde5f20) Expanded type checking
* [```1c7bd73a```](https://github.com/ktmeaton/NCBImeta/commit/1c7bd73a) StringElement items troubleshooting
* [```62aceba4```](https://github.com/ktmeaton/NCBImeta/commit/62aceba4) Extra type checking
* [```30423509```](https://github.com/ktmeaton/NCBImeta/commit/30423509) Fix internal recursion function name
* [```5b4de6b2```](https://github.com/ktmeaton/NCBImeta/commit/5b4de6b2) Correct and document flatten_dict implementation
* [```4bdedf32```](https://github.com/ktmeaton/NCBImeta/commit/4bdedf32) More testing
* [```efb3b86f```](https://github.com/ktmeaton/NCBImeta/commit/efb3b86f) Utilities testing
* [```269ca821```](https://github.com/ktmeaton/NCBImeta/commit/269ca821) Recode unicode and remove sys
* [```a20daf6a```](https://github.com/ktmeaton/NCBImeta/commit/a20daf6a) Description in header [skip ci]
* [```05f693ec```](https://github.com/ktmeaton/NCBImeta/commit/05f693ec) Remove sys module and flushprint [skip ci]
* [```e9efc279```](https://github.com/ktmeaton/NCBImeta/commit/e9efc279) cleanup [skip ci]
* [```d92dc091```](https://github.com/ktmeaton/NCBImeta/commit/d92dc091) NCBImetaErrors namespace
* [```fc983047```](https://github.com/ktmeaton/NCBImeta/commit/fc983047) Remove io module
* [```78032d93```](https://github.com/ktmeaton/NCBImeta/commit/78032d93) Remove sys module [skip ci]
* [```f76c6dc1```](https://github.com/ktmeaton/NCBImeta/commit/f76c6dc1) Error classes documented
* [```45ddc1bb```](https://github.com/ktmeaton/NCBImeta/commit/45ddc1bb) Description in header [skip ci]
* [```1f86ec3d```](https://github.com/ktmeaton/NCBImeta/commit/1f86ec3d) Remove sys module and flushprint
* [```b1b668a0```](https://github.com/ktmeaton/NCBImeta/commit/b1b668a0) Remove sys module
* [```8400f5f2```](https://github.com/ktmeaton/NCBImeta/commit/8400f5f2) Code doc first pass [skip ci]
* [```b681b08d```](https://github.com/ktmeaton/NCBImeta/commit/b681b08d) Reduce code redundancy
* [```86be6028```](https://github.com/ktmeaton/NCBImeta/commit/86be6028) Document UpdateDB [skip CI]
* [```d7a22f2e```](https://github.com/ktmeaton/NCBImeta/commit/d7a22f2e) Force unicode str recode
* [```691111ab```](https://github.com/ktmeaton/NCBImeta/commit/691111ab) Cleanup, document HTTPErrorCatch
* [```778817bf```](https://github.com/ktmeaton/NCBImeta/commit/778817bf) import slim down
* [```470e3422```](https://github.com/ktmeaton/NCBImeta/commit/470e3422) flush print true case sensitive
* [```4addaaa3```](https://github.com/ktmeaton/NCBImeta/commit/4addaaa3) Travis CI for dev
* [```5ae1c7af```](https://github.com/ktmeaton/NCBImeta/commit/5ae1c7af) Replace flushprint method python3
* [```90bedfba```](https://github.com/ktmeaton/NCBImeta/commit/90bedfba) Merge branch 'master' into dev
* [```5f3548e1```](https://github.com/ktmeaton/NCBImeta/commit/5f3548e1) Header change
* [```12722915```](https://github.com/ktmeaton/NCBImeta/commit/12722915) Gif reupload and typo fix
* [```359b0ee6```](https://github.com/ktmeaton/NCBImeta/commit/359b0ee6) Nucleotidet typo fix
* [```5e56fee8```](https://github.com/ktmeaton/NCBImeta/commit/5e56fee8) gif reupload
* [```6f6d74c5```](https://github.com/ktmeaton/NCBImeta/commit/6f6d74c5) Delete broken DB gif
* [```4e53ed01```](https://github.com/ktmeaton/NCBImeta/commit/4e53ed01) Delete broken CLI gif
* [```546edd68```](https://github.com/ktmeaton/NCBImeta/commit/546edd68) Move requirements higher [skip ci]

## v0.4.2

### Commits

* [```37adf6ec```](https://github.com/ktmeaton/NCBImeta/commit/37adf6ec) v0.4.2 finish
* [```296f34fa```](https://github.com/ktmeaton/NCBImeta/commit/296f34fa) Require cloning github [skip ci]
* [```2847a272```](https://github.com/ktmeaton/NCBImeta/commit/2847a272) section break [skip ci]
* [```b1cb2ce3```](https://github.com/ktmeaton/NCBImeta/commit/b1cb2ce3) Remove test warning [skip ci]
* [```994f4b1e```](https://github.com/ktmeaton/NCBImeta/commit/994f4b1e) Web portal links, src dir fix [skip ci]
* [```359f572f```](https://github.com/ktmeaton/NCBImeta/commit/359f572f) Extra quotes remove [skip ci]
* [```8edd0a2a```](https://github.com/ktmeaton/NCBImeta/commit/8edd0a2a) Fix module import and author credit
* [```01f45929```](https://github.com/ktmeaton/NCBImeta/commit/01f45929) fix quote
* [```b477b425```](https://github.com/ktmeaton/NCBImeta/commit/b477b425) missing comma
* [```02f81dcf```](https://github.com/ktmeaton/NCBImeta/commit/02f81dcf) pip install match PyPI Repo [skip ci]
* [```2c856e20```](https://github.com/ktmeaton/NCBImeta/commit/2c856e20) Remove underscores
* [```c6cd147b```](https://github.com/ktmeaton/NCBImeta/commit/c6cd147b) Remove src dir links [skip ci]
* [```6c5ba1e4```](https://github.com/ktmeaton/NCBImeta/commit/6c5ba1e4) Install update [skip ci]
* [```473ca4c1```](https://github.com/ktmeaton/NCBImeta/commit/473ca4c1) Home url and markdown desc [skip ci]
* [```b8b71bd9```](https://github.com/ktmeaton/NCBImeta/commit/b8b71bd9) Short description [skip CI]
* [```82a05e65```](https://github.com/ktmeaton/NCBImeta/commit/82a05e65) pip install included
* [```05cb9a0e```](https://github.com/ktmeaton/NCBImeta/commit/05cb9a0e) branch restrict rearrange
* [```f3510328```](https://github.com/ktmeaton/NCBImeta/commit/f3510328) travis troubleshooting part 2
* [```d9f58431```](https://github.com/ktmeaton/NCBImeta/commit/d9f58431) travis troubleshooting
* [```f574720d```](https://github.com/ktmeaton/NCBImeta/commit/f574720d) requirements note
* [```a411aea1```](https://github.com/ktmeaton/NCBImeta/commit/a411aea1) v0.4.2 starts, PyPI
* [```bc572124```](https://github.com/ktmeaton/NCBImeta/commit/bc572124) Merge branch 'dev'
* [```9909a130```](https://github.com/ktmeaton/NCBImeta/commit/9909a130) PyPI Preparation
* [```27568510```](https://github.com/ktmeaton/NCBImeta/commit/27568510) setup.py install works!
* [```b549268d```](https://github.com/ktmeaton/NCBImeta/commit/b549268d) Rename src dir to ncbimeta
* [```d25112a9```](https://github.com/ktmeaton/NCBImeta/commit/d25112a9) More setup.py config
* [```4236ceca```](https://github.com/ktmeaton/NCBImeta/commit/4236ceca) requirements.txt control for setup
* [```0978441e```](https://github.com/ktmeaton/NCBImeta/commit/0978441e) Added setup.py
* [```ff158f0e```](https://github.com/ktmeaton/NCBImeta/commit/ff158f0e) Remove underscores in program
* [```f074685b```](https://github.com/ktmeaton/NCBImeta/commit/f074685b) Merge branch 'master' into dev
* [```0f33450f```](https://github.com/ktmeaton/NCBImeta/commit/0f33450f) Run Travis CI only on master
* [```a7d90383```](https://github.com/ktmeaton/NCBImeta/commit/a7d90383) v0.4.1 version update

## v0.4.1

### Commits

* [```69941ada```](https://github.com/ktmeaton/NCBImeta/commit/69941ada) v0.4.1 release
* [```abf1e032```](https://github.com/ktmeaton/NCBImeta/commit/abf1e032) Rearrange and simplify [skip ci]
* [```92c7b601```](https://github.com/ktmeaton/NCBImeta/commit/92c7b601) Update CLI gif to v0.4.0 [skip ci]
* [```86b5eb58```](https://github.com/ktmeaton/NCBImeta/commit/86b5eb58) relative path change gif [skip ci]
* [```1479a5ac```](https://github.com/ktmeaton/NCBImeta/commit/1479a5ac) Merge branch 'dev'
* [```f8112bd1```](https://github.com/ktmeaton/NCBImeta/commit/f8112bd1) url update [skip ci]
* [```770f7e5b```](https://github.com/ktmeaton/NCBImeta/commit/770f7e5b) gh pages url [skip ci]
* [```00c72f40```](https://github.com/ktmeaton/NCBImeta/commit/00c72f40) Merge branch 'master' into dev
* [```46931e02```](https://github.com/ktmeaton/NCBImeta/commit/46931e02) Dependency clarification
* [```83ceb7fb```](https://github.com/ktmeaton/NCBImeta/commit/83ceb7fb) Merge branch 'dev'
* [```ac09f722```](https://github.com/ktmeaton/NCBImeta/commit/ac09f722) Multi-python linux
* [```6731c998```](https://github.com/ktmeaton/NCBImeta/commit/6731c998) Remove accessory program section [skip ci]
* [```50121e98```](https://github.com/ktmeaton/NCBImeta/commit/50121e98) Spacing [skip ci]
* [```409dfcc9```](https://github.com/ktmeaton/NCBImeta/commit/409dfcc9) Describe OS support [skip ci]
* [```6defc7b7```](https://github.com/ktmeaton/NCBImeta/commit/6defc7b7) Description simplify [skip ci]
* [```790dd776```](https://github.com/ktmeaton/NCBImeta/commit/790dd776) Links and issues doc
* [```9c2de04e```](https://github.com/ktmeaton/NCBImeta/commit/9c2de04e) Update program description
* [```da2c18c3```](https://github.com/ktmeaton/NCBImeta/commit/da2c18c3) Simple Linux and MacOS Build
* [```04151edb```](https://github.com/ktmeaton/NCBImeta/commit/04151edb) Remove URL Error code print
* [```c67c9d5f```](https://github.com/ktmeaton/NCBImeta/commit/c67c9d5f) Travis-CI MacOS Supported!
* [```f22cedcb```](https://github.com/ktmeaton/NCBImeta/commit/f22cedcb) esearch cmd rearrange
* [```c3b34440```](https://github.com/ktmeaton/NCBImeta/commit/c3b34440) kwargs typo
* [```df09ab4b```](https://github.com/ktmeaton/NCBImeta/commit/df09ab4b) Typo http_metthod
* [```4e38e926```](https://github.com/ktmeaton/NCBImeta/commit/4e38e926) Try to catch URL Error
* [```e8b43d98```](https://github.com/ktmeaton/NCBImeta/commit/e8b43d98) Retry macos testing
* [```fb9c1e77```](https://github.com/ktmeaton/NCBImeta/commit/fb9c1e77) Macos test first
* [```20762ef8```](https://github.com/ktmeaton/NCBImeta/commit/20762ef8) Pip3 change and multi-python test
* [```19eb7e98```](https://github.com/ktmeaton/NCBImeta/commit/19eb7e98) Bugfix plus exclude windwos from travis
* [```8571960e```](https://github.com/ktmeaton/NCBImeta/commit/8571960e) Test with python3 prefix
* [```d4a1b519```](https://github.com/ktmeaton/NCBImeta/commit/d4a1b519) Windows travis-ci test
* [```5e8e305c```](https://github.com/ktmeaton/NCBImeta/commit/5e8e305c) Py 3.7.4 multi-dist test
* [```65be703b```](https://github.com/ktmeaton/NCBImeta/commit/65be703b) Add Travis-CI Badge [skip ci]
* [```0f8c65f3```](https://github.com/ktmeaton/NCBImeta/commit/0f8c65f3) Just test python 3.6 on dev
* [```be9fc9b0```](https://github.com/ktmeaton/NCBImeta/commit/be9fc9b0) v0.4.1 preparation [skip ci]
* [```53006214```](https://github.com/ktmeaton/NCBImeta/commit/53006214) Annotation script bugfixes [skip ci]
* [```c0218440```](https://github.com/ktmeaton/NCBImeta/commit/c0218440) Updated biopython requirement to 1.74 [skip ci]
* [```535aa7fc```](https://github.com/ktmeaton/NCBImeta/commit/535aa7fc) Add numpy to requirements although packaged with biopython
* [```b837b59f```](https://github.com/ktmeaton/NCBImeta/commit/b837b59f) Only dev branches for travis-ci
* [```56105f96```](https://github.com/ktmeaton/NCBImeta/commit/56105f96) Add execute permissions to scripts
* [```d4f8c295```](https://github.com/ktmeaton/NCBImeta/commit/d4f8c295) Remove --user flag from travis-ci pip
* [```9a436e42```](https://github.com/ktmeaton/NCBImeta/commit/9a436e42) Travis CI Integration
* [```df2ddfd0```](https://github.com/ktmeaton/NCBImeta/commit/df2ddfd0) Merge branch 'master' into dev
* [```3b6276d7```](https://github.com/ktmeaton/NCBImeta/commit/3b6276d7) Simplify dependency written desc

## v0.4.0

### Commits

* [```c7b320df```](https://github.com/ktmeaton/NCBImeta/commit/c7b320df) v0.4.0 release changes
* [```5a4a799c```](https://github.com/ktmeaton/NCBImeta/commit/5a4a799c) Merge branch 'dev'
* [```17527168```](https://github.com/ktmeaton/NCBImeta/commit/17527168) Variable formatting
* [```3a910cb9```](https://github.com/ktmeaton/NCBImeta/commit/3a910cb9) Formatting test
* [```70ee5cb3```](https://github.com/ktmeaton/NCBImeta/commit/70ee5cb3) Improved parameter description
* [```882140af```](https://github.com/ktmeaton/NCBImeta/commit/882140af) annot file updates
* [```3c257cd6```](https://github.com/ktmeaton/NCBImeta/commit/3c257cd6) Fixed program description.
* [```92bee130```](https://github.com/ktmeaton/NCBImeta/commit/92bee130) Improved description
* [```6471de76```](https://github.com/ktmeaton/NCBImeta/commit/6471de76) Additional explanation for join and export
* [```9a0e1957```](https://github.com/ktmeaton/NCBImeta/commit/9a0e1957) Section rearrange
* [```0ccfe8a9```](https://github.com/ktmeaton/NCBImeta/commit/0ccfe8a9) Fix incorrect release version links
* [```8ae3a9c5```](https://github.com/ktmeaton/NCBImeta/commit/8ae3a9c5) more typos
* [```d1df6f8c```](https://github.com/ktmeaton/NCBImeta/commit/d1df6f8c) Typos etc
* [```87802934```](https://github.com/ktmeaton/NCBImeta/commit/87802934) Release links
* [```82eacccb```](https://github.com/ktmeaton/NCBImeta/commit/82eacccb) Main README simplify
* [```28c28d74```](https://github.com/ktmeaton/NCBImeta/commit/28c28d74) Force metadata concatenation even if same value
* [```acb33f82```](https://github.com/ktmeaton/NCBImeta/commit/acb33f82) Simply annotation file name
* [```32de6625```](https://github.com/ktmeaton/NCBImeta/commit/32de6625) Print msg upon successful match
* [```edca4c0d```](https://github.com/ktmeaton/NCBImeta/commit/edca4c0d) Rmv records outside time boundary
* [```f42d0a61```](https://github.com/ktmeaton/NCBImeta/commit/f42d0a61) Rmv annot file outside time boundary
* [```792d683b```](https://github.com/ktmeaton/NCBImeta/commit/792d683b) Output msg for succesful match
* [```047124b1```](https://github.com/ktmeaton/NCBImeta/commit/047124b1) Schema updates
* [```a33cc251```](https://github.com/ktmeaton/NCBImeta/commit/a33cc251) Removed scripts folder with R code
* [```930f2e52```](https://github.com/ktmeaton/NCBImeta/commit/930f2e52) Clarify whitespace importance
* [```bb86c43d```](https://github.com/ktmeaton/NCBImeta/commit/bb86c43d) Remove unsupported schema txt files
* [```5e2a608a```](https://github.com/ktmeaton/NCBImeta/commit/5e2a608a) Doc formatting
* [```7008c663```](https://github.com/ktmeaton/NCBImeta/commit/7008c663) Doc formatting and spacing
* [```a1c5bed5```](https://github.com/ktmeaton/NCBImeta/commit/a1c5bed5) Schema documentation update for yaml
* [```b7edcea2```](https://github.com/ktmeaton/NCBImeta/commit/b7edcea2) schema files in yaml format
* [```4762bb94```](https://github.com/ktmeaton/NCBImeta/commit/4762bb94) Added comment field to Pubmed table
* [```e431a2b0```](https://github.com/ktmeaton/NCBImeta/commit/e431a2b0) badge typo and gif to be updated
* [```52ad0019```](https://github.com/ktmeaton/NCBImeta/commit/52ad0019) Sanity commit
* [```23005247```](https://github.com/ktmeaton/NCBImeta/commit/23005247) update, MWE back to plague
* [```0014b8c5```](https://github.com/ktmeaton/NCBImeta/commit/0014b8c5) MWE organism change and direct execution without python3 call
* [```4ebc22b9```](https://github.com/ktmeaton/NCBImeta/commit/4ebc22b9) Shebang updates for source files
* [```e81f17e0```](https://github.com/ktmeaton/NCBImeta/commit/e81f17e0) Annotation file rename
* [```a883af4f```](https://github.com/ktmeaton/NCBImeta/commit/a883af4f) Switch MWE back to plague for shorter time
* [```7dc0cb4a```](https://github.com/ktmeaton/NCBImeta/commit/7dc0cb4a) Remove old py config file
* [```5320f671```](https://github.com/ktmeaton/NCBImeta/commit/5320f671) Slowdown fetch as default
* [```c8973be6```](https://github.com/ktmeaton/NCBImeta/commit/c8973be6) Shebang typo
* [```58829cce```](https://github.com/ktmeaton/NCBImeta/commit/58829cce) Implement read handle error checking
* [```035d0fdb```](https://github.com/ktmeaton/NCBImeta/commit/035d0fdb) Error class for read errors
* [```1df91ae4```](https://github.com/ktmeaton/NCBImeta/commit/1df91ae4) add db name
* [```25aff7c9```](https://github.com/ktmeaton/NCBImeta/commit/25aff7c9) Switch example config file to P. aeruginosa to match paper
* [```d50e012d```](https://github.com/ktmeaton/NCBImeta/commit/d50e012d) v0.4.0 update rename
* [```91b3d626```](https://github.com/ktmeaton/NCBImeta/commit/91b3d626) HTTP 429 error checking for efetch
* [```906b10ce```](https://github.com/ktmeaton/NCBImeta/commit/906b10ce) Example command now references config.yaml
* [```99ed046b```](https://github.com/ktmeaton/NCBImeta/commit/99ed046b) Updates for v0.3.4 and v0.3.5
* [```d6483289```](https://github.com/ktmeaton/NCBImeta/commit/d6483289) Implement hiearchical fields
* [```36639a2d```](https://github.com/ktmeaton/NCBImeta/commit/36639a2d) Update hierarchical docs
* [```66ab5c06```](https://github.com/ktmeaton/NCBImeta/commit/66ab5c06) comma space sep hierarchical fields
* [```6da14858```](https://github.com/ktmeaton/NCBImeta/commit/6da14858) Removed extra config files
* [```e4f36d3d```](https://github.com/ktmeaton/NCBImeta/commit/e4f36d3d) Add requirements file
* [```d772064e```](https://github.com/ktmeaton/NCBImeta/commit/d772064e) Implement config file as yaml
* [```91cd02da```](https://github.com/ktmeaton/NCBImeta/commit/91cd02da) Error for incorrect config parameters
* [```b6966d80```](https://github.com/ktmeaton/NCBImeta/commit/b6966d80) config file in yaml format
* [```a16f9699```](https://github.com/ktmeaton/NCBImeta/commit/a16f9699) Merge branch 'master' into dev

## v0.3.4

### Commits

* [```da80053d```](https://github.com/ktmeaton/NCBImeta/commit/da80053d) Update
* [```f53bdac5```](https://github.com/ktmeaton/NCBImeta/commit/f53bdac5) Merge branch 'dev'
* [```297f12fe```](https://github.com/ktmeaton/NCBImeta/commit/297f12fe) Document AnnotateReplace and AnnotateConcatenate
* [```18cf30ab```](https://github.com/ktmeaton/NCBImeta/commit/18cf30ab) Document FORCE_PAUSE_SECONDS parameter
* [```8bf17376```](https://github.com/ktmeaton/NCBImeta/commit/8bf17376) Include all 6 tables in config files
* [```424e4a63```](https://github.com/ktmeaton/NCBImeta/commit/424e4a63) Create a parameter to allow force pausing in seconds
* [```09425f43```](https://github.com/ktmeaton/NCBImeta/commit/09425f43) Implemented NCBI API key acceptance
* [```9927043a```](https://github.com/ktmeaton/NCBImeta/commit/9927043a) Add P. aeruginosa config file
* [```3bc04b38```](https://github.com/ktmeaton/NCBImeta/commit/3bc04b38) Added error when max fetch records is exceeded
* [```5464bc40```](https://github.com/ktmeaton/NCBImeta/commit/5464bc40) Bugfix for HTTPError 429 checking
* [```8bb603e9```](https://github.com/ktmeaton/NCBImeta/commit/8bb603e9) HTTP Error catching for Entrez esummary
* [```aad5f391```](https://github.com/ktmeaton/NCBImeta/commit/aad5f391) shebang change to python3
* [```3ebd8e0c```](https://github.com/ktmeaton/NCBImeta/commit/3ebd8e0c) Set theme jekyll-theme-slate
* [```578e473b```](https://github.com/ktmeaton/NCBImeta/commit/578e473b) Set theme jekyll-theme-minimal
* [```b1981633```](https://github.com/ktmeaton/NCBImeta/commit/b1981633) Goal: Config file as yaml format.

## v0.3.3

### Commits

* [```2c8f1d5e```](https://github.com/ktmeaton/NCBImeta/commit/2c8f1d5e) v0.3.3 Release
* [```b792547b```](https://github.com/ktmeaton/NCBImeta/commit/b792547b) Merge branch 'dev'
* [```27974444```](https://github.com/ktmeaton/NCBImeta/commit/27974444) README update
* [```3c557c6d```](https://github.com/ktmeaton/NCBImeta/commit/3c557c6d) Force python3 command
* [```5e10baa8```](https://github.com/ktmeaton/NCBImeta/commit/5e10baa8) Changelog includes Pubmed support
* [```7e2196cd```](https://github.com/ktmeaton/NCBImeta/commit/7e2196cd) Pubmed add and successful example test
* [```f12a9201```](https://github.com/ktmeaton/NCBImeta/commit/f12a9201) Schema update and Pubmed addition
* [```83e321e6```](https://github.com/ktmeaton/NCBImeta/commit/83e321e6) Pubmed Table functionality added
* [```ea5e6768```](https://github.com/ktmeaton/NCBImeta/commit/ea5e6768) Added a str conversion to xml parsing step for Pubmed table
* [```d08dc500```](https://github.com/ktmeaton/NCBImeta/commit/d08dc500) v0.3.3 development begins

## v0.3.2

### Commits

* [```f800a8a3```](https://github.com/ktmeaton/NCBImeta/commit/f800a8a3) Changelog update
* [```39ab0e2b```](https://github.com/ktmeaton/NCBImeta/commit/39ab0e2b) Delete my_organism_annot.txt
* [```86abec38```](https://github.com/ktmeaton/NCBImeta/commit/86abec38) New annotation files for exampel
* [```8eeb4be8```](https://github.com/ktmeaton/NCBImeta/commit/8eeb4be8) NCBImeta_Export converted to python3
* [```6ec68501```](https://github.com/ktmeaton/NCBImeta/commit/6ec68501) Missing DB Error added
* [```315a3889```](https://github.com/ktmeaton/NCBImeta/commit/315a3889) Last commit before v0.3.3
* [```63141bd9```](https://github.com/ktmeaton/NCBImeta/commit/63141bd9) Nucleotide annotations fully fixed
* [```0d3db343```](https://github.com/ktmeaton/NCBImeta/commit/0d3db343) Nucleotide schema CDSs pluralize
* [```65b03ee5```](https://github.com/ktmeaton/NCBImeta/commit/65b03ee5) Nucleotide Table accession fix
* [```1fe7dc10```](https://github.com/ktmeaton/NCBImeta/commit/1fe7dc10) Fix Nucleotide annotations for py3
* [```36027f7d```](https://github.com/ktmeaton/NCBImeta/commit/36027f7d) Remove utf encode for py3 functionality
* [```c3923ca4```](https://github.com/ktmeaton/NCBImeta/commit/c3923ca4) Unicode decoding for Python3
* [```359a2f7d```](https://github.com/ktmeaton/NCBImeta/commit/359a2f7d) dict_items error fix
* [```7c10cc81```](https://github.com/ktmeaton/NCBImeta/commit/7c10cc81) dev branch reinstated
* [```d7603d54```](https://github.com/ktmeaton/NCBImeta/commit/d7603d54) Add XPath and XLST to wishlist
* [```c2838e84```](https://github.com/ktmeaton/NCBImeta/commit/c2838e84) Update README.md
* [```d21c83ef```](https://github.com/ktmeaton/NCBImeta/commit/d21c83ef) Updated instructions to reflect no dev branch yet
* [```4f37c238```](https://github.com/ktmeaton/NCBImeta/commit/4f37c238) Update development branch names to dev
* [```0472be92```](https://github.com/ktmeaton/NCBImeta/commit/0472be92) Update README.md
* [```689a5452```](https://github.com/ktmeaton/NCBImeta/commit/689a5452) First successful completion of the master join table
* [```a21171da```](https://github.com/ktmeaton/NCBImeta/commit/a21171da) Join script beginning functional
* [```28f47249```](https://github.com/ktmeaton/NCBImeta/commit/28f47249) Major fix in Nucleotide metadata retrieval, properly recovering annotations
* [```22aa44d7```](https://github.com/ktmeaton/NCBImeta/commit/22aa44d7) Little things in README
* [```8d2f0cf9```](https://github.com/ktmeaton/NCBImeta/commit/8d2f0cf9) example update
* [```1a852f1e```](https://github.com/ktmeaton/NCBImeta/commit/1a852f1e) Rename default column names so that they are all unique
* [```e53deb13```](https://github.com/ktmeaton/NCBImeta/commit/e53deb13) Metadata annotation scripts now allow concatenation
* [```5cdd9bb5```](https://github.com/ktmeaton/NCBImeta/commit/5cdd9bb5) Change annotation script to 2 options, replace or concatenate
* [```6f9247c8```](https://github.com/ktmeaton/NCBImeta/commit/6f9247c8) Use harmonized_name attribute for BioSample fields
* [```bffefeef```](https://github.com/ktmeaton/NCBImeta/commit/bffefeef) Remove excess config files and example output
* [```c45fc56c```](https://github.com/ktmeaton/NCBImeta/commit/c45fc56c) Slow down NCBI fetch requests to 1 per second
* [```a0459535```](https://github.com/ktmeaton/NCBImeta/commit/a0459535) Update README.md
* [```c90c644f```](https://github.com/ktmeaton/NCBImeta/commit/c90c644f) Update README.md
* [```b7315650```](https://github.com/ktmeaton/NCBImeta/commit/b7315650) Update README.md
* [```5993154b```](https://github.com/ktmeaton/NCBImeta/commit/5993154b) GIFs
* [```b21cbfde```](https://github.com/ktmeaton/NCBImeta/commit/b21cbfde) Add files via upload
* [```8cd4cbe7```](https://github.com/ktmeaton/NCBImeta/commit/8cd4cbe7) Encode rather than unicode, experiment

## v0.3.1

### Pull Requests

* [```pull/1```](https://github.com/ktmeaton/NCBImeta/pull/1) 1.1

### Commits

* [```cd8fda64```](https://github.com/ktmeaton/NCBImeta/commit/cd8fda64) Update README.md
* [```e6da6252```](https://github.com/ktmeaton/NCBImeta/commit/e6da6252) Update README_schema.md
* [```2d96701b```](https://github.com/ktmeaton/NCBImeta/commit/2d96701b) Update README_schema.md
* [```57ff5348```](https://github.com/ktmeaton/NCBImeta/commit/57ff5348) Update README_schema.md
* [```38e48611```](https://github.com/ktmeaton/NCBImeta/commit/38e48611) Update README_schema.md
* [```13d25e7b```](https://github.com/ktmeaton/NCBImeta/commit/13d25e7b) Update README_schema.md
* [```93b830da```](https://github.com/ktmeaton/NCBImeta/commit/93b830da) Update and rename README_schema.txt to README_schema.md
* [```4b806836```](https://github.com/ktmeaton/NCBImeta/commit/4b806836) Update README_config.md
* [```dc775f57```](https://github.com/ktmeaton/NCBImeta/commit/dc775f57) Update README_config.md
* [```766c815b```](https://github.com/ktmeaton/NCBImeta/commit/766c815b) Update README_config.md
* [```fd2a4f4e```](https://github.com/ktmeaton/NCBImeta/commit/fd2a4f4e) Update README_config.md
* [```328a571b```](https://github.com/ktmeaton/NCBImeta/commit/328a571b) Update README_config.md
* [```ee9eedc4```](https://github.com/ktmeaton/NCBImeta/commit/ee9eedc4) Update and rename README_config.txt to README_config.md
* [```b6d4bac2```](https://github.com/ktmeaton/NCBImeta/commit/b6d4bac2) Update README.md
* [```047b0afa```](https://github.com/ktmeaton/NCBImeta/commit/047b0afa) Move images
* [```e5c3f027```](https://github.com/ktmeaton/NCBImeta/commit/e5c3f027) Moved images
* [```87a24a72```](https://github.com/ktmeaton/NCBImeta/commit/87a24a72) Update README.md
* [```77148b58```](https://github.com/ktmeaton/NCBImeta/commit/77148b58) Update README.md
* [```104ef14a```](https://github.com/ktmeaton/NCBImeta/commit/104ef14a) Update README.md
* [```57349436```](https://github.com/ktmeaton/NCBImeta/commit/57349436) Update README.md
* [```09ec0dc4```](https://github.com/ktmeaton/NCBImeta/commit/09ec0dc4) Update README.md
* [```40858177```](https://github.com/ktmeaton/NCBImeta/commit/40858177) Update README.md
* [```69347c18```](https://github.com/ktmeaton/NCBImeta/commit/69347c18) Update README.md
* [```c3488df1```](https://github.com/ktmeaton/NCBImeta/commit/c3488df1) Update README.md
* [```e3ac31c6```](https://github.com/ktmeaton/NCBImeta/commit/e3ac31c6) Update README.md
* [```3551c725```](https://github.com/ktmeaton/NCBImeta/commit/3551c725) Add files via upload
* [```3963afca```](https://github.com/ktmeaton/NCBImeta/commit/3963afca) Update README.md
* [```d49de82a```](https://github.com/ktmeaton/NCBImeta/commit/d49de82a) Update README.md
* [```25bf63e1```](https://github.com/ktmeaton/NCBImeta/commit/25bf63e1) Add files via upload
* [```ac30fd21```](https://github.com/ktmeaton/NCBImeta/commit/ac30fd21) Delete NCBImeta_screenshot.png
* [```f3afa94e```](https://github.com/ktmeaton/NCBImeta/commit/f3afa94e) Update README.md
* [```3a393d08```](https://github.com/ktmeaton/NCBImeta/commit/3a393d08) higher
* [```5349f1a1```](https://github.com/ktmeaton/NCBImeta/commit/5349f1a1) merge
* [```e2a132bb```](https://github.com/ktmeaton/NCBImeta/commit/e2a132bb) higher res
* [```d65ee345```](https://github.com/ktmeaton/NCBImeta/commit/d65ee345) Delete NCBImeta_snapshot.jpg
* [```82a1350e```](https://github.com/ktmeaton/NCBImeta/commit/82a1350e) Delete NCBImeta_screenshot.png
* [```a15b8041```](https://github.com/ktmeaton/NCBImeta/commit/a15b8041) Update README.md
* [```693012ca```](https://github.com/ktmeaton/NCBImeta/commit/693012ca) Update README.md
* [```f9b1169c```](https://github.com/ktmeaton/NCBImeta/commit/f9b1169c) Update README.md
* [```97e58579```](https://github.com/ktmeaton/NCBImeta/commit/97e58579) Update README.md
* [```489bd94b```](https://github.com/ktmeaton/NCBImeta/commit/489bd94b) Update README.md
* [```01670d9d```](https://github.com/ktmeaton/NCBImeta/commit/01670d9d) Update README.md
* [```7118d274```](https://github.com/ktmeaton/NCBImeta/commit/7118d274) Update README.md
* [```58dce4cd```](https://github.com/ktmeaton/NCBImeta/commit/58dce4cd) Update README.md
* [```dd7c5a0d```](https://github.com/ktmeaton/NCBImeta/commit/dd7c5a0d) Update README.md
* [```ddd78577```](https://github.com/ktmeaton/NCBImeta/commit/ddd78577) Update README.md
* [```df34ed8a```](https://github.com/ktmeaton/NCBImeta/commit/df34ed8a) Update README.md
* [```94d1827e```](https://github.com/ktmeaton/NCBImeta/commit/94d1827e) Update README.md
* [```3bd6f24c```](https://github.com/ktmeaton/NCBImeta/commit/3bd6f24c) Update README.md
* [```c84af731```](https://github.com/ktmeaton/NCBImeta/commit/c84af731) Update README.md
* [```037833f4```](https://github.com/ktmeaton/NCBImeta/commit/037833f4) Update README.md
* [```24a010b8```](https://github.com/ktmeaton/NCBImeta/commit/24a010b8) Update README.md
* [```830ca965```](https://github.com/ktmeaton/NCBImeta/commit/830ca965) Image add
* [```e906acc4```](https://github.com/ktmeaton/NCBImeta/commit/e906acc4) Update README.md
* [```719e986b```](https://github.com/ktmeaton/NCBImeta/commit/719e986b) Update README.md
* [```1f590d11```](https://github.com/ktmeaton/NCBImeta/commit/1f590d11) Update README.md
* [```35eb86ec```](https://github.com/ktmeaton/NCBImeta/commit/35eb86ec) Update README.md
* [```331489ac```](https://github.com/ktmeaton/NCBImeta/commit/331489ac) Update README.md
* [```a6485891```](https://github.com/ktmeaton/NCBImeta/commit/a6485891) Update README.md
* [```a6be4056```](https://github.com/ktmeaton/NCBImeta/commit/a6be4056) Update README.md
* [```18cb81a8```](https://github.com/ktmeaton/NCBImeta/commit/18cb81a8) Update README.md
* [```0ca9547a```](https://github.com/ktmeaton/NCBImeta/commit/0ca9547a) Update README.md
* [```65b290f7```](https://github.com/ktmeaton/NCBImeta/commit/65b290f7) Update README.md
* [```00817806```](https://github.com/ktmeaton/NCBImeta/commit/00817806) Update README.md
* [```ede0ebe4```](https://github.com/ktmeaton/NCBImeta/commit/ede0ebe4) example dir, more unicode testing
* [```13b657b1```](https://github.com/ktmeaton/NCBImeta/commit/13b657b1) Bug fixing more unicode
* [```de6e9bba```](https://github.com/ktmeaton/NCBImeta/commit/de6e9bba) Update NCBImeta.py
* [```c4051788```](https://github.com/ktmeaton/NCBImeta/commit/c4051788) Update README.md
* [```7c136555```](https://github.com/ktmeaton/NCBImeta/commit/7c136555) Update README.md
* [```aa71c05c```](https://github.com/ktmeaton/NCBImeta/commit/aa71c05c) Update README.md
* [```4df6f8f3```](https://github.com/ktmeaton/NCBImeta/commit/4df6f8f3) Update README.md
* [```46041edc```](https://github.com/ktmeaton/NCBImeta/commit/46041edc) Update README.md
* [```fdc00290```](https://github.com/ktmeaton/NCBImeta/commit/fdc00290) Update README.md
* [```9ffcf42a```](https://github.com/ktmeaton/NCBImeta/commit/9ffcf42a) Update README.md
* [```cf734b6f```](https://github.com/ktmeaton/NCBImeta/commit/cf734b6f) Update README.md
* [```e2937a41```](https://github.com/ktmeaton/NCBImeta/commit/e2937a41) Update README.md
* [```95228b11```](https://github.com/ktmeaton/NCBImeta/commit/95228b11) Update README.md
* [```39217244```](https://github.com/ktmeaton/NCBImeta/commit/39217244) Update README.md
* [```e74a81a1```](https://github.com/ktmeaton/NCBImeta/commit/e74a81a1) Update README.md
* [```3fc5925a```](https://github.com/ktmeaton/NCBImeta/commit/3fc5925a) Update README.md
* [```88dd3d8a```](https://github.com/ktmeaton/NCBImeta/commit/88dd3d8a) Update README.md
* [```c4064540```](https://github.com/ktmeaton/NCBImeta/commit/c4064540) Update README.md
* [```611501b9```](https://github.com/ktmeaton/NCBImeta/commit/611501b9) Update README.md
* [```06193bf2```](https://github.com/ktmeaton/NCBImeta/commit/06193bf2) Update README.md
* [```65594c59```](https://github.com/ktmeaton/NCBImeta/commit/65594c59) Update README.md
* [```a33ef410```](https://github.com/ktmeaton/NCBImeta/commit/a33ef410) Update README.md
* [```ee70c20c```](https://github.com/ktmeaton/NCBImeta/commit/ee70c20c) Update README.md
* [```c0ba54c0```](https://github.com/ktmeaton/NCBImeta/commit/c0ba54c0) Update README.md
* [```b9563742```](https://github.com/ktmeaton/NCBImeta/commit/b9563742) Update README.md
* [```f8bd00d9```](https://github.com/ktmeaton/NCBImeta/commit/f8bd00d9) Update README.md
* [```5d504ec8```](https://github.com/ktmeaton/NCBImeta/commit/5d504ec8) Update README.md
* [```a4a7a42d```](https://github.com/ktmeaton/NCBImeta/commit/a4a7a42d) Update README.md
* [```576f7df7```](https://github.com/ktmeaton/NCBImeta/commit/576f7df7) Update README.md
* [```6c3c5f6b```](https://github.com/ktmeaton/NCBImeta/commit/6c3c5f6b) Update README.md
* [```4716032e```](https://github.com/ktmeaton/NCBImeta/commit/4716032e) Update README.md
* [```5360f6db```](https://github.com/ktmeaton/NCBImeta/commit/5360f6db) Badge Testing
* [```c268d121```](https://github.com/ktmeaton/NCBImeta/commit/c268d121) Update geocode.R
* [```21dcccda```](https://github.com/ktmeaton/NCBImeta/commit/21dcccda) Update CHANGELOG.md
* [```995b7027```](https://github.com/ktmeaton/NCBImeta/commit/995b7027) Update CHANGELOG.md
* [```7cdae9c5```](https://github.com/ktmeaton/NCBImeta/commit/7cdae9c5) Update CHANGELOG.md
* [```e2d363de```](https://github.com/ktmeaton/NCBImeta/commit/e2d363de) Merge branch 'master' of https://github.com/ktmeaton/NCBImeta
* [```2b24140c```](https://github.com/ktmeaton/NCBImeta/commit/2b24140c) Remove R tmp
* [```d98765f8```](https://github.com/ktmeaton/NCBImeta/commit/d98765f8) Delete .Rhistory
* [```97df6dc6```](https://github.com/ktmeaton/NCBImeta/commit/97df6dc6) Ignore .pyc
* [```b759f156```](https://github.com/ktmeaton/NCBImeta/commit/b759f156) Cleanup
* [```ab545b4b```](https://github.com/ktmeaton/NCBImeta/commit/ab545b4b) cleanup .pyc
* [```7a133696```](https://github.com/ktmeaton/NCBImeta/commit/7a133696) License and setup
* [```7118d109```](https://github.com/ktmeaton/NCBImeta/commit/7118d109) Set theme jekyll-theme-slate
* [```b946a249```](https://github.com/ktmeaton/NCBImeta/commit/b946a249) Delete _config.yml
* [```17af0ace```](https://github.com/ktmeaton/NCBImeta/commit/17af0ace) Update _config.yml
* [```c9c73316```](https://github.com/ktmeaton/NCBImeta/commit/c9c73316) Set theme jekyll-theme-architect
* [```eb87a53c```](https://github.com/ktmeaton/NCBImeta/commit/eb87a53c) Update
* [```5fa31d50```](https://github.com/ktmeaton/NCBImeta/commit/5fa31d50) Update _config.yml
* [```09b4b07e```](https://github.com/ktmeaton/NCBImeta/commit/09b4b07e) Update _config.yml
* [```2471e051```](https://github.com/ktmeaton/NCBImeta/commit/2471e051) Update _config.yml
* [```63caeaaa```](https://github.com/ktmeaton/NCBImeta/commit/63caeaaa) Update _config.yml
* [```50b00455```](https://github.com/ktmeaton/NCBImeta/commit/50b00455) Initial config
* [```3dd8ba90```](https://github.com/ktmeaton/NCBImeta/commit/3dd8ba90) GH Pages
* [```0d15ea95```](https://github.com/ktmeaton/NCBImeta/commit/0d15ea95) Set theme jekyll-theme-architect
* [```7aa3cf5d```](https://github.com/ktmeaton/NCBImeta/commit/7aa3cf5d) Set theme jekyll-theme-architect
* [```02115ff9```](https://github.com/ktmeaton/NCBImeta/commit/02115ff9) Python 2/3 compatability fix for unicode
* [```e3faf17c```](https://github.com/ktmeaton/NCBImeta/commit/e3faf17c) Update README.md
* [```42d82a16```](https://github.com/ktmeaton/NCBImeta/commit/42d82a16) Update requirements
* [```aed30d83```](https://github.com/ktmeaton/NCBImeta/commit/aed30d83) v0.3.2 Begins
* [```0eb820cb```](https://github.com/ktmeaton/NCBImeta/commit/0eb820cb) Merge branch 'v0.3.1'
* [```7a23c244```](https://github.com/ktmeaton/NCBImeta/commit/7a23c244) Ignore pycache directories
* [```d5701128```](https://github.com/ktmeaton/NCBImeta/commit/d5701128) Ignore
* [```4e90e7f7```](https://github.com/ktmeaton/NCBImeta/commit/4e90e7f7) Cleanup
* [```a606819c```](https://github.com/ktmeaton/NCBImeta/commit/a606819c) Cleanup
* [```c541c9bb```](https://github.com/ktmeaton/NCBImeta/commit/c541c9bb) READMEs, bug fix on error raising
* [```6bddcdef```](https://github.com/ktmeaton/NCBImeta/commit/6bddcdef) remove kits
* [```9d21d465```](https://github.com/ktmeaton/NCBImeta/commit/9d21d465) v0.3.1 begins
* [```1c97b6b4```](https://github.com/ktmeaton/NCBImeta/commit/1c97b6b4) Merge branch 'v0.3.0'
* [```d78b88bf```](https://github.com/ktmeaton/NCBImeta/commit/d78b88bf) merge with master
* [```ed269111```](https://github.com/ktmeaton/NCBImeta/commit/ed269111) Final commit of v0.3.0
* [```df8cedba```](https://github.com/ktmeaton/NCBImeta/commit/df8cedba) Fixed node-attribute confliction
* [```c57bafaf```](https://github.com/ktmeaton/NCBImeta/commit/c57bafaf) Cleanup old source files
* [```b069d4fb```](https://github.com/ktmeaton/NCBImeta/commit/b069d4fb) Proper unicode re-encoding
* [```5f673b7b```](https://github.com/ktmeaton/NCBImeta/commit/5f673b7b) Function Assembly, Biosample, BioProject, SRA
* [```e49648e2```](https://github.com/ktmeaton/NCBImeta/commit/e49648e2) Functional state SRA
* [```ace7c457```](https://github.com/ktmeaton/NCBImeta/commit/ace7c457) Temp save
* [```5693e3ef```](https://github.com/ktmeaton/NCBImeta/commit/5693e3ef) Functional automation of Assembly and BioSample
* [```7790c658```](https://github.com/ktmeaton/NCBImeta/commit/7790c658) Remove some files
* [```700bded3```](https://github.com/ktmeaton/NCBImeta/commit/700bded3) Rename and reorganize
* [```3bf8e94c```](https://github.com/ktmeaton/NCBImeta/commit/3bf8e94c) Dynamic play
* [```f39e9c6c```](https://github.com/ktmeaton/NCBImeta/commit/f39e9c6c) Annotate update
* [```2ff7b251```](https://github.com/ktmeaton/NCBImeta/commit/2ff7b251) Annotations
* [```8cdc64c7```](https://github.com/ktmeaton/NCBImeta/commit/8cdc64c7) Duplicate removal
* [```7edfbc17```](https://github.com/ktmeaton/NCBImeta/commit/7edfbc17) XML frustration
* [```26b50189```](https://github.com/ktmeaton/NCBImeta/commit/26b50189) On the fly, database creation
* [```27eb8a32```](https://github.com/ktmeaton/NCBImeta/commit/27eb8a32) v0.3.0 rewrite main script
* [```efa56cee```](https://github.com/ktmeaton/NCBImeta/commit/efa56cee) v0.3.0 Development Beings - Automated
* [```00b6f1ba```](https://github.com/ktmeaton/NCBImeta/commit/00b6f1ba) v0.3.0 Development Begins - Automated
* [```7e77ac20```](https://github.com/ktmeaton/NCBImeta/commit/7e77ac20) Merge branch 'v0.2.1'
* [```9383f6c5```](https://github.com/ktmeaton/NCBImeta/commit/9383f6c5) More annotations
* [```97f6419c```](https://github.com/ktmeaton/NCBImeta/commit/97f6419c) More annotation
* [```49654725```](https://github.com/ktmeaton/NCBImeta/commit/49654725) Massive annotation
* [```a09d4012```](https://github.com/ktmeaton/NCBImeta/commit/a09d4012) CURATED annotation files
* [```5df0f48e```](https://github.com/ktmeaton/NCBImeta/commit/5df0f48e) bioproject addition to annotate
* [```4b3a27f3```](https://github.com/ktmeaton/NCBImeta/commit/4b3a27f3) Check multiple matches in annotation, removed BioSample columns
* [```82661a4b```](https://github.com/ktmeaton/NCBImeta/commit/82661a4b) Consensus file after plotting
* [```710d6b57```](https://github.com/ktmeaton/NCBImeta/commit/710d6b57) First pass at plotting geographic distribution, Pandemic, PubUse
* [```d851f704```](https://github.com/ktmeaton/NCBImeta/commit/d851f704) Plotting
* [```fc2a0712```](https://github.com/ktmeaton/NCBImeta/commit/fc2a0712) Plotting
* [```1248978b```](https://github.com/ktmeaton/NCBImeta/commit/1248978b) Consensus lat,long,pandemic,established
* [```2716c6ba```](https://github.com/ktmeaton/NCBImeta/commit/2716c6ba) Geocoding switched to Doogal
* [```d8161aee```](https://github.com/ktmeaton/NCBImeta/commit/d8161aee) geocoding
* [```9a59fd28```](https://github.com/ktmeaton/NCBImeta/commit/9a59fd28) Yp database for plotting
* [```e5036c05```](https://github.com/ktmeaton/NCBImeta/commit/e5036c05) More annotations
* [```69a3865a```](https://github.com/ktmeaton/NCBImeta/commit/69a3865a) Annotation files
* [```1f095613```](https://github.com/ktmeaton/NCBImeta/commit/1f095613) Dare I say.. functional v0.2?
* [```2d4399d1```](https://github.com/ktmeaton/NCBImeta/commit/2d4399d1) Heuristisc approach to merge
* [```fcce7556```](https://github.com/ktmeaton/NCBImeta/commit/fcce7556) Biosample improve
* [```4450cd1a```](https://github.com/ktmeaton/NCBImeta/commit/4450cd1a) README Title Fix
* [```665b47c3```](https://github.com/ktmeaton/NCBImeta/commit/665b47c3) v0.2.1 README Update
* [```64467aee```](https://github.com/ktmeaton/NCBImeta/commit/64467aee) Assembly, SRA, Bioproject functional
* [```d33fc246```](https://github.com/ktmeaton/NCBImeta/commit/d33fc246) Update CHANGELOG.md
* [```38be909e```](https://github.com/ktmeaton/NCBImeta/commit/38be909e) v0.2.0 In Development Begins, Rename NCBInfect
* [```b0103143```](https://github.com/ktmeaton/NCBImeta/commit/b0103143) v0.2.0 Begins Development
* [```9dfe4d7e```](https://github.com/ktmeaton/NCBImeta/commit/9dfe4d7e) Merge branch 'v0.1.2'
* [```594b4bc4```](https://github.com/ktmeaton/NCBImeta/commit/594b4bc4) Sanity commit, Jan 2018
* [```0d6b9312```](https://github.com/ktmeaton/NCBImeta/commit/0d6b9312) Update CHANGELOG.md
* [```ae8d3b31```](https://github.com/ktmeaton/NCBImeta/commit/ae8d3b31) Update CHANGELOG.md
* [```ee485710```](https://github.com/ktmeaton/NCBImeta/commit/ee485710) Update CHANGELOG.md
* [```0e3e2f33```](https://github.com/ktmeaton/NCBImeta/commit/0e3e2f33) v0.1.2 changes
* [```973d416f```](https://github.com/ktmeaton/NCBImeta/commit/973d416f) Update README.md
* [```b7326d9a```](https://github.com/ktmeaton/NCBImeta/commit/b7326d9a) Update README.md
* [```6e1a7000```](https://github.com/ktmeaton/NCBImeta/commit/6e1a7000) Create CHANGELOG.md
* [```08053019```](https://github.com/ktmeaton/NCBImeta/commit/08053019) Moved annotate file
* [```54f194f7```](https://github.com/ktmeaton/NCBImeta/commit/54f194f7) Merge branch 'master' of https://github.com/ktmeaton/GenomeCollector
* [```b8af9498```](https://github.com/ktmeaton/NCBImeta/commit/b8af9498) Big Changes
* [```43c68690```](https://github.com/ktmeaton/NCBImeta/commit/43c68690) Annotate File
* [```3e0e006e```](https://github.com/ktmeaton/NCBImeta/commit/3e0e006e) python3 version fix
* [```e8e40418```](https://github.com/ktmeaton/NCBImeta/commit/e8e40418) python3 print statement fix
* [```a80db3fd```](https://github.com/ktmeaton/NCBImeta/commit/a80db3fd) Added extra fields in SRA table
* [```5453c8dc```](https://github.com/ktmeaton/NCBImeta/commit/5453c8dc) Version 2.0 fully functional
* [```2a1c06e0```](https://github.com/ktmeaton/NCBImeta/commit/2a1c06e0) SRA functionality as separate module
* [```7280508a```](https://github.com/ktmeaton/NCBImeta/commit/7280508a) Partially working version of SRA Table
* [```5c8f3e98```](https://github.com/ktmeaton/NCBImeta/commit/5c8f3e98) Added SRA function
* [```85654bcf```](https://github.com/ktmeaton/NCBImeta/commit/85654bcf) Added more to do in README.md
* [```2e859b11```](https://github.com/ktmeaton/NCBImeta/commit/2e859b11) Added genomeannotator log file writing
* [```c6ecaf04```](https://github.com/ktmeaton/NCBImeta/commit/c6ecaf04) Added the GenomeAnnotator tool
* [```5de9f53f```](https://github.com/ktmeaton/NCBImeta/commit/5de9f53f) Removed tmp file writing for record_dict testing
* [```fdbe0768```](https://github.com/ktmeaton/NCBImeta/commit/fdbe0768) Changed biosample Entrez parsing to validate=false
* [```bac6eb7a```](https://github.com/ktmeaton/NCBImeta/commit/bac6eb7a) Field AsmReleaseDate replaced with SubmissionDate
* [```642c2a2a```](https://github.com/ktmeaton/NCBImeta/commit/642c2a2a) Updated README.md to version 2.0
* [```506b9768```](https://github.com/ktmeaton/NCBImeta/commit/506b9768) Merge pull request #1 from ktmeaton/1.1
* [```fdc45fb2```](https://github.com/ktmeaton/NCBImeta/commit/fdc45fb2) Update README.md
* [```cde62fab```](https://github.com/ktmeaton/NCBImeta/commit/cde62fab) Update README.md
* [```a521f23f```](https://github.com/ktmeaton/NCBImeta/commit/a521f23f) Update README.md
* [```fc5bea46```](https://github.com/ktmeaton/NCBImeta/commit/fc5bea46) Update genomecollector.py
* [```c433ef22```](https://github.com/ktmeaton/NCBImeta/commit/c433ef22) Update genomecollector.py
* [```7736c0ea```](https://github.com/ktmeaton/NCBImeta/commit/7736c0ea) Add files via upload
* [```dfff42c1```](https://github.com/ktmeaton/NCBImeta/commit/dfff42c1) Delete genomeconcatenate.py
* [```ac643db0```](https://github.com/ktmeaton/NCBImeta/commit/ac643db0) Add files via upload
