# CHANGELOG

## Development

### Commits

* [```3d79271```](https://github.com/ktmeaton/NCBImeta/commit/3d79271) test numpy-1.19.5
* [```3bbbe59```](https://github.com/ktmeaton/NCBImeta/commit/3bbbe59) bugfix in example workflow trigger
* [```05d1960```](https://github.com/ktmeaton/NCBImeta/commit/05d1960) remove print statement
* [```d302811```](https://github.com/ktmeaton/NCBImeta/commit/d302811) update dependencies lxml pyaml and numpy
* [```d5eaf85```](https://github.com/ktmeaton/NCBImeta/commit/d5eaf85) fix bug in test_export biosample date
* [```4a1d4a7```](https://github.com/ktmeaton/NCBImeta/commit/4a1d4a7) add python 3.9 to the testing matrix
* [```31bfcaf```](https://github.com/ktmeaton/NCBImeta/commit/31bfcaf) update ver number to v0.8.1
* [```4272d9e```](https://github.com/ktmeaton/NCBImeta/commit/4272d9e) update changelog for v0.8.0

## v0.8.0

### Notes

1. Start new dev branch.
1. Update miniconda actions and use mamba.
1. Update lxml for security vulnerability.
1. Add autologs as a submodule.
1. Create CHANGELOG with autologs.
1. Simplify test config with less fields to check.
1. Remove .py extension from executable scripts.

### Commits

* [```3b3d630```](https://github.com/ktmeaton/NCBImeta/commit/3b3d630) update autologs for branch and tag link
* [```b5ee965```](https://github.com/ktmeaton/NCBImeta/commit/b5ee965) update submodules for release
* [```2290cfc```](https://github.com/ktmeaton/NCBImeta/commit/2290cfc) restrict testing workflows to master and dev
* [```3efed45```](https://github.com/ktmeaton/NCBImeta/commit/3efed45) Merge pull request #19 from ktmeaton/dev
* [```b325be9```](https://github.com/ktmeaton/NCBImeta/commit/b325be9) disable fail fast and restrict codecov upload
* [```5cb358e```](https://github.com/ktmeaton/NCBImeta/commit/5cb358e) workflow overhaul
* [```cc99998```](https://github.com/ktmeaton/NCBImeta/commit/cc99998) restrict python versions to >=3.6,<3.9
* [```4eaa9bc```](https://github.com/ktmeaton/NCBImeta/commit/4eaa9bc) overhaul test workflow
* [```c174036```](https://github.com/ktmeaton/NCBImeta/commit/c174036) overhaul build workflow
* [```cb0ba1f```](https://github.com/ktmeaton/NCBImeta/commit/cb0ba1f) add execute permissions to the newly renamed files
* [```0a519e0```](https://github.com/ktmeaton/NCBImeta/commit/0a519e0) replace .py extensions for Utilities script
* [```bf4f608```](https://github.com/ktmeaton/NCBImeta/commit/bf4f608) update autologs
* [```db06fd2```](https://github.com/ktmeaton/NCBImeta/commit/db06fd2) update to v0.8.0
* [```6b85877```](https://github.com/ktmeaton/NCBImeta/commit/6b85877) Remove .py extension from executable scripts
* [```c73cb41```](https://github.com/ktmeaton/NCBImeta/commit/c73cb41) Update dev changelog
* [```7c22198```](https://github.com/ktmeaton/NCBImeta/commit/7c22198) simplify database columns for export
* [```9aacd23```](https://github.com/ktmeaton/NCBImeta/commit/9aacd23) add autolog notes for v0.7.0
* [```34ba549```](https://github.com/ktmeaton/NCBImeta/commit/34ba549) update lxml for security fixes
* [```2e203c1```](https://github.com/ktmeaton/NCBImeta/commit/2e203c1) update miniconda action and use mamba for macos
* [```78b4229```](https://github.com/ktmeaton/NCBImeta/commit/78b4229) lint new mammuthus config
* [```6dd43d8```](https://github.com/ktmeaton/NCBImeta/commit/6dd43d8) update miniconda action and use mamba
* [```d275545```](https://github.com/ktmeaton/NCBImeta/commit/d275545) start new dev branch for v0.7.1
* [```c96c186```](https://github.com/ktmeaton/NCBImeta/commit/c96c186) resolve merges before v0.7.1
* [```42239bc```](https://github.com/ktmeaton/NCBImeta/commit/42239bc) add mammuthus config
* [```97fcd95```](https://github.com/ktmeaton/NCBImeta/commit/97fcd95) add new logo to README
* [```e466196```](https://github.com/ktmeaton/NCBImeta/commit/e466196) add new logo

## v0.7.0

### Notes

1. Allow config parameters to be specified at run-time.
1. Restrict biopython to >=1.74,<1.77 because of Issue #13
1. Specify versions for all user and dev dependencies.
1. Updated and moved PR template.
1. Updated Contributor's Guideline
1. Remove Python 3.5 support because of incompatibility with black.
1. Remove user email and API key from stdout.

### Commits

* [```18bf0c8```](https://github.com/ktmeaton/NCBImeta/commit/18bf0c8) remove dev suffix from version number
* [```fc0a25d```](https://github.com/ktmeaton/NCBImeta/commit/fc0a25d) Merge pull request #16 from ktmeaton/cli-params
* [```1592e45```](https://github.com/ktmeaton/NCBImeta/commit/1592e45) update help output
* [```2823318```](https://github.com/ktmeaton/NCBImeta/commit/2823318) add tests for cli param and bad api key
* [```2710af1```](https://github.com/ktmeaton/NCBImeta/commit/2710af1) cli parameters for NCBI api
* [```b5d4811```](https://github.com/ktmeaton/NCBImeta/commit/b5d4811) reorganize community and additional contributors
* [```34cf225```](https://github.com/ktmeaton/NCBImeta/commit/34cf225) update docs with final changes from typeerror
* [```adb51d7```](https://github.com/ktmeaton/NCBImeta/commit/adb51d7) Merge pull request #14 from ktmeaton/typeerror
* [```2baf5ba```](https://github.com/ktmeaton/NCBImeta/commit/2baf5ba) update the PR checklist
* [```1faf1bc```](https://github.com/ktmeaton/NCBImeta/commit/1faf1bc) change ver to 0.7.0dev and restrict all ver
* [```894f428```](https://github.com/ktmeaton/NCBImeta/commit/894f428) restrict pytest and cov ver
* [```8c490bc```](https://github.com/ktmeaton/NCBImeta/commit/8c490bc) fix lint workflow desc
* [```6dd0033```](https://github.com/ktmeaton/NCBImeta/commit/6dd0033) remove python 3.5 support
* [```d3a5a06```](https://github.com/ktmeaton/NCBImeta/commit/d3a5a06) update changelog for v0.6.7dev
* [```e9ce2d0```](https://github.com/ktmeaton/NCBImeta/commit/e9ce2d0) fix build ver broken link
* [```ffe033b```](https://github.com/ktmeaton/NCBImeta/commit/ffe033b) merge pr changes for lint workflow
* [```4146f54```](https://github.com/ktmeaton/NCBImeta/commit/4146f54) Merge branch 'dev' into typeerror
* [```5e69d66```](https://github.com/ktmeaton/NCBImeta/commit/5e69d66) change PR template
* [```bf2e903```](https://github.com/ktmeaton/NCBImeta/commit/bf2e903) restrict pre-commit to <= 2.6.0
* [```6da67ac```](https://github.com/ktmeaton/NCBImeta/commit/6da67ac) reformat workflows and use pip dev install
* [```fa35c1e```](https://github.com/ktmeaton/NCBImeta/commit/fa35c1e) restrict biopython ver to <1.77
* [```b637930```](https://github.com/ktmeaton/NCBImeta/commit/b637930) update to 0.6.7dev
* [```dca6bb2```](https://github.com/ktmeaton/NCBImeta/commit/dca6bb2) remove outdated project ver number
* [```f7d748f```](https://github.com/ktmeaton/NCBImeta/commit/f7d748f) add python support and black linting
* [```aa1f84c```](https://github.com/ktmeaton/NCBImeta/commit/aa1f84c) update test table assembly values
* [```b709450```](https://github.com/ktmeaton/NCBImeta/commit/b709450) add black and flake8 to dev dependencies
* [```4129006```](https://github.com/ktmeaton/NCBImeta/commit/4129006) fix list numbering indentation
* [```f302cce```](https://github.com/ktmeaton/NCBImeta/commit/f302cce) add logo asset attribution
* [```5acd4b7```](https://github.com/ktmeaton/NCBImeta/commit/5acd4b7) add draft logo
* [```d9ce045```](https://github.com/ktmeaton/NCBImeta/commit/d9ce045) add logo draft
* [```2723e18```](https://github.com/ktmeaton/NCBImeta/commit/2723e18) add url link for ncbi ref
* [```8677627```](https://github.com/ktmeaton/NCBImeta/commit/8677627) note add contributors
* [```d9b4723```](https://github.com/ktmeaton/NCBImeta/commit/d9b4723) conda link

## v0.6.6.post1

### Commits

* [```e931c8b```](https://github.com/ktmeaton/NCBImeta/commit/e931c8b) test post suffix again
* [```c4b682b```](https://github.com/ktmeaton/NCBImeta/commit/c4b682b) update to a1 ver
* [```5aafdbc```](https://github.com/ktmeaton/NCBImeta/commit/5aafdbc) run example workflow when src py is changed
* [```ffe86dd```](https://github.com/ktmeaton/NCBImeta/commit/ffe86dd) just increment the minor ver for bug fix
* [```d4e906b```](https://github.com/ktmeaton/NCBImeta/commit/d4e906b) commit to post suffix after test
* [```5be6eb6```](https://github.com/ktmeaton/NCBImeta/commit/5be6eb6) update ver use the post suffix
* [```7103b41```](https://github.com/ktmeaton/NCBImeta/commit/7103b41) remove tag hook causes problems for pypi dup
* [```32dcb8a```](https://github.com/ktmeaton/NCBImeta/commit/32dcb8a) update ver
* [```cb2997e```](https://github.com/ktmeaton/NCBImeta/commit/cb2997e) delete bad position tags
* [```047e96b```](https://github.com/ktmeaton/NCBImeta/commit/047e96b) move tag hook under push action
* [```93d6b81```](https://github.com/ktmeaton/NCBImeta/commit/93d6b81) also run workflows on tags
* [```72c0585```](https://github.com/ktmeaton/NCBImeta/commit/72c0585) update to ver v0.6.7a
* [```43747be```](https://github.com/ktmeaton/NCBImeta/commit/43747be) update pip and setuptools
* [```6116cdb```](https://github.com/ktmeaton/NCBImeta/commit/6116cdb) remove mention of setuptools install
* [```792e28c```](https://github.com/ktmeaton/NCBImeta/commit/792e28c) mark project v0.6.6 as RELEASED
* [```9582d78```](https://github.com/ktmeaton/NCBImeta/commit/9582d78) make pypi install plain no ver

## v0.6.6

### Commits

* [```739cad3```](https://github.com/ktmeaton/NCBImeta/commit/739cad3) switch release trigger to published
* [```c68e540```](https://github.com/ktmeaton/NCBImeta/commit/c68e540) run all workflows on release
* [```ceb5b36```](https://github.com/ktmeaton/NCBImeta/commit/ceb5b36) fix github links
* [```b5681ab```](https://github.com/ktmeaton/NCBImeta/commit/b5681ab) Merge pull request #12 from ktmeaton/dev
* [```1d40f9a```](https://github.com/ktmeaton/NCBImeta/commit/1d40f9a) change ver back to v0.6.6
* [```6d99f14```](https://github.com/ktmeaton/NCBImeta/commit/6d99f14) update links
* [```1d939de```](https://github.com/ktmeaton/NCBImeta/commit/1d939de) update pypi workflows
* [```6f383e3```](https://github.com/ktmeaton/NCBImeta/commit/6f383e3) update version numbers with dev suffix
* [```1c30f04```](https://github.com/ktmeaton/NCBImeta/commit/1c30f04) change version to v0.6.6dev
* [```74dffef```](https://github.com/ktmeaton/NCBImeta/commit/74dffef) remove test dir from package
* [```89875fa```](https://github.com/ktmeaton/NCBImeta/commit/89875fa) tidy up pull request template
* [```ebceb16```](https://github.com/ktmeaton/NCBImeta/commit/ebceb16) tidy up pull request template
* [```c0cc2e2```](https://github.com/ktmeaton/NCBImeta/commit/c0cc2e2) force linters to wait for pre-commit
* [```dcf8bf8```](https://github.com/ktmeaton/NCBImeta/commit/dcf8bf8) change flake8 py find file
* [```94762f7```](https://github.com/ktmeaton/NCBImeta/commit/94762f7) trial yaml lint
* [```90e1087```](https://github.com/ktmeaton/NCBImeta/commit/90e1087) contributors guide file update
* [```6add9ae```](https://github.com/ktmeaton/NCBImeta/commit/6add9ae) contributors guide file update
* [```f11e68f```](https://github.com/ktmeaton/NCBImeta/commit/f11e68f) install notes and project update
* [```5a515e4```](https://github.com/ktmeaton/NCBImeta/commit/5a515e4) separate install and run for flake8
* [```764cbac```](https://github.com/ktmeaton/NCBImeta/commit/764cbac) add back in flake8
* [```69d8a7c```](https://github.com/ktmeaton/NCBImeta/commit/69d8a7c) add back in python linting
* [```d96c0d0```](https://github.com/ktmeaton/NCBImeta/commit/d96c0d0) remove python linting
* [```4ad65ae```](https://github.com/ktmeaton/NCBImeta/commit/4ad65ae) remove python linting
* [```445c69a```](https://github.com/ktmeaton/NCBImeta/commit/445c69a) remove python linting
* [```8f51169```](https://github.com/ktmeaton/NCBImeta/commit/8f51169) just test flake8 install
* [```6637634```](https://github.com/ktmeaton/NCBImeta/commit/6637634) temp remove pr hook and test flake8
* [```522e124```](https://github.com/ktmeaton/NCBImeta/commit/522e124) add macos conda install test
* [```56b3ec1```](https://github.com/ktmeaton/NCBImeta/commit/56b3ec1) add macos build workflow
* [```2750e06```](https://github.com/ktmeaton/NCBImeta/commit/2750e06) use conda setup action
* [```36b88da```](https://github.com/ktmeaton/NCBImeta/commit/36b88da) comment out PyPI repo install for now
* [```48b15d5```](https://github.com/ktmeaton/NCBImeta/commit/48b15d5) remove setuptools install
* [```d44db94```](https://github.com/ktmeaton/NCBImeta/commit/d44db94) MAJOR bug fix for build workflow
* [```86c9065```](https://github.com/ktmeaton/NCBImeta/commit/86c9065) try to force build workflow with python3
* [```d7c1f9f```](https://github.com/ktmeaton/NCBImeta/commit/d7c1f9f) test lint workflow no req
* [```9870d23```](https://github.com/ktmeaton/NCBImeta/commit/9870d23) update markdownlint-cli to ignore paper dir
* [```062928c```](https://github.com/ktmeaton/NCBImeta/commit/062928c) lint setup
* [```c684f97```](https://github.com/ktmeaton/NCBImeta/commit/c684f97) lint bug_report
* [```f7a62b5```](https://github.com/ktmeaton/NCBImeta/commit/f7a62b5) new pull request template
* [```9cdc4ce```](https://github.com/ktmeaton/NCBImeta/commit/9cdc4ce) more detailed contibute guide
* [```0c1fcd1```](https://github.com/ktmeaton/NCBImeta/commit/0c1fcd1) more paths to lint
* [```12078ef```](https://github.com/ktmeaton/NCBImeta/commit/12078ef) lint conftest
* [```4976309```](https://github.com/ktmeaton/NCBImeta/commit/4976309) lint schema README
* [```cd724f9```](https://github.com/ktmeaton/NCBImeta/commit/cd724f9) lint config README
* [```93035ff```](https://github.com/ktmeaton/NCBImeta/commit/93035ff) update linting workflow
* [```aab8b90```](https://github.com/ktmeaton/NCBImeta/commit/aab8b90) lint conftest max 80
* [```3c4f9a8```](https://github.com/ktmeaton/NCBImeta/commit/3c4f9a8) lint conftest max 80
* [```7d6217d```](https://github.com/ktmeaton/NCBImeta/commit/7d6217d) lint test_ncbimeta max 80
* [```f8952f5```](https://github.com/ktmeaton/NCBImeta/commit/f8952f5) lint test_ncbimeta max 80
* [```267b401```](https://github.com/ktmeaton/NCBImeta/commit/267b401) lint test_xml max 80
* [```534f397```](https://github.com/ktmeaton/NCBImeta/commit/534f397) lint test_utilities max 80
* [```5cc8065```](https://github.com/ktmeaton/NCBImeta/commit/5cc8065) lint test_utilities max 80
* [```499c75f```](https://github.com/ktmeaton/NCBImeta/commit/499c75f) lint test_join max 80
* [```1dd3536```](https://github.com/ktmeaton/NCBImeta/commit/1dd3536) lint test_annotate max 80
* [```0693806```](https://github.com/ktmeaton/NCBImeta/commit/0693806) lint test_export max 80
* [```c8c67cc```](https://github.com/ktmeaton/NCBImeta/commit/c8c67cc) lint test_errors
* [```db12ca1```](https://github.com/ktmeaton/NCBImeta/commit/db12ca1) lint NCBImetaAnnotate max 80
* [```b8030ef```](https://github.com/ktmeaton/NCBImeta/commit/b8030ef) lint NCBImetaJoin max 80
* [```6b78c19```](https://github.com/ktmeaton/NCBImeta/commit/6b78c19) lint NCBImetaJoin max 80
* [```49a4e91```](https://github.com/ktmeaton/NCBImeta/commit/49a4e91) lint NCBImetaErrors max 80
* [```8e51118```](https://github.com/ktmeaton/NCBImeta/commit/8e51118) lint NCBImetaErrors max 80
* [```1c90912```](https://github.com/ktmeaton/NCBImeta/commit/1c90912) redo py linting max 80
* [```3a912e5```](https://github.com/ktmeaton/NCBImeta/commit/3a912e5) lint NCBImetaUtilities max 80
* [```69abba6```](https://github.com/ktmeaton/NCBImeta/commit/69abba6) lint NCBImetaUtilities max 80
* [```c635f83```](https://github.com/ktmeaton/NCBImeta/commit/c635f83) lint NCBImetaUtilities max 80
* [```f0818fd```](https://github.com/ktmeaton/NCBImeta/commit/f0818fd) lint NCBImetaUtilities max 80
* [```4997e76```](https://github.com/ktmeaton/NCBImeta/commit/4997e76) lint NCBImetaUtilities max 80
* [```d2c2f6d```](https://github.com/ktmeaton/NCBImeta/commit/d2c2f6d) lint NCBImetaExport max 80
* [```6cc3f06```](https://github.com/ktmeaton/NCBImeta/commit/6cc3f06) lint NCBImetaExport max 80
* [```8f3f898```](https://github.com/ktmeaton/NCBImeta/commit/8f3f898) update badges and requirements
* [```5432727```](https://github.com/ktmeaton/NCBImeta/commit/5432727) more informative names for build jobs
* [```229fc56```](https://github.com/ktmeaton/NCBImeta/commit/229fc56) lint CHANGELOG with line len
* [```35a6482```](https://github.com/ktmeaton/NCBImeta/commit/35a6482) update build badge to gh actions
* [```64f5145```](https://github.com/ktmeaton/NCBImeta/commit/64f5145) simplify workflow names
* [```9625959```](https://github.com/ktmeaton/NCBImeta/commit/9625959) simplify workflow names
* [```8bf6dbe```](https://github.com/ktmeaton/NCBImeta/commit/8bf6dbe) lint CHANGELOG
* [```3158f78```](https://github.com/ktmeaton/NCBImeta/commit/3158f78) add markdownlint-cli to pre-commit lint
* [```863aaf2```](https://github.com/ktmeaton/NCBImeta/commit/863aaf2) ncbimeta adhere to B950 max 80
* [```f6e58b2```](https://github.com/ktmeaton/NCBImeta/commit/f6e58b2) ignore B950 in conftest
* [```a609a1e```](https://github.com/ktmeaton/NCBImeta/commit/a609a1e) fix i var in loop
* [```212fd6e```](https://github.com/ktmeaton/NCBImeta/commit/212fd6e) revert to len 80
* [```eb1b9f0```](https://github.com/ktmeaton/NCBImeta/commit/eb1b9f0) revert to len 80
* [```a3c2b93```](https://github.com/ktmeaton/NCBImeta/commit/a3c2b93) lint ncbimeta with bugbear
* [```abdb2eb```](https://github.com/ktmeaton/NCBImeta/commit/abdb2eb) lint ncbimeta with bugbear
* [```b223578```](https://github.com/ktmeaton/NCBImeta/commit/b223578) lint ncbimeta with bugbear
* [```f60b425```](https://github.com/ktmeaton/NCBImeta/commit/f60b425) lint ncbimeta with bugbear
* [```ec29ac9```](https://github.com/ktmeaton/NCBImeta/commit/ec29ac9) lint ncbimeta with bugbear
* [```0e1b519```](https://github.com/ktmeaton/NCBImeta/commit/0e1b519) lint ncbimeta
* [```a790dd6```](https://github.com/ktmeaton/NCBImeta/commit/a790dd6) attempte bugbear for flake8
* [```3e93d0f```](https://github.com/ktmeaton/NCBImeta/commit/3e93d0f) more str format for black
* [```37beb31```](https://github.com/ktmeaton/NCBImeta/commit/37beb31) reformat str to make black happy
* [```caf0312```](https://github.com/ktmeaton/NCBImeta/commit/caf0312) lint with new rules
* [```15b5e41```](https://github.com/ktmeaton/NCBImeta/commit/15b5e41) lint with new rules
* [```a7a8e32```](https://github.com/ktmeaton/NCBImeta/commit/a7a8e32) lint with new rules
* [```e63ddc9```](https://github.com/ktmeaton/NCBImeta/commit/e63ddc9) lint with new rules
* [```048855c```](https://github.com/ktmeaton/NCBImeta/commit/048855c) lint with new rules
* [```f69c9bf```](https://github.com/ktmeaton/NCBImeta/commit/f69c9bf) lint with new rules
* [```9ccea9c```](https://github.com/ktmeaton/NCBImeta/commit/9ccea9c) fix test module import
* [```55da36c```](https://github.com/ktmeaton/NCBImeta/commit/55da36c) lint conftest
* [```e41452e```](https://github.com/ktmeaton/NCBImeta/commit/e41452e) test file specific rule exclude
* [```cad97a4```](https://github.com/ktmeaton/NCBImeta/commit/cad97a4) change E501 exclude to conftest
* [```294bc50```](https://github.com/ktmeaton/NCBImeta/commit/294bc50) flake8 config in setup.cfg
* [```ab23797```](https://github.com/ktmeaton/NCBImeta/commit/ab23797) test flake8 ignore conftest
* [```4563f59```](https://github.com/ktmeaton/NCBImeta/commit/4563f59) lint test_ncbimeta
* [```9fa14aa```](https://github.com/ktmeaton/NCBImeta/commit/9fa14aa) lint test_ncbimeta
* [```b22f491```](https://github.com/ktmeaton/NCBImeta/commit/b22f491) lint test_utilities
* [```2be3896```](https://github.com/ktmeaton/NCBImeta/commit/2be3896) lint test_utilities
* [```082f915```](https://github.com/ktmeaton/NCBImeta/commit/082f915) lint test_join
* [```b5742fd```](https://github.com/ktmeaton/NCBImeta/commit/b5742fd) lint test_export
* [```7d4f29c```](https://github.com/ktmeaton/NCBImeta/commit/7d4f29c) lint test_join
* [```efa7c75```](https://github.com/ktmeaton/NCBImeta/commit/efa7c75) lint test_export
* [```93cffa9```](https://github.com/ktmeaton/NCBImeta/commit/93cffa9) lint test_errors
* [```c25df97```](https://github.com/ktmeaton/NCBImeta/commit/c25df97) lint test_errors
* [```cde103f```](https://github.com/ktmeaton/NCBImeta/commit/cde103f) fix comment spacing
* [```5d57430```](https://github.com/ktmeaton/NCBImeta/commit/5d57430) fix comment spacing
* [```25ee239```](https://github.com/ktmeaton/NCBImeta/commit/25ee239) lint test_xml
* [```000d906```](https://github.com/ktmeaton/NCBImeta/commit/000d906) lint test_xml
* [```810a70c```](https://github.com/ktmeaton/NCBImeta/commit/810a70c) lint test_annotateconcatenate
* [```c76738a```](https://github.com/ktmeaton/NCBImeta/commit/c76738a) lint test_annotateconcatenate
* [```04dd099```](https://github.com/ktmeaton/NCBImeta/commit/04dd099) lint test_annotatereplace
* [```ce0f9e0```](https://github.com/ktmeaton/NCBImeta/commit/ce0f9e0) lint test_annotatereplace
* [```a6acbc0```](https://github.com/ktmeaton/NCBImeta/commit/a6acbc0) lint test_annotatereplace
* [```f039442```](https://github.com/ktmeaton/NCBImeta/commit/f039442) lint NCBImeta
* [```489fc1b```](https://github.com/ktmeaton/NCBImeta/commit/489fc1b) lint NCBImeta
* [```cf980d1```](https://github.com/ktmeaton/NCBImeta/commit/cf980d1) lint NCBImetaAnnotate
* [```2829ebf```](https://github.com/ktmeaton/NCBImeta/commit/2829ebf) lint NCBImetaAnnotate
* [```8388eb7```](https://github.com/ktmeaton/NCBImeta/commit/8388eb7) lint NCBImetaAnnotateReplace
* [```4ccfc22```](https://github.com/ktmeaton/NCBImeta/commit/4ccfc22) lint NCBImetaExport
* [```c2aab36```](https://github.com/ktmeaton/NCBImeta/commit/c2aab36) lint NCBImetaExport
* [```b00abcd```](https://github.com/ktmeaton/NCBImeta/commit/b00abcd) add lint checklist for py files
* [```cc8243d```](https://github.com/ktmeaton/NCBImeta/commit/cc8243d) lint NCBImetaJoin
* [```7614502```](https://github.com/ktmeaton/NCBImeta/commit/7614502) lint single quotes
* [```3fe5436```](https://github.com/ktmeaton/NCBImeta/commit/3fe5436) pull request template
* [```6d2c1b7```](https://github.com/ktmeaton/NCBImeta/commit/6d2c1b7) remove old dev require file
* [```58b0669```](https://github.com/ktmeaton/NCBImeta/commit/58b0669) draft contributors guide
* [```6fdf6f2```](https://github.com/ktmeaton/NCBImeta/commit/6fdf6f2) update build comments
* [```348dd6d```](https://github.com/ktmeaton/NCBImeta/commit/348dd6d) linting and dev depend update
* [```854a437```](https://github.com/ktmeaton/NCBImeta/commit/854a437) remove old contributor notes
* [```ce36f02```](https://github.com/ktmeaton/NCBImeta/commit/ce36f02) remove poetry toml
* [```db9889b```](https://github.com/ktmeaton/NCBImeta/commit/db9889b) update Lat and Lon for BioSample schema yaml
* [```e325985```](https://github.com/ktmeaton/NCBImeta/commit/e325985) notes on replacing old schema
* [```ebaf5eb```](https://github.com/ktmeaton/NCBImeta/commit/ebaf5eb) dev testing toml file
* [```e83d45c```](https://github.com/ktmeaton/NCBImeta/commit/e83d45c) formatting notes
* [```4c9de0a```](https://github.com/ktmeaton/NCBImeta/commit/4c9de0a) fix line wrap
* [```ff0bb82```](https://github.com/ktmeaton/NCBImeta/commit/ff0bb82) add python file trigger
* [```cb07360```](https://github.com/ktmeaton/NCBImeta/commit/cb07360) black and flake8 format
* [```1b1454f```](https://github.com/ktmeaton/NCBImeta/commit/1b1454f) poetry experiment
* [```3f1e911```](https://github.com/ktmeaton/NCBImeta/commit/3f1e911) force python 3.7
* [```dc9da60```](https://github.com/ktmeaton/NCBImeta/commit/dc9da60) rename dev requirements
* [```323627e```](https://github.com/ktmeaton/NCBImeta/commit/323627e) basic deploy workflow for pypi
* [```b9d387a```](https://github.com/ktmeaton/NCBImeta/commit/b9d387a) fix indendation
* [```506d2c3```](https://github.com/ktmeaton/NCBImeta/commit/506d2c3) create deploy.yaml
* [```22652e0```](https://github.com/ktmeaton/NCBImeta/commit/22652e0) line endings fix
* [```b072928```](https://github.com/ktmeaton/NCBImeta/commit/b072928) developer dependencies install guide
* [```a35875f```](https://github.com/ktmeaton/NCBImeta/commit/a35875f) fix mixed line endings
* [```5b32ea9```](https://github.com/ktmeaton/NCBImeta/commit/5b32ea9) pre-commit configuration
* [```52e9e79```](https://github.com/ktmeaton/NCBImeta/commit/52e9e79) remove trailing whitespace and file endings
* [```094ff71```](https://github.com/ktmeaton/NCBImeta/commit/094ff71) gh actions and linting
* [```1cc59ea```](https://github.com/ktmeaton/NCBImeta/commit/1cc59ea) fix indentation
* [```5bdb96e```](https://github.com/ktmeaton/NCBImeta/commit/5bdb96e) correct lint extension
* [```0dc356f```](https://github.com/ktmeaton/NCBImeta/commit/0dc356f) restrict paths and add more workflow
* [```6a5590c```](https://github.com/ktmeaton/NCBImeta/commit/6a5590c) remove codecov command
* [```58d10a9```](https://github.com/ktmeaton/NCBImeta/commit/58d10a9) update all names
* [```e648561```](https://github.com/ktmeaton/NCBImeta/commit/e648561) path in quotes
* [```8abaea0```](https://github.com/ktmeaton/NCBImeta/commit/8abaea0) linting yaml with markdown
* [```d1a446a```](https://github.com/ktmeaton/NCBImeta/commit/d1a446a) add workflow to path
* [```fce336d```](https://github.com/ktmeaton/NCBImeta/commit/fce336d) wildcard paths in quotes
* [```e5332c8```](https://github.com/ktmeaton/NCBImeta/commit/e5332c8) all branches restrict paths
* [```ff2e1c6```](https://github.com/ktmeaton/NCBImeta/commit/ff2e1c6) restrict paths remove example
* [```76d1a77```](https://github.com/ktmeaton/NCBImeta/commit/76d1a77) change repo name to action
* [```ba74603```](https://github.com/ktmeaton/NCBImeta/commit/ba74603) remove version number
* [```c2d46a1```](https://github.com/ktmeaton/NCBImeta/commit/c2d46a1) change repo tarball ver
* [```906fc27```](https://github.com/ktmeaton/NCBImeta/commit/906fc27) quick start example
* [```e5e516c```](https://github.com/ktmeaton/NCBImeta/commit/e5e516c) codecov upload
* [```01b938f```](https://github.com/ktmeaton/NCBImeta/commit/01b938f) fix indentation
* [```c65775f```](https://github.com/ktmeaton/NCBImeta/commit/c65775f) ensure unique names
* [```e49c4b1```](https://github.com/ktmeaton/NCBImeta/commit/e49c4b1) formatting and new test
* [```49a8008```](https://github.com/ktmeaton/NCBImeta/commit/49a8008) rename to build
* [```07f4bc7```](https://github.com/ktmeaton/NCBImeta/commit/07f4bc7) add main installation
* [```744a0d4```](https://github.com/ktmeaton/NCBImeta/commit/744a0d4) use the pip tutorial gh actions
* [```d6c38cd```](https://github.com/ktmeaton/NCBImeta/commit/d6c38cd) sudo for pip installs
* [```b31a524```](https://github.com/ktmeaton/NCBImeta/commit/b31a524) try sudo instead
* [```3bec231```](https://github.com/ktmeaton/NCBImeta/commit/3bec231) remember to activate
* [```d37982a```](https://github.com/ktmeaton/NCBImeta/commit/d37982a) try in conda env
* [```eb7ea00```](https://github.com/ktmeaton/NCBImeta/commit/eb7ea00) Try github actions
* [```16799dc```](https://github.com/ktmeaton/NCBImeta/commit/16799dc) test out github actions
* [```3cc1572```](https://github.com/ktmeaton/NCBImeta/commit/3cc1572) fix line ending
* [```0e90105```](https://github.com/ktmeaton/NCBImeta/commit/0e90105) MANIFEST for pypi packaging
* [```228c369```](https://github.com/ktmeaton/NCBImeta/commit/228c369) update v0.6.5 to v0.6.6
* [```7f92720```](https://github.com/ktmeaton/NCBImeta/commit/7f92720) minor update
* [```7740335```](https://github.com/ktmeaton/NCBImeta/commit/7740335) config for manifest
* [```5c4479b```](https://github.com/ktmeaton/NCBImeta/commit/5c4479b) document dev dependencies
* [```42ba5ea```](https://github.com/ktmeaton/NCBImeta/commit/42ba5ea) ignore fusion table
* [```3c9abf3```](https://github.com/ktmeaton/NCBImeta/commit/3c9abf3) Convert BioProjectTitle to XPath query
* [```a150cd6```](https://github.com/ktmeaton/NCBImeta/commit/a150cd6) SQL Update fix for char escape and security
* [```bb52506```](https://github.com/ktmeaton/NCBImeta/commit/bb52506) fix slack PR notification bool
* [```e55544e```](https://github.com/ktmeaton/NCBImeta/commit/e55544e) temp pip instructions
* [```3a88be8```](https://github.com/ktmeaton/NCBImeta/commit/3a88be8) pypi install note

## v0.6.5

### Commits

* [```ccfabe3```](https://github.com/ktmeaton/NCBImeta/commit/ccfabe3) Merge branch 'master' into dev
* [```4531fdf```](https://github.com/ktmeaton/NCBImeta/commit/4531fdf) update dates, ver, community guidelines
* [```a455008```](https://github.com/ktmeaton/NCBImeta/commit/a455008) Move XPath query explanation to README_schema doc
* [```b61f931```](https://github.com/ktmeaton/NCBImeta/commit/b61f931) update version number to v0.6.5
* [```138cc5f```](https://github.com/ktmeaton/NCBImeta/commit/138cc5f) Merge pull request #9 from hellothisisMatt/feature/add-xpath-support
* [```a741506```](https://github.com/ktmeaton/NCBImeta/commit/a741506) Update all other issues template [skip ci]
* [```7d37574```](https://github.com/ktmeaton/NCBImeta/commit/7d37574) Update feature request template
* [```86ccf01```](https://github.com/ktmeaton/NCBImeta/commit/86ccf01) Update bug report templates
* [```66dfb73```](https://github.com/ktmeaton/NCBImeta/commit/66dfb73) remove branch pr9 from travis testing
* [```fd16615```](https://github.com/ktmeaton/NCBImeta/commit/fd16615) update unittests to include new values retrieved with advanced Xpath
* [```564c18f```](https://github.com/ktmeaton/NCBImeta/commit/564c18f) test implementation and error class when Xpath query is empty and unspecified
* [```9b0b11f```](https://github.com/ktmeaton/NCBImeta/commit/9b0b11f) define an error class for when the Xpath query is empty and unspecified
* [```ade807c```](https://github.com/ktmeaton/NCBImeta/commit/ade807c) check if the Xpath query was left empty and unspecified
* [```2c5e9a2```](https://github.com/ktmeaton/NCBImeta/commit/2c5e9a2) add XPath queries to match the example config
* [```35b7b95```](https://github.com/ktmeaton/NCBImeta/commit/35b7b95) Merge pull request #1 from ktmeaton/pr9
* [```659238a```](https://github.com/ktmeaton/NCBImeta/commit/659238a) update changelog with error classes, run CI
* [```c8b6584```](https://github.com/ktmeaton/NCBImeta/commit/c8b6584) New error classes that can be raised in adv_xml_search
* [```c04e981```](https://github.com/ktmeaton/NCBImeta/commit/c04e981) revise adv_xml_search to check type of result
* [```1817705```](https://github.com/ktmeaton/NCBImeta/commit/1817705) remove outdated test xpath commands
* [```a5c2fff```](https://github.com/ktmeaton/NCBImeta/commit/a5c2fff) fix non-specific BioSampleBioProjectAccession
* [```c390ff6```](https://github.com/ktmeaton/NCBImeta/commit/c390ff6) refine xpath search for attr value
* [```8418002```](https://github.com/ktmeaton/NCBImeta/commit/8418002) Enable attribute returns and error checking
* [```a9b8057```](https://github.com/ktmeaton/NCBImeta/commit/a9b8057) enable dev and py9 for CI
* [```ddc6a06```](https://github.com/ktmeaton/NCBImeta/commit/ddc6a06) correct NucleotideBioSampleAccession xpath query
* [```8f0e2ea```](https://github.com/ktmeaton/NCBImeta/commit/8f0e2ea) simplified test example for tip search
* [```0a5fb69```](https://github.com/ktmeaton/NCBImeta/commit/0a5fb69) add contrib info
* [```bffe829```](https://github.com/ktmeaton/NCBImeta/commit/bffe829) XPATH use, remove nonspecific NucleotideAssemblyAccession
* [```16451e7```](https://github.com/ktmeaton/NCBImeta/commit/16451e7) clarify comma separated with space
* [```6e972c2```](https://github.com/ktmeaton/NCBImeta/commit/6e972c2) 2 tests of adv_xml_search starting from root or starting from tip
* [```5158486```](https://github.com/ktmeaton/NCBImeta/commit/5158486) xpath param of adv_xml_search as preformatted xml
* [```c22c859```](https://github.com/ktmeaton/NCBImeta/commit/c22c859) document the space and comma
* [```65353dc```](https://github.com/ktmeaton/NCBImeta/commit/65353dc) Merge upstream dev changes into pr9
* [```5585964```](https://github.com/ktmeaton/NCBImeta/commit/5585964) Add XPATH Information to config readme
* [```56e0b3c```](https://github.com/ktmeaton/NCBImeta/commit/56e0b3c) Add ability to use full XPath for XML searching
* [```3206959```](https://github.com/ktmeaton/NCBImeta/commit/3206959) new to do fixes
* [```62ab0e7```](https://github.com/ktmeaton/NCBImeta/commit/62ab0e7) SRABioProjectAccession edge cases
* [```0cadcec```](https://github.com/ktmeaton/NCBImeta/commit/0cadcec) SRABioSampleAccession typo
* [```5232f35```](https://github.com/ktmeaton/NCBImeta/commit/5232f35) update SRABioProjectAccession
* [```24bdfe0```](https://github.com/ktmeaton/NCBImeta/commit/24bdfe0) v0.6.5 init changes
* [```cd5efb4```](https://github.com/ktmeaton/NCBImeta/commit/cd5efb4) update BioSample schema and example for SRABioSampleAccession
* [```fb57950```](https://github.com/ktmeaton/NCBImeta/commit/fb57950) badge update and rearrange [skip ci]

## v0.6.4

### Commits

* [```33cdca3```](https://github.com/ktmeaton/NCBImeta/commit/33cdca3) switch travis back to master only
* [```03481fe```](https://github.com/ktmeaton/NCBImeta/commit/03481fe) expand issues section [skip ci]
* [```1b94dfe```](https://github.com/ktmeaton/NCBImeta/commit/1b94dfe) upcoming features section [skip ci]
* [```a55df18```](https://github.com/ktmeaton/NCBImeta/commit/a55df18) doc update [skip ci]
* [```8088943```](https://github.com/ktmeaton/NCBImeta/commit/8088943) test multi match
* [```7e5cabd```](https://github.com/ktmeaton/NCBImeta/commit/7e5cabd) update nucleotide data
* [```575b153```](https://github.com/ktmeaton/NCBImeta/commit/575b153) switch travis to dev branch
* [```9e56cec```](https://github.com/ktmeaton/NCBImeta/commit/9e56cec) small progress on nucleotide data
* [```a81357e```](https://github.com/ktmeaton/NCBImeta/commit/a81357e) v0.6.4 notes
* [```c7dd8fa```](https://github.com/ktmeaton/NCBImeta/commit/c7dd8fa) ignore some example and test files
* [```d3fa6cc```](https://github.com/ktmeaton/NCBImeta/commit/d3fa6cc) ver update
* [```a4a3d9e```](https://github.com/ktmeaton/NCBImeta/commit/a4a3d9e) commit to encode removal
* [```d3abba5```](https://github.com/ktmeaton/NCBImeta/commit/d3abba5) working on nucleotide export
* [```300055a```](https://github.com/ktmeaton/NCBImeta/commit/300055a) encode decode purge
* [```cdf506f```](https://github.com/ktmeaton/NCBImeta/commit/cdf506f) uncomment encode decode calls
* [```33192fa```](https://github.com/ktmeaton/NCBImeta/commit/33192fa) SQL parameter and table col check
* [```31f9203```](https://github.com/ktmeaton/NCBImeta/commit/31f9203) Extra param checking Issue #7
* [```585a200```](https://github.com/ktmeaton/NCBImeta/commit/585a200) table and column name catching
* [```5b46aaa```](https://github.com/ktmeaton/NCBImeta/commit/5b46aaa) SQL parameter complete
* [```b11f401```](https://github.com/ktmeaton/NCBImeta/commit/b11f401) sanitize output format again
* [```4eca45d```](https://github.com/ktmeaton/NCBImeta/commit/4eca45d) sanitize output format
* [```69b1f76```](https://github.com/ktmeaton/NCBImeta/commit/69b1f76) test for table and column name format
* [```8104b94```](https://github.com/ktmeaton/NCBImeta/commit/8104b94) SQL table name add
* [```73af61a```](https://github.com/ktmeaton/NCBImeta/commit/73af61a) before more sql select change
* [```e4c17c7```](https://github.com/ktmeaton/NCBImeta/commit/e4c17c7) slack travis-ci integration

## v0.6.3

### Commits

* [```77e6fe1```](https://github.com/ktmeaton/NCBImeta/commit/77e6fe1) enable master branch CI
* [```8850b76```](https://github.com/ktmeaton/NCBImeta/commit/8850b76) codecov spacing syntax error
* [```523687b```](https://github.com/ktmeaton/NCBImeta/commit/523687b) allow dev CI
* [```16a4f56```](https://github.com/ktmeaton/NCBImeta/commit/16a4f56) bioconda autobump bot rely
* [```4446f94```](https://github.com/ktmeaton/NCBImeta/commit/4446f94) remove bioconda scripts
* [```b92681f```](https://github.com/ktmeaton/NCBImeta/commit/b92681f) autobump script update
* [```49e4bab```](https://github.com/ktmeaton/NCBImeta/commit/49e4bab) hidden script names and autobump script [skip ci]
* [```502338b```](https://github.com/ktmeaton/NCBImeta/commit/502338b) resume full testing
* [```a76ff01```](https://github.com/ktmeaton/NCBImeta/commit/a76ff01) use the bioconda recommended install code [skip ci]
* [```052a82c```](https://github.com/ktmeaton/NCBImeta/commit/052a82c) solved mystery of branch bc-tbd [skip ci]
* [```324b668```](https://github.com/ktmeaton/NCBImeta/commit/324b668) automate ver number update [skip ci]
* [```41ed42e```](https://github.com/ktmeaton/NCBImeta/commit/41ed42e) switch zenodo citation to the all ver DOI [skip ci]
* [```c00f478```](https://github.com/ktmeaton/NCBImeta/commit/c00f478) remove unnecessary tag condition
* [```2d5dda8```](https://github.com/ktmeaton/NCBImeta/commit/2d5dda8) if statement spacing
* [```f139b22```](https://github.com/ktmeaton/NCBImeta/commit/f139b22) simplify commit message [skip ci]
* [```d6c3c02```](https://github.com/ktmeaton/NCBImeta/commit/d6c3c02) travis tag match v
* [```553124d```](https://github.com/ktmeaton/NCBImeta/commit/553124d) tag check
* [```fb349ca```](https://github.com/ktmeaton/NCBImeta/commit/fb349ca) file format linux
* [```086851a```](https://github.com/ktmeaton/NCBImeta/commit/086851a) executable mode
* [```4053814```](https://github.com/ktmeaton/NCBImeta/commit/4053814) move conda update commands to script
* [```2329288```](https://github.com/ktmeaton/NCBImeta/commit/2329288) set branch upstream origin
* [```6b3d224```](https://github.com/ktmeaton/NCBImeta/commit/6b3d224) add, commit, push recipe
* [```20b58f6```](https://github.com/ktmeaton/NCBImeta/commit/20b58f6) enforce semicolon
* [```98609ac```](https://github.com/ktmeaton/NCBImeta/commit/98609ac) split up branch creation and check out
* [```c9ac4a5```](https://github.com/ktmeaton/NCBImeta/commit/c9ac4a5) no tag enforcement yet
* [```03b675e```](https://github.com/ktmeaton/NCBImeta/commit/03b675e) forgotten symbols
* [```5d742c2```](https://github.com/ktmeaton/NCBImeta/commit/5d742c2) sed replacement check
* [```d6268fe```](https://github.com/ktmeaton/NCBImeta/commit/d6268fe) src sha256 refine
* [```4c11045```](https://github.com/ktmeaton/NCBImeta/commit/4c11045) path and url fix
* [```d8fc06b```](https://github.com/ktmeaton/NCBImeta/commit/d8fc06b) conda var first pass
* [```5f6d8ca```](https://github.com/ktmeaton/NCBImeta/commit/5f6d8ca) hide user ref [skip ci]
* [```e48c1e6```](https://github.com/ktmeaton/NCBImeta/commit/e48c1e6) remove dir checks
* [```87266f3```](https://github.com/ktmeaton/NCBImeta/commit/87266f3) remove comments
* [```4e1606b```](https://github.com/ktmeaton/NCBImeta/commit/4e1606b) travis-ci test post local verify
* [```4554da9```](https://github.com/ktmeaton/NCBImeta/commit/4554da9) fix incorrect slug usage
* [```a0fa57e```](https://github.com/ktmeaton/NCBImeta/commit/a0fa57e) show git config
* [```664a2aa```](https://github.com/ktmeaton/NCBImeta/commit/664a2aa) directory checking
* [```3d711f3```](https://github.com/ktmeaton/NCBImeta/commit/3d711f3) config param fix
* [```a1a9fa0```](https://github.com/ktmeaton/NCBImeta/commit/a1a9fa0) bioconda update attempt 2
* [```dfbd7e3```](https://github.com/ktmeaton/NCBImeta/commit/dfbd7e3) bioconda update attempt
* [```8f0c064```](https://github.com/ktmeaton/NCBImeta/commit/8f0c064) more pwd checks
* [```95801e5```](https://github.com/ktmeaton/NCBImeta/commit/95801e5) pwd test
* [```b6b0856```](https://github.com/ktmeaton/NCBImeta/commit/b6b0856) condition testing
* [```cb49712```](https://github.com/ktmeaton/NCBImeta/commit/cb49712) Document the bioconda repo master branch switch

## v0.6.2

### Commits

* [```3f8fac5```](https://github.com/ktmeaton/NCBImeta/commit/3f8fac5) PR and Issue linking
* [```4c94207```](https://github.com/ktmeaton/NCBImeta/commit/4c94207) branch match ver release name [skip ci]
* [```d421128```](https://github.com/ktmeaton/NCBImeta/commit/d421128) date update v0.6.2 [skip ci]
* [```7dc34dd```](https://github.com/ktmeaton/NCBImeta/commit/7dc34dd) automate pypi deploy [skip ci]
* [```e636c6a```](https://github.com/ktmeaton/NCBImeta/commit/e636c6a) CI branch on master
* [```37aedb5```](https://github.com/ktmeaton/NCBImeta/commit/37aedb5) Merge branch 'dev'
* [```71ab0e8```](https://github.com/ktmeaton/NCBImeta/commit/71ab0e8) merge dev into master
* [```3de43a8```](https://github.com/ktmeaton/NCBImeta/commit/3de43a8) Merge branch 'dev' of https://github.com/ktmeaton/NCBImeta into dev
* [```7a2a3f3```](https://github.com/ktmeaton/NCBImeta/commit/7a2a3f3) version update
* [```a9605de```](https://github.com/ktmeaton/NCBImeta/commit/a9605de) version update [skip ci]
* [```2dc5920```](https://github.com/ktmeaton/NCBImeta/commit/2dc5920) ver update to 0.6.2 [skip ci]
* [```92d810d```](https://github.com/ktmeaton/NCBImeta/commit/92d810d) Ouput dir is created instead of raising error [skip ci]
* [```bd10fbd```](https://github.com/ktmeaton/NCBImeta/commit/bd10fbd) remove slash [skip ci]
* [```f820068```](https://github.com/ktmeaton/NCBImeta/commit/f820068) bioconda installation and media as raw links [skip ci]
* [```f994d97```](https://github.com/ktmeaton/NCBImeta/commit/f994d97) Merge branch 'dev' of https://github.com/ktmeaton/NCBImeta into dev
* [```f89ccc5```](https://github.com/ktmeaton/NCBImeta/commit/f89ccc5) CI on dev
* [```5be403b```](https://github.com/ktmeaton/NCBImeta/commit/5be403b) v0.6.2 bioconda and output dir
* [```619e1db```](https://github.com/ktmeaton/NCBImeta/commit/619e1db) v0.6.2 bioconda and output dir
* [```c640928```](https://github.com/ktmeaton/NCBImeta/commit/c640928) Remove the output dir error class [skip ci]
* [```a68d17a```](https://github.com/ktmeaton/NCBImeta/commit/a68d17a) Make output dir instead of raising error
* [```a99a4bf```](https://github.com/ktmeaton/NCBImeta/commit/a99a4bf) Test if output dir is created [skip ci]
* [```c0c1019```](https://github.com/ktmeaton/NCBImeta/commit/c0c1019) Merge pull request #5 from druvus/patch-1
* [```c3fc9a4```](https://github.com/ktmeaton/NCBImeta/commit/c3fc9a4) Paper remove unnecessary reference ncbi website [skip ci]
* [```55e8bf7```](https://github.com/ktmeaton/NCBImeta/commit/55e8bf7) Adding Bioconda as an installation option
* [```17d6ba1```](https://github.com/ktmeaton/NCBImeta/commit/17d6ba1) Zenodo citation and badge [skip ci]

## v0.6.1

### Commits

* [```0ab131a```](https://github.com/ktmeaton/NCBImeta/commit/0ab131a) update schema doc for new xml parsing [skip ci]
* [```680ddb2```](https://github.com/ktmeaton/NCBImeta/commit/680ddb2) ref and capitalization [skip ci]
* [```edaa3b0```](https://github.com/ktmeaton/NCBImeta/commit/edaa3b0) python 3 spacing [skip ci]
* [```49132f9```](https://github.com/ktmeaton/NCBImeta/commit/49132f9) re-enable master CI
* [```fb3a01c```](https://github.com/ktmeaton/NCBImeta/commit/fb3a01c) aeruginosa text db
* [```7dd9804```](https://github.com/ktmeaton/NCBImeta/commit/7dd9804) aeruginosa paper db
* [```c50677f```](https://github.com/ktmeaton/NCBImeta/commit/c50677f) add release names
* [```a93d1a1```](https://github.com/ktmeaton/NCBImeta/commit/a93d1a1) Name for upcoming v0.6.1
* [```a24a82f```](https://github.com/ktmeaton/NCBImeta/commit/a24a82f) release to development
* [```32d162b```](https://github.com/ktmeaton/NCBImeta/commit/32d162b) prematurely add v0.6.1 links
* [```88d6d0c```](https://github.com/ktmeaton/NCBImeta/commit/88d6d0c) format and typos
* [```5ccd379```](https://github.com/ktmeaton/NCBImeta/commit/5ccd379) conceptual rearrange
* [```be94004```](https://github.com/ktmeaton/NCBImeta/commit/be94004) Initial commit desc
* [```37068fc```](https://github.com/ktmeaton/NCBImeta/commit/37068fc) complicated compare
* [```8637310```](https://github.com/ktmeaton/NCBImeta/commit/8637310) link fix
* [```900c6a5```](https://github.com/ktmeaton/NCBImeta/commit/900c6a5) test commit compare
* [```8efd081```](https://github.com/ktmeaton/NCBImeta/commit/8efd081) hyperlink exp
* [```3b17a51```](https://github.com/ktmeaton/NCBImeta/commit/3b17a51) aeruginosa config file
* [```fdc4672```](https://github.com/ktmeaton/NCBImeta/commit/fdc4672) list format tables
* [```4bbeb10```](https://github.com/ktmeaton/NCBImeta/commit/4bbeb10) rename
* [```f7f5311```](https://github.com/ktmeaton/NCBImeta/commit/f7f5311) remove old gif
* [```27d1a2d```](https://github.com/ktmeaton/NCBImeta/commit/27d1a2d) biosample gif
* [```e5716aa```](https://github.com/ktmeaton/NCBImeta/commit/e5716aa) gif larger
* [```79b82f1```](https://github.com/ktmeaton/NCBImeta/commit/79b82f1) gif rename
* [```867775a```](https://github.com/ktmeaton/NCBImeta/commit/867775a) image swap
* [```37dd9c9```](https://github.com/ktmeaton/NCBImeta/commit/37dd9c9) gif full path, reformat print
* [```107fa14```](https://github.com/ktmeaton/NCBImeta/commit/107fa14) State CLI Application
* [```69f1765```](https://github.com/ktmeaton/NCBImeta/commit/69f1765) specify CLI application
* [```8a72fc8```](https://github.com/ktmeaton/NCBImeta/commit/8a72fc8) asciicast
* [```7d0a4cc```](https://github.com/ktmeaton/NCBImeta/commit/7d0a4cc) dependency reference to file
* [```e2674d7```](https://github.com/ktmeaton/NCBImeta/commit/e2674d7) JOSS file additions
* [```bec9abe```](https://github.com/ktmeaton/NCBImeta/commit/bec9abe) update record number
* [```08b5e58```](https://github.com/ktmeaton/NCBImeta/commit/08b5e58) uncomment xml print
* [```084f200```](https://github.com/ktmeaton/NCBImeta/commit/084f200) jpg rename
* [```d1bd7ae```](https://github.com/ktmeaton/NCBImeta/commit/d1bd7ae) doc xml printing, disable master ci
* [```bd92205```](https://github.com/ktmeaton/NCBImeta/commit/bd92205) doc update API optional [skip ci]
* [```2c72df8```](https://github.com/ktmeaton/NCBImeta/commit/2c72df8) doc update node parse [skip ci]
* [```c7a77e4```](https://github.com/ktmeaton/NCBImeta/commit/c7a77e4) new prettyprint func [skip ci]
* [```c78b8c2```](https://github.com/ktmeaton/NCBImeta/commit/c78b8c2) update version [skip ci]
* [```ec6ebf8```](https://github.com/ktmeaton/NCBImeta/commit/ec6ebf8) edit description [skip ci]
* [```daa228b```](https://github.com/ktmeaton/NCBImeta/commit/daa228b) JOSS Paper
* [```82b7b4f```](https://github.com/ktmeaton/NCBImeta/commit/82b7b4f) Reformat headings [skip ci]
* [```7cbabc9```](https://github.com/ktmeaton/NCBImeta/commit/7cbabc9) config file name update [skip ci]
* [```9808a6b```](https://github.com/ktmeaton/NCBImeta/commit/9808a6b) Clarify and formatting [skip ci]

## v0.6.0

### Commits

* [```44805c7```](https://github.com/ktmeaton/NCBImeta/commit/44805c7) remove exp code
* [```717431a```](https://github.com/ktmeaton/NCBImeta/commit/717431a) license badge update [skip ci]
* [```05b3aea```](https://github.com/ktmeaton/NCBImeta/commit/05b3aea) Merge branch 'dev'
* [```cb575ae```](https://github.com/ktmeaton/NCBImeta/commit/cb575ae) conflict resolve
* [```b489ca6```](https://github.com/ktmeaton/NCBImeta/commit/b489ca6) v0.6.0 updates [skip ci]
* [```778ab45```](https://github.com/ktmeaton/NCBImeta/commit/778ab45) ignore test database and log files [skip ci]
* [```b3a310c```](https://github.com/ktmeaton/NCBImeta/commit/b3a310c) Single quote sql query fix
* [```0e0be9c```](https://github.com/ktmeaton/NCBImeta/commit/0e0be9c) XML overhaul test update
* [```b0dc3bc```](https://github.com/ktmeaton/NCBImeta/commit/b0dc3bc) successful run
* [```848cb53```](https://github.com/ktmeaton/NCBImeta/commit/848cb53) quotation fix
* [```e03597a```](https://github.com/ktmeaton/NCBImeta/commit/e03597a) bugfixes for xml_search
* [```e0ff551```](https://github.com/ktmeaton/NCBImeta/commit/e0ff551) remove old search and flatten functions
* [```546c368```](https://github.com/ktmeaton/NCBImeta/commit/546c368) Before function move
* [```9d97bbe```](https://github.com/ktmeaton/NCBImeta/commit/9d97bbe) efetch part functional for sra
* [```ab88a2c```](https://github.com/ktmeaton/NCBImeta/commit/ab88a2c) before nucleotide switch to efetch
* [```393466d```](https://github.com/ktmeaton/NCBImeta/commit/393466d) bioproject with efetch
* [```be00272```](https://github.com/ktmeaton/NCBImeta/commit/be00272) efetch switch
* [```ee6cf17```](https://github.com/ktmeaton/NCBImeta/commit/ee6cf17) lxml overhaul
* [```31de27e```](https://github.com/ktmeaton/NCBImeta/commit/31de27e) before efetch switch
* [```1b720a8```](https://github.com/ktmeaton/NCBImeta/commit/1b720a8) xml in xml parsing
* [```d33d58e```](https://github.com/ktmeaton/NCBImeta/commit/d33d58e) cdata testing
* [```940ba36```](https://github.com/ktmeaton/NCBImeta/commit/940ba36) lxml experiment
* [```88b7c0d```](https://github.com/ktmeaton/NCBImeta/commit/88b7c0d) Config and Schema README as full path [skip ci]
* [```0b1030d```](https://github.com/ktmeaton/NCBImeta/commit/0b1030d) ignore build dir [skip ci]
* [```e4be83d```](https://github.com/ktmeaton/NCBImeta/commit/e4be83d) change codecov to branch master [skip ci]
* [```f3b30de```](https://github.com/ktmeaton/NCBImeta/commit/f3b30de) v0.5.0 version update [skip ci]

## v0.5.0

### Commits

* [```046f68b```](https://github.com/ktmeaton/NCBImeta/commit/046f68b) travis-ci only master
* [```cecb1aa```](https://github.com/ktmeaton/NCBImeta/commit/cecb1aa) Merge branch 'dev'
* [```0d84457```](https://github.com/ktmeaton/NCBImeta/commit/0d84457) v0.5.0 updates [skip ci]
* [```849f733```](https://github.com/ktmeaton/NCBImeta/commit/849f733) Explain default slow download [skip ci]
* [```b810c2a```](https://github.com/ktmeaton/NCBImeta/commit/b810c2a) Fix codecov badge link [skip ci]
* [```5d4c83b```](https://github.com/ktmeaton/NCBImeta/commit/5d4c83b) Switch release badge to PyPI [skip ci]
* [```4112027```](https://github.com/ktmeaton/NCBImeta/commit/4112027) codecov url fix [skip ci]
* [```29a80f5```](https://github.com/ktmeaton/NCBImeta/commit/29a80f5) try codecov badge with dev [skip ci]
* [```1488da7```](https://github.com/ktmeaton/NCBImeta/commit/1488da7) Clear debugging output
* [```199559b```](https://github.com/ktmeaton/NCBImeta/commit/199559b) formatting [skip ci]
* [```e6bd3b7```](https://github.com/ktmeaton/NCBImeta/commit/e6bd3b7) document column_index bugfix [skip ci]
* [```f560dfe```](https://github.com/ktmeaton/NCBImeta/commit/f560dfe) Correct bioproject accession
* [```10e1d2a```](https://github.com/ktmeaton/NCBImeta/commit/10e1d2a) column_index reposition
* [```75edf10```](https://github.com/ktmeaton/NCBImeta/commit/75edf10) http error catching for esearch
* [```1bc5317```](https://github.com/ktmeaton/NCBImeta/commit/1bc5317) strpath for tmpdir [skip ci]
* [```fd509f4```](https://github.com/ktmeaton/NCBImeta/commit/fd509f4) Fix sys path and list comparison
* [```3eb8367```](https://github.com/ktmeaton/NCBImeta/commit/3eb8367) Remove Python 3.4 linux build [skip ci]
* [```6bd4d8d```](https://github.com/ktmeaton/NCBImeta/commit/6bd4d8d) require time module
* [```8df4245```](https://github.com/ktmeaton/NCBImeta/commit/8df4245) test db path fix [skip ci]
* [```1488745```](https://github.com/ktmeaton/NCBImeta/commit/1488745) Python 3.5+ required (no more 3.4) [skip ci]
* [```d5624e1```](https://github.com/ktmeaton/NCBImeta/commit/d5624e1) travis troubleshooting
* [```763c54f```](https://github.com/ktmeaton/NCBImeta/commit/763c54f) bash uploader for codecov
* [```b807bb0```](https://github.com/ktmeaton/NCBImeta/commit/b807bb0) example config reset
* [```595959b```](https://github.com/ktmeaton/NCBImeta/commit/595959b) Merge branch 'dev' of https://github.com/ktmeaton/NCBImeta into dev
* [```4c9b15b```](https://github.com/ktmeaton/NCBImeta/commit/4c9b15b) full travis-ci test build
* [```8f5f38b```](https://github.com/ktmeaton/NCBImeta/commit/8f5f38b) full travis-ci test build
* [```e2ede2c```](https://github.com/ktmeaton/NCBImeta/commit/e2ede2c) pytest for non-flat mode
* [```7d4ec41```](https://github.com/ktmeaton/NCBImeta/commit/7d4ec41) prototype master table check
* [```ea642b8```](https://github.com/ktmeaton/NCBImeta/commit/ea642b8) Move the HTTPErrorCatch method to utilities [skip ci]
* [```dd55698```](https://github.com/ktmeaton/NCBImeta/commit/dd55698) pytest db verify SRA
* [```8de4feb```](https://github.com/ktmeaton/NCBImeta/commit/8de4feb) pytesting and db verification [skip ci]
* [```bc7f6b2```](https://github.com/ktmeaton/NCBImeta/commit/bc7f6b2) Update test biosample metadata
* [```6771d87```](https://github.com/ktmeaton/NCBImeta/commit/6771d87) Exclude temporary test files [skip ci]
* [```b2df525```](https://github.com/ktmeaton/NCBImeta/commit/b2df525) pytest pubmed [skip ci]
* [```a88ba92```](https://github.com/ktmeaton/NCBImeta/commit/a88ba92) pytest export nucleotide table values
* [```a21df9e```](https://github.com/ktmeaton/NCBImeta/commit/a21df9e) pytest conftest assembly and bioproject
* [```fe28e8b```](https://github.com/ktmeaton/NCBImeta/commit/fe28e8b) pytest export assembly and bioproject
* [```918f1e3```](https://github.com/ktmeaton/NCBImeta/commit/918f1e3) remove annot file [skip ci]
* [```26b8d15```](https://github.com/ktmeaton/NCBImeta/commit/26b8d15) pytest join
* [```87ca364```](https://github.com/ktmeaton/NCBImeta/commit/87ca364) Comment clarification [skip ci]
* [```e84af09```](https://github.com/ktmeaton/NCBImeta/commit/e84af09) pytest annotatereplace
* [```90da750```](https://github.com/ktmeaton/NCBImeta/commit/90da750) annotateconcat 90% cov
* [```34725c9```](https://github.com/ktmeaton/NCBImeta/commit/34725c9) fix indentation
* [```4eea6b3```](https://github.com/ktmeaton/NCBImeta/commit/4eea6b3) annotation file for pytest [skip ci]
* [```1a4f068```](https://github.com/ktmeaton/NCBImeta/commit/1a4f068) specify pytest files in order
* [```98a4221```](https://github.com/ktmeaton/NCBImeta/commit/98a4221) module import in test dir [skip ci]
* [```8827992```](https://github.com/ktmeaton/NCBImeta/commit/8827992) sra table metadata re-fix [skip ci]
* [```debd22a```](https://github.com/ktmeaton/NCBImeta/commit/debd22a) debug help points [skip ci]
* [```7150dc3```](https://github.com/ktmeaton/NCBImeta/commit/7150dc3) generic test db name [skip ci]
* [```cac8a12```](https://github.com/ktmeaton/NCBImeta/commit/cac8a12) pytest annotateconcatenate [skip ci]
* [```e16d6e8```](https://github.com/ktmeaton/NCBImeta/commit/e16d6e8) formatting [skip ci]
* [```469bb6a```](https://github.com/ktmeaton/NCBImeta/commit/469bb6a) proper support for SRA [skip ci]
* [```b5e5753```](https://github.com/ktmeaton/NCBImeta/commit/b5e5753) codecov badge [skip ci]
* [```2622124```](https://github.com/ktmeaton/NCBImeta/commit/2622124) remove extra script
* [```081a76e```](https://github.com/ktmeaton/NCBImeta/commit/081a76e) execute permissions [skip ci]
* [```ef5b084```](https://github.com/ktmeaton/NCBImeta/commit/ef5b084) remove the test1 ref
* [```14133ba```](https://github.com/ktmeaton/NCBImeta/commit/14133ba) update execute permissions again
* [```90d9bfb```](https://github.com/ktmeaton/NCBImeta/commit/90d9bfb) test for output dir
* [```9a1533e```](https://github.com/ktmeaton/NCBImeta/commit/9a1533e) config data yaml test [skip ci]
* [```12b6736```](https://github.com/ktmeaton/NCBImeta/commit/12b6736) old cov [skip ci]
* [```14ce1bf```](https://github.com/ktmeaton/NCBImeta/commit/14ce1bf) remove extraneous [skip ci]
* [```7d8037b```](https://github.com/ktmeaton/NCBImeta/commit/7d8037b) simplify name [skip ci]
* [```f7b4502```](https://github.com/ktmeaton/NCBImeta/commit/f7b4502) yaml file error test [skip ci]
* [```727e547```](https://github.com/ktmeaton/NCBImeta/commit/727e547) main test [skip ci]
* [```1efd279```](https://github.com/ktmeaton/NCBImeta/commit/1efd279) give unique names [skip ci]
* [```dac651b```](https://github.com/ktmeaton/NCBImeta/commit/dac651b) print config file path [skip ci]
* [```e2de196```](https://github.com/ktmeaton/NCBImeta/commit/e2de196) yaml file test [skip ci]
* [```0eb1681```](https://github.com/ktmeaton/NCBImeta/commit/0eb1681) yaml files for testing [skip ci]
* [```805efa6```](https://github.com/ktmeaton/NCBImeta/commit/805efa6) super simple config [skip ci]
* [```de49481```](https://github.com/ktmeaton/NCBImeta/commit/de49481) Proper os join
* [```7a9be44```](https://github.com/ktmeaton/NCBImeta/commit/7a9be44) pytest Errors complete
* [```8afd1ce```](https://github.com/ktmeaton/NCBImeta/commit/8afd1ce) pytest errors
* [```ebfc8af```](https://github.com/ktmeaton/NCBImeta/commit/ebfc8af) duplicate column debugging [skip ci]
* [```b663a97```](https://github.com/ktmeaton/NCBImeta/commit/b663a97) str repr of errors proper return [skip ci]
* [```c51cd55```](https://github.com/ktmeaton/NCBImeta/commit/c51cd55) check for duplicate column names
* [```1725bc6```](https://github.com/ktmeaton/NCBImeta/commit/1725bc6) pytest and coverage [skip ci]
* [```35ffd83```](https://github.com/ktmeaton/NCBImeta/commit/35ffd83) cleanup docstring [skip ci]
* [```640ed49```](https://github.com/ktmeaton/NCBImeta/commit/640ed49) use tmp dir and files
* [```13b0edb```](https://github.com/ktmeaton/NCBImeta/commit/13b0edb) remove unnecessary code [skip ci]
* [```fb3bff1```](https://github.com/ktmeaton/NCBImeta/commit/fb3bff1) Re-enable full travis-ci [skip ci]
* [```ac82f6b```](https://github.com/ktmeaton/NCBImeta/commit/ac82f6b) remove unnecessary type checking
* [```1a801bc```](https://github.com/ktmeaton/NCBImeta/commit/1a801bc) codecov take 2
* [```e6b8977```](https://github.com/ktmeaton/NCBImeta/commit/e6b8977) codecov coverage test upload [skip ci]
* [```40a601f```](https://github.com/ktmeaton/NCBImeta/commit/40a601f) remove testing modules to travis requirements
* [```6021380```](https://github.com/ktmeaton/NCBImeta/commit/6021380) codecov require and test
* [```a35e251```](https://github.com/ktmeaton/NCBImeta/commit/a35e251) pytest-cov requirement
* [```ef91447```](https://github.com/ktmeaton/NCBImeta/commit/ef91447) pytest troubleshooting
* [```3485096```](https://github.com/ktmeaton/NCBImeta/commit/3485096) --cov-report= proper parameter
* [```2bbd425```](https://github.com/ktmeaton/NCBImeta/commit/2bbd425) pytest and codecov
* [```ba90d3f```](https://github.com/ktmeaton/NCBImeta/commit/ba90d3f) Longer description [skip ci]
* [```971a633```](https://github.com/ktmeaton/NCBImeta/commit/971a633) Comment out debugging [skip ci]
* [```eb75e65```](https://github.com/ktmeaton/NCBImeta/commit/eb75e65) v0.4.3 changelog updates [skip ci]
* [```3fd7cec```](https://github.com/ktmeaton/NCBImeta/commit/3fd7cec) Remove unicode ref
* [```e908625```](https://github.com/ktmeaton/NCBImeta/commit/e908625) Remove XmlXXXConfig Functions
* [```120cfa5```](https://github.com/ktmeaton/NCBImeta/commit/120cfa5) Extended description, remove troubleshooting printout [skip ci]
* [```a462b23```](https://github.com/ktmeaton/NCBImeta/commit/a462b23) Namespace clarify
* [```6bde5f2```](https://github.com/ktmeaton/NCBImeta/commit/6bde5f2) Expanded type checking
* [```1c7bd73```](https://github.com/ktmeaton/NCBImeta/commit/1c7bd73) StringElement items troubleshooting
* [```62aceba```](https://github.com/ktmeaton/NCBImeta/commit/62aceba) Extra type checking
* [```3042350```](https://github.com/ktmeaton/NCBImeta/commit/3042350) Fix internal recursion function name
* [```5b4de6b```](https://github.com/ktmeaton/NCBImeta/commit/5b4de6b) Correct and document flatten_dict implementation
* [```4bdedf3```](https://github.com/ktmeaton/NCBImeta/commit/4bdedf3) More testing
* [```efb3b86```](https://github.com/ktmeaton/NCBImeta/commit/efb3b86) Utilities testing
* [```269ca82```](https://github.com/ktmeaton/NCBImeta/commit/269ca82) Recode unicode and remove sys
* [```a20daf6```](https://github.com/ktmeaton/NCBImeta/commit/a20daf6) Description in header [skip ci]
* [```05f693e```](https://github.com/ktmeaton/NCBImeta/commit/05f693e) Remove sys module and flushprint [skip ci]
* [```e9efc27```](https://github.com/ktmeaton/NCBImeta/commit/e9efc27) cleanup [skip ci]
* [```d92dc09```](https://github.com/ktmeaton/NCBImeta/commit/d92dc09) NCBImetaErrors namespace
* [```fc98304```](https://github.com/ktmeaton/NCBImeta/commit/fc98304) Remove io module
* [```78032d9```](https://github.com/ktmeaton/NCBImeta/commit/78032d9) Remove sys module [skip ci]
* [```f76c6dc```](https://github.com/ktmeaton/NCBImeta/commit/f76c6dc) Error classes documented
* [```45ddc1b```](https://github.com/ktmeaton/NCBImeta/commit/45ddc1b) Description in header [skip ci]
* [```1f86ec3```](https://github.com/ktmeaton/NCBImeta/commit/1f86ec3) Remove sys module and flushprint
* [```b1b668a```](https://github.com/ktmeaton/NCBImeta/commit/b1b668a) Remove sys module
* [```8400f5f```](https://github.com/ktmeaton/NCBImeta/commit/8400f5f) Code doc first pass [skip ci]
* [```b681b08```](https://github.com/ktmeaton/NCBImeta/commit/b681b08) Reduce code redundancy
* [```86be602```](https://github.com/ktmeaton/NCBImeta/commit/86be602) Document UpdateDB [skip CI]
* [```d7a22f2```](https://github.com/ktmeaton/NCBImeta/commit/d7a22f2) Force unicode str recode
* [```691111a```](https://github.com/ktmeaton/NCBImeta/commit/691111a) Cleanup, document HTTPErrorCatch
* [```778817b```](https://github.com/ktmeaton/NCBImeta/commit/778817b) import slim down
* [```470e342```](https://github.com/ktmeaton/NCBImeta/commit/470e342) flush print true case sensitive
* [```4addaaa```](https://github.com/ktmeaton/NCBImeta/commit/4addaaa) Travis CI for dev
* [```5ae1c7a```](https://github.com/ktmeaton/NCBImeta/commit/5ae1c7a) Replace flushprint method python3
* [```90bedfb```](https://github.com/ktmeaton/NCBImeta/commit/90bedfb) Merge branch 'master' into dev
* [```5f3548e```](https://github.com/ktmeaton/NCBImeta/commit/5f3548e) Header change
* [```1272291```](https://github.com/ktmeaton/NCBImeta/commit/1272291) Gif reupload and typo fix
* [```359b0ee```](https://github.com/ktmeaton/NCBImeta/commit/359b0ee) Nucleotidet typo fix
* [```5e56fee```](https://github.com/ktmeaton/NCBImeta/commit/5e56fee) gif reupload
* [```6f6d74c```](https://github.com/ktmeaton/NCBImeta/commit/6f6d74c) Delete broken DB gif
* [```4e53ed0```](https://github.com/ktmeaton/NCBImeta/commit/4e53ed0) Delete broken CLI gif
* [```546edd6```](https://github.com/ktmeaton/NCBImeta/commit/546edd6) Move requirements higher [skip ci]

## v0.4.2

### Commits

* [```37adf6e```](https://github.com/ktmeaton/NCBImeta/commit/37adf6e) v0.4.2 finish
* [```296f34f```](https://github.com/ktmeaton/NCBImeta/commit/296f34f) Require cloning github [skip ci]
* [```2847a27```](https://github.com/ktmeaton/NCBImeta/commit/2847a27) section break [skip ci]
* [```b1cb2ce```](https://github.com/ktmeaton/NCBImeta/commit/b1cb2ce) Remove test warning [skip ci]
* [```994f4b1```](https://github.com/ktmeaton/NCBImeta/commit/994f4b1) Web portal links, src dir fix [skip ci]
* [```359f572```](https://github.com/ktmeaton/NCBImeta/commit/359f572) Extra quotes remove [skip ci]
* [```8edd0a2```](https://github.com/ktmeaton/NCBImeta/commit/8edd0a2) Fix module import and author credit
* [```01f4592```](https://github.com/ktmeaton/NCBImeta/commit/01f4592) fix quote
* [```b477b42```](https://github.com/ktmeaton/NCBImeta/commit/b477b42) missing comma
* [```02f81dc```](https://github.com/ktmeaton/NCBImeta/commit/02f81dc) pip install match PyPI Repo [skip ci]
* [```2c856e2```](https://github.com/ktmeaton/NCBImeta/commit/2c856e2) Remove underscores
* [```c6cd147```](https://github.com/ktmeaton/NCBImeta/commit/c6cd147) Remove src dir links [skip ci]
* [```6c5ba1e```](https://github.com/ktmeaton/NCBImeta/commit/6c5ba1e) Install update [skip ci]
* [```473ca4c```](https://github.com/ktmeaton/NCBImeta/commit/473ca4c) Home url and markdown desc [skip ci]
* [```b8b71bd```](https://github.com/ktmeaton/NCBImeta/commit/b8b71bd) Short description [skip CI]
* [```82a05e6```](https://github.com/ktmeaton/NCBImeta/commit/82a05e6) pip install included
* [```05cb9a0```](https://github.com/ktmeaton/NCBImeta/commit/05cb9a0) branch restrict rearrange
* [```f351032```](https://github.com/ktmeaton/NCBImeta/commit/f351032) travis troubleshooting part 2
* [```d9f5843```](https://github.com/ktmeaton/NCBImeta/commit/d9f5843) travis troubleshooting
* [```f574720```](https://github.com/ktmeaton/NCBImeta/commit/f574720) requirements note
* [```a411aea```](https://github.com/ktmeaton/NCBImeta/commit/a411aea) v0.4.2 starts, PyPI
* [```bc57212```](https://github.com/ktmeaton/NCBImeta/commit/bc57212) Merge branch 'dev'
* [```9909a13```](https://github.com/ktmeaton/NCBImeta/commit/9909a13) PyPI Preparation
* [```2756851```](https://github.com/ktmeaton/NCBImeta/commit/2756851) setup.py install works!
* [```b549268```](https://github.com/ktmeaton/NCBImeta/commit/b549268) Rename src dir to ncbimeta
* [```d25112a```](https://github.com/ktmeaton/NCBImeta/commit/d25112a) More setup.py config
* [```4236cec```](https://github.com/ktmeaton/NCBImeta/commit/4236cec) requirements.txt control for setup
* [```0978441```](https://github.com/ktmeaton/NCBImeta/commit/0978441) Added setup.py
* [```ff158f0```](https://github.com/ktmeaton/NCBImeta/commit/ff158f0) Remove underscores in program
* [```f074685```](https://github.com/ktmeaton/NCBImeta/commit/f074685) Merge branch 'master' into dev
* [```0f33450```](https://github.com/ktmeaton/NCBImeta/commit/0f33450) Run Travis CI only on master
* [```a7d9038```](https://github.com/ktmeaton/NCBImeta/commit/a7d9038) v0.4.1 version update

## v0.4.1

### Commits

* [```69941ad```](https://github.com/ktmeaton/NCBImeta/commit/69941ad) v0.4.1 release
* [```abf1e03```](https://github.com/ktmeaton/NCBImeta/commit/abf1e03) Rearrange and simplify [skip ci]
* [```92c7b60```](https://github.com/ktmeaton/NCBImeta/commit/92c7b60) Update CLI gif to v0.4.0 [skip ci]
* [```86b5eb5```](https://github.com/ktmeaton/NCBImeta/commit/86b5eb5) relative path change gif [skip ci]
* [```1479a5a```](https://github.com/ktmeaton/NCBImeta/commit/1479a5a) Merge branch 'dev'
* [```f8112bd```](https://github.com/ktmeaton/NCBImeta/commit/f8112bd) url update [skip ci]
* [```770f7e5```](https://github.com/ktmeaton/NCBImeta/commit/770f7e5) gh pages url [skip ci]
* [```00c72f4```](https://github.com/ktmeaton/NCBImeta/commit/00c72f4) Merge branch 'master' into dev
* [```46931e0```](https://github.com/ktmeaton/NCBImeta/commit/46931e0) Dependency clarification
* [```83ceb7f```](https://github.com/ktmeaton/NCBImeta/commit/83ceb7f) Merge branch 'dev'
* [```ac09f72```](https://github.com/ktmeaton/NCBImeta/commit/ac09f72) Multi-python linux
* [```6731c99```](https://github.com/ktmeaton/NCBImeta/commit/6731c99) Remove accessory program section [skip ci]
* [```50121e9```](https://github.com/ktmeaton/NCBImeta/commit/50121e9) Spacing [skip ci]
* [```409dfcc```](https://github.com/ktmeaton/NCBImeta/commit/409dfcc) Describe OS support [skip ci]
* [```6defc7b```](https://github.com/ktmeaton/NCBImeta/commit/6defc7b) Description simplify [skip ci]
* [```790dd77```](https://github.com/ktmeaton/NCBImeta/commit/790dd77) Links and issues doc
* [```9c2de04```](https://github.com/ktmeaton/NCBImeta/commit/9c2de04) Update program description
* [```da2c18c```](https://github.com/ktmeaton/NCBImeta/commit/da2c18c) Simple Linux and MacOS Build
* [```04151ed```](https://github.com/ktmeaton/NCBImeta/commit/04151ed) Remove URL Error code print
* [```c67c9d5```](https://github.com/ktmeaton/NCBImeta/commit/c67c9d5) Travis-CI MacOS Supported!
* [```f22cedc```](https://github.com/ktmeaton/NCBImeta/commit/f22cedc) esearch cmd rearrange
* [```c3b3444```](https://github.com/ktmeaton/NCBImeta/commit/c3b3444) kwargs typo
* [```df09ab4```](https://github.com/ktmeaton/NCBImeta/commit/df09ab4) Typo http_metthod
* [```4e38e92```](https://github.com/ktmeaton/NCBImeta/commit/4e38e92) Try to catch URL Error
* [```e8b43d9```](https://github.com/ktmeaton/NCBImeta/commit/e8b43d9) Retry macos testing
* [```fb9c1e7```](https://github.com/ktmeaton/NCBImeta/commit/fb9c1e7) Macos test first
* [```20762ef```](https://github.com/ktmeaton/NCBImeta/commit/20762ef) Pip3 change and multi-python test
* [```19eb7e9```](https://github.com/ktmeaton/NCBImeta/commit/19eb7e9) Bugfix plus exclude windwos from travis
* [```8571960```](https://github.com/ktmeaton/NCBImeta/commit/8571960) Test with python3 prefix
* [```d4a1b51```](https://github.com/ktmeaton/NCBImeta/commit/d4a1b51) Windows travis-ci test
* [```5e8e305```](https://github.com/ktmeaton/NCBImeta/commit/5e8e305) Py 3.7.4 multi-dist test
* [```65be703```](https://github.com/ktmeaton/NCBImeta/commit/65be703) Add Travis-CI Badge [skip ci]
* [```0f8c65f```](https://github.com/ktmeaton/NCBImeta/commit/0f8c65f) Just test python 3.6 on dev
* [```be9fc9b```](https://github.com/ktmeaton/NCBImeta/commit/be9fc9b) v0.4.1 preparation [skip ci]
* [```5300621```](https://github.com/ktmeaton/NCBImeta/commit/5300621) Annotation script bugfixes [skip ci]
* [```c021844```](https://github.com/ktmeaton/NCBImeta/commit/c021844) Updated biopython requirement to 1.74 [skip ci]
* [```535aa7f```](https://github.com/ktmeaton/NCBImeta/commit/535aa7f) Add numpy to requirements although packaged with biopython
* [```b837b59```](https://github.com/ktmeaton/NCBImeta/commit/b837b59) Only dev branches for travis-ci
* [```56105f9```](https://github.com/ktmeaton/NCBImeta/commit/56105f9) Add execute permissions to scripts
* [```d4f8c29```](https://github.com/ktmeaton/NCBImeta/commit/d4f8c29) Remove --user flag from travis-ci pip
* [```9a436e4```](https://github.com/ktmeaton/NCBImeta/commit/9a436e4) Travis CI Integration
* [```df2ddfd```](https://github.com/ktmeaton/NCBImeta/commit/df2ddfd) Merge branch 'master' into dev
* [```3b6276d```](https://github.com/ktmeaton/NCBImeta/commit/3b6276d) Simplify dependency written desc

## v0.4.0

### Commits

* [```c7b320d```](https://github.com/ktmeaton/NCBImeta/commit/c7b320d) v0.4.0 release changes
* [```5a4a799```](https://github.com/ktmeaton/NCBImeta/commit/5a4a799) Merge branch 'dev'
* [```1752716```](https://github.com/ktmeaton/NCBImeta/commit/1752716) Variable formatting
* [```3a910cb```](https://github.com/ktmeaton/NCBImeta/commit/3a910cb) Formatting test
* [```70ee5cb```](https://github.com/ktmeaton/NCBImeta/commit/70ee5cb) Improved parameter description
* [```882140a```](https://github.com/ktmeaton/NCBImeta/commit/882140a) annot file updates
* [```3c257cd```](https://github.com/ktmeaton/NCBImeta/commit/3c257cd) Fixed program description.
* [```92bee13```](https://github.com/ktmeaton/NCBImeta/commit/92bee13) Improved description
* [```6471de7```](https://github.com/ktmeaton/NCBImeta/commit/6471de7) Additional explanation for join and export
* [```9a0e195```](https://github.com/ktmeaton/NCBImeta/commit/9a0e195) Section rearrange
* [```0ccfe8a```](https://github.com/ktmeaton/NCBImeta/commit/0ccfe8a) Fix incorrect release version links
* [```8ae3a9c```](https://github.com/ktmeaton/NCBImeta/commit/8ae3a9c) more typos
* [```d1df6f8```](https://github.com/ktmeaton/NCBImeta/commit/d1df6f8) Typos etc
* [```8780293```](https://github.com/ktmeaton/NCBImeta/commit/8780293) Release links
* [```82eaccc```](https://github.com/ktmeaton/NCBImeta/commit/82eaccc) Main README simplify
* [```28c28d7```](https://github.com/ktmeaton/NCBImeta/commit/28c28d7) Force metadata concatenation even if same value
* [```acb33f8```](https://github.com/ktmeaton/NCBImeta/commit/acb33f8) Simply annotation file name
* [```32de662```](https://github.com/ktmeaton/NCBImeta/commit/32de662) Print msg upon successful match
* [```edca4c0```](https://github.com/ktmeaton/NCBImeta/commit/edca4c0) Rmv records outside time boundary
* [```f42d0a6```](https://github.com/ktmeaton/NCBImeta/commit/f42d0a6) Rmv annot file outside time boundary
* [```792d683```](https://github.com/ktmeaton/NCBImeta/commit/792d683) Output msg for succesful match
* [```047124b```](https://github.com/ktmeaton/NCBImeta/commit/047124b) Schema updates
* [```a33cc25```](https://github.com/ktmeaton/NCBImeta/commit/a33cc25) Removed scripts folder with R code
* [```930f2e5```](https://github.com/ktmeaton/NCBImeta/commit/930f2e5) Clarify whitespace importance
* [```bb86c43```](https://github.com/ktmeaton/NCBImeta/commit/bb86c43) Remove unsupported schema txt files
* [```5e2a608```](https://github.com/ktmeaton/NCBImeta/commit/5e2a608) Doc formatting
* [```7008c66```](https://github.com/ktmeaton/NCBImeta/commit/7008c66) Doc formatting and spacing
* [```a1c5bed```](https://github.com/ktmeaton/NCBImeta/commit/a1c5bed) Schema documentation update for yaml
* [```b7edcea```](https://github.com/ktmeaton/NCBImeta/commit/b7edcea) schema files in yaml format
* [```4762bb9```](https://github.com/ktmeaton/NCBImeta/commit/4762bb9) Added comment field to Pubmed table
* [```e431a2b```](https://github.com/ktmeaton/NCBImeta/commit/e431a2b) badge typo and gif to be updated
* [```52ad001```](https://github.com/ktmeaton/NCBImeta/commit/52ad001) Sanity commit
* [```2300524```](https://github.com/ktmeaton/NCBImeta/commit/2300524) update, MWE back to plague
* [```0014b8c```](https://github.com/ktmeaton/NCBImeta/commit/0014b8c) MWE organism change and direct execution without python3 call
* [```4ebc22b```](https://github.com/ktmeaton/NCBImeta/commit/4ebc22b) Shebang updates for source files
* [```e81f17e```](https://github.com/ktmeaton/NCBImeta/commit/e81f17e) Annotation file rename
* [```a883af4```](https://github.com/ktmeaton/NCBImeta/commit/a883af4) Switch MWE back to plague for shorter time
* [```7dc0cb4```](https://github.com/ktmeaton/NCBImeta/commit/7dc0cb4) Remove old py config file
* [```5320f67```](https://github.com/ktmeaton/NCBImeta/commit/5320f67) Slowdown fetch as default
* [```c8973be```](https://github.com/ktmeaton/NCBImeta/commit/c8973be) Shebang typo
* [```58829cc```](https://github.com/ktmeaton/NCBImeta/commit/58829cc) Implement read handle error checking
* [```035d0fd```](https://github.com/ktmeaton/NCBImeta/commit/035d0fd) Error class for read errors
* [```1df91ae```](https://github.com/ktmeaton/NCBImeta/commit/1df91ae) add db name
* [```25aff7c```](https://github.com/ktmeaton/NCBImeta/commit/25aff7c) Switch example config file to P. aeruginosa to match paper
* [```d50e012```](https://github.com/ktmeaton/NCBImeta/commit/d50e012) v0.4.0 update rename
* [```91b3d62```](https://github.com/ktmeaton/NCBImeta/commit/91b3d62) HTTP 429 error checking for efetch
* [```906b10c```](https://github.com/ktmeaton/NCBImeta/commit/906b10c) Example command now references config.yaml
* [```99ed046```](https://github.com/ktmeaton/NCBImeta/commit/99ed046) Updates for v0.3.4 and v0.3.5
* [```d648328```](https://github.com/ktmeaton/NCBImeta/commit/d648328) Implement hiearchical fields
* [```36639a2```](https://github.com/ktmeaton/NCBImeta/commit/36639a2) Update hierarchical docs
* [```66ab5c0```](https://github.com/ktmeaton/NCBImeta/commit/66ab5c0) comma space sep hierarchical fields
* [```6da1485```](https://github.com/ktmeaton/NCBImeta/commit/6da1485) Removed extra config files
* [```e4f36d3```](https://github.com/ktmeaton/NCBImeta/commit/e4f36d3) Add requirements file
* [```d772064```](https://github.com/ktmeaton/NCBImeta/commit/d772064) Implement config file as yaml
* [```91cd02d```](https://github.com/ktmeaton/NCBImeta/commit/91cd02d) Error for incorrect config parameters
* [```b6966d8```](https://github.com/ktmeaton/NCBImeta/commit/b6966d8) config file in yaml format
* [```a16f969```](https://github.com/ktmeaton/NCBImeta/commit/a16f969) Merge branch 'master' into dev

## v0.3.4

### Commits

* [```da80053```](https://github.com/ktmeaton/NCBImeta/commit/da80053) Update
* [```f53bdac```](https://github.com/ktmeaton/NCBImeta/commit/f53bdac) Merge branch 'dev'
* [```297f12f```](https://github.com/ktmeaton/NCBImeta/commit/297f12f) Document AnnotateReplace and AnnotateConcatenate
* [```18cf30a```](https://github.com/ktmeaton/NCBImeta/commit/18cf30a) Document FORCE_PAUSE_SECONDS parameter
* [```8bf1737```](https://github.com/ktmeaton/NCBImeta/commit/8bf1737) Include all 6 tables in config files
* [```424e4a6```](https://github.com/ktmeaton/NCBImeta/commit/424e4a6) Create a parameter to allow force pausing in seconds
* [```09425f4```](https://github.com/ktmeaton/NCBImeta/commit/09425f4) Implemented NCBI API key acceptance
* [```9927043```](https://github.com/ktmeaton/NCBImeta/commit/9927043) Add P. aeruginosa config file
* [```3bc04b3```](https://github.com/ktmeaton/NCBImeta/commit/3bc04b3) Added error when max fetch records is exceeded
* [```5464bc4```](https://github.com/ktmeaton/NCBImeta/commit/5464bc4) Bugfix for HTTPError 429 checking
* [```8bb603e```](https://github.com/ktmeaton/NCBImeta/commit/8bb603e) HTTP Error catching for Entrez esummary
* [```aad5f39```](https://github.com/ktmeaton/NCBImeta/commit/aad5f39) shebang change to python3
* [```3ebd8e0```](https://github.com/ktmeaton/NCBImeta/commit/3ebd8e0) Set theme jekyll-theme-slate
* [```578e473```](https://github.com/ktmeaton/NCBImeta/commit/578e473) Set theme jekyll-theme-minimal
* [```b198163```](https://github.com/ktmeaton/NCBImeta/commit/b198163) Goal: Config file as yaml format.

## v0.3.3

### Commits

* [```2c8f1d5```](https://github.com/ktmeaton/NCBImeta/commit/2c8f1d5) v0.3.3 Release
* [```b792547```](https://github.com/ktmeaton/NCBImeta/commit/b792547) Merge branch 'dev'
* [```2797444```](https://github.com/ktmeaton/NCBImeta/commit/2797444) README update
* [```3c557c6```](https://github.com/ktmeaton/NCBImeta/commit/3c557c6) Force python3 command
* [```5e10baa```](https://github.com/ktmeaton/NCBImeta/commit/5e10baa) Changelog includes Pubmed support
* [```7e2196c```](https://github.com/ktmeaton/NCBImeta/commit/7e2196c) Pubmed add and successful example test
* [```f12a920```](https://github.com/ktmeaton/NCBImeta/commit/f12a920) Schema update and Pubmed addition
* [```83e321e```](https://github.com/ktmeaton/NCBImeta/commit/83e321e) Pubmed Table functionality added
* [```ea5e676```](https://github.com/ktmeaton/NCBImeta/commit/ea5e676) Added a str conversion to xml parsing step for Pubmed table
* [```d08dc50```](https://github.com/ktmeaton/NCBImeta/commit/d08dc50) v0.3.3 development begins

## v0.3.2

### Commits

* [```f800a8a```](https://github.com/ktmeaton/NCBImeta/commit/f800a8a) Changelog update
* [```39ab0e2```](https://github.com/ktmeaton/NCBImeta/commit/39ab0e2) Delete my_organism_annot.txt
* [```86abec3```](https://github.com/ktmeaton/NCBImeta/commit/86abec3) New annotation files for exampel
* [```8eeb4be```](https://github.com/ktmeaton/NCBImeta/commit/8eeb4be) NCBImeta_Export converted to python3
* [```6ec6850```](https://github.com/ktmeaton/NCBImeta/commit/6ec6850) Missing DB Error added
* [```315a388```](https://github.com/ktmeaton/NCBImeta/commit/315a388) Last commit before v0.3.3
* [```63141bd```](https://github.com/ktmeaton/NCBImeta/commit/63141bd) Nucleotide annotations fully fixed
* [```0d3db34```](https://github.com/ktmeaton/NCBImeta/commit/0d3db34) Nucleotide schema CDSs pluralize
* [```65b03ee```](https://github.com/ktmeaton/NCBImeta/commit/65b03ee) Nucleotide Table accession fix
* [```1fe7dc1```](https://github.com/ktmeaton/NCBImeta/commit/1fe7dc1) Fix Nucleotide annotations for py3
* [```36027f7```](https://github.com/ktmeaton/NCBImeta/commit/36027f7) Remove utf encode for py3 functionality
* [```c3923ca```](https://github.com/ktmeaton/NCBImeta/commit/c3923ca) Unicode decoding for Python3
* [```359a2f7```](https://github.com/ktmeaton/NCBImeta/commit/359a2f7) dict_items error fix
* [```7c10cc8```](https://github.com/ktmeaton/NCBImeta/commit/7c10cc8) dev branch reinstated
* [```d7603d5```](https://github.com/ktmeaton/NCBImeta/commit/d7603d5) Add XPath and XLST to wishlist
* [```c2838e8```](https://github.com/ktmeaton/NCBImeta/commit/c2838e8) Update README.md
* [```d21c83e```](https://github.com/ktmeaton/NCBImeta/commit/d21c83e) Updated instructions to reflect no dev branch yet
* [```4f37c23```](https://github.com/ktmeaton/NCBImeta/commit/4f37c23) Update development branch names to dev
* [```0472be9```](https://github.com/ktmeaton/NCBImeta/commit/0472be9) Update README.md
* [```689a545```](https://github.com/ktmeaton/NCBImeta/commit/689a545) First successful completion of the master join table
* [```a21171d```](https://github.com/ktmeaton/NCBImeta/commit/a21171d) Join script beginning functional
* [```28f4724```](https://github.com/ktmeaton/NCBImeta/commit/28f4724) Major fix in Nucleotide metadata retrieval, properly recovering annotations
* [```22aa44d```](https://github.com/ktmeaton/NCBImeta/commit/22aa44d) Little things in README
* [```8d2f0cf```](https://github.com/ktmeaton/NCBImeta/commit/8d2f0cf) example update
* [```1a852f1```](https://github.com/ktmeaton/NCBImeta/commit/1a852f1) Rename default column names so that they are all unique
* [```e53deb1```](https://github.com/ktmeaton/NCBImeta/commit/e53deb1) Metadata annotation scripts now allow concatenation
* [```5cdd9bb```](https://github.com/ktmeaton/NCBImeta/commit/5cdd9bb) Change annotation script to 2 options, replace or concatenate
* [```6f9247c```](https://github.com/ktmeaton/NCBImeta/commit/6f9247c) Use harmonized_name attribute for BioSample fields
* [```bffefee```](https://github.com/ktmeaton/NCBImeta/commit/bffefee) Remove excess config files and example output
* [```c45fc56```](https://github.com/ktmeaton/NCBImeta/commit/c45fc56) Slow down NCBI fetch requests to 1 per second
* [```a045953```](https://github.com/ktmeaton/NCBImeta/commit/a045953) Update README.md
* [```c90c644```](https://github.com/ktmeaton/NCBImeta/commit/c90c644) Update README.md
* [```b731565```](https://github.com/ktmeaton/NCBImeta/commit/b731565) Update README.md
* [```5993154```](https://github.com/ktmeaton/NCBImeta/commit/5993154) GIFs
* [```b21cbfd```](https://github.com/ktmeaton/NCBImeta/commit/b21cbfd) Add files via upload
* [```8cd4cbe```](https://github.com/ktmeaton/NCBImeta/commit/8cd4cbe) Encode rather than unicode, experiment

## v0.3.1

### Commits

* [```cd8fda6```](https://github.com/ktmeaton/NCBImeta/commit/cd8fda6) Update README.md
* [```e6da625```](https://github.com/ktmeaton/NCBImeta/commit/e6da625) Update README_schema.md
* [```2d96701```](https://github.com/ktmeaton/NCBImeta/commit/2d96701) Update README_schema.md
* [```57ff534```](https://github.com/ktmeaton/NCBImeta/commit/57ff534) Update README_schema.md
* [```38e4861```](https://github.com/ktmeaton/NCBImeta/commit/38e4861) Update README_schema.md
* [```13d25e7```](https://github.com/ktmeaton/NCBImeta/commit/13d25e7) Update README_schema.md
* [```93b830d```](https://github.com/ktmeaton/NCBImeta/commit/93b830d) Update and rename README_schema.txt to README_schema.md
* [```4b80683```](https://github.com/ktmeaton/NCBImeta/commit/4b80683) Update README_config.md
* [```dc775f5```](https://github.com/ktmeaton/NCBImeta/commit/dc775f5) Update README_config.md
* [```766c815```](https://github.com/ktmeaton/NCBImeta/commit/766c815) Update README_config.md
* [```fd2a4f4```](https://github.com/ktmeaton/NCBImeta/commit/fd2a4f4) Update README_config.md
* [```328a571```](https://github.com/ktmeaton/NCBImeta/commit/328a571) Update README_config.md
* [```ee9eedc```](https://github.com/ktmeaton/NCBImeta/commit/ee9eedc) Update and rename README_config.txt to README_config.md
* [```b6d4bac```](https://github.com/ktmeaton/NCBImeta/commit/b6d4bac) Update README.md
* [```047b0af```](https://github.com/ktmeaton/NCBImeta/commit/047b0af) Move images
* [```e5c3f02```](https://github.com/ktmeaton/NCBImeta/commit/e5c3f02) Moved images
* [```87a24a7```](https://github.com/ktmeaton/NCBImeta/commit/87a24a7) Update README.md
* [```77148b5```](https://github.com/ktmeaton/NCBImeta/commit/77148b5) Update README.md
* [```104ef14```](https://github.com/ktmeaton/NCBImeta/commit/104ef14) Update README.md
* [```5734943```](https://github.com/ktmeaton/NCBImeta/commit/5734943) Update README.md
* [```09ec0dc```](https://github.com/ktmeaton/NCBImeta/commit/09ec0dc) Update README.md
* [```4085817```](https://github.com/ktmeaton/NCBImeta/commit/4085817) Update README.md
* [```69347c1```](https://github.com/ktmeaton/NCBImeta/commit/69347c1) Update README.md
* [```c3488df```](https://github.com/ktmeaton/NCBImeta/commit/c3488df) Update README.md
* [```e3ac31c```](https://github.com/ktmeaton/NCBImeta/commit/e3ac31c) Update README.md
* [```3551c72```](https://github.com/ktmeaton/NCBImeta/commit/3551c72) Add files via upload
* [```3963afc```](https://github.com/ktmeaton/NCBImeta/commit/3963afc) Update README.md
* [```d49de82```](https://github.com/ktmeaton/NCBImeta/commit/d49de82) Update README.md
* [```25bf63e```](https://github.com/ktmeaton/NCBImeta/commit/25bf63e) Add files via upload
* [```ac30fd2```](https://github.com/ktmeaton/NCBImeta/commit/ac30fd2) Delete NCBImeta_screenshot.png
* [```f3afa94```](https://github.com/ktmeaton/NCBImeta/commit/f3afa94) Update README.md
* [```3a393d0```](https://github.com/ktmeaton/NCBImeta/commit/3a393d0) higher
* [```5349f1a```](https://github.com/ktmeaton/NCBImeta/commit/5349f1a) merge
* [```e2a132b```](https://github.com/ktmeaton/NCBImeta/commit/e2a132b) higher res
* [```d65ee34```](https://github.com/ktmeaton/NCBImeta/commit/d65ee34) Delete NCBImeta_snapshot.jpg
* [```82a1350```](https://github.com/ktmeaton/NCBImeta/commit/82a1350) Delete NCBImeta_screenshot.png
* [```a15b804```](https://github.com/ktmeaton/NCBImeta/commit/a15b804) Update README.md
* [```693012c```](https://github.com/ktmeaton/NCBImeta/commit/693012c) Update README.md
* [```f9b1169```](https://github.com/ktmeaton/NCBImeta/commit/f9b1169) Update README.md
* [```97e5857```](https://github.com/ktmeaton/NCBImeta/commit/97e5857) Update README.md
* [```489bd94```](https://github.com/ktmeaton/NCBImeta/commit/489bd94) Update README.md
* [```01670d9```](https://github.com/ktmeaton/NCBImeta/commit/01670d9) Update README.md
* [```7118d27```](https://github.com/ktmeaton/NCBImeta/commit/7118d27) Update README.md
* [```58dce4c```](https://github.com/ktmeaton/NCBImeta/commit/58dce4c) Update README.md
* [```dd7c5a0```](https://github.com/ktmeaton/NCBImeta/commit/dd7c5a0) Update README.md
* [```ddd7857```](https://github.com/ktmeaton/NCBImeta/commit/ddd7857) Update README.md
* [```df34ed8```](https://github.com/ktmeaton/NCBImeta/commit/df34ed8) Update README.md
* [```94d1827```](https://github.com/ktmeaton/NCBImeta/commit/94d1827) Update README.md
* [```3bd6f24```](https://github.com/ktmeaton/NCBImeta/commit/3bd6f24) Update README.md
* [```c84af73```](https://github.com/ktmeaton/NCBImeta/commit/c84af73) Update README.md
* [```037833f```](https://github.com/ktmeaton/NCBImeta/commit/037833f) Update README.md
* [```24a010b```](https://github.com/ktmeaton/NCBImeta/commit/24a010b) Update README.md
* [```830ca96```](https://github.com/ktmeaton/NCBImeta/commit/830ca96) Image add
* [```e906acc```](https://github.com/ktmeaton/NCBImeta/commit/e906acc) Update README.md
* [```719e986```](https://github.com/ktmeaton/NCBImeta/commit/719e986) Update README.md
* [```1f590d1```](https://github.com/ktmeaton/NCBImeta/commit/1f590d1) Update README.md
* [```35eb86e```](https://github.com/ktmeaton/NCBImeta/commit/35eb86e) Update README.md
* [```331489a```](https://github.com/ktmeaton/NCBImeta/commit/331489a) Update README.md
* [```a648589```](https://github.com/ktmeaton/NCBImeta/commit/a648589) Update README.md
* [```a6be405```](https://github.com/ktmeaton/NCBImeta/commit/a6be405) Update README.md
* [```18cb81a```](https://github.com/ktmeaton/NCBImeta/commit/18cb81a) Update README.md
* [```0ca9547```](https://github.com/ktmeaton/NCBImeta/commit/0ca9547) Update README.md
* [```65b290f```](https://github.com/ktmeaton/NCBImeta/commit/65b290f) Update README.md
* [```0081780```](https://github.com/ktmeaton/NCBImeta/commit/0081780) Update README.md
* [```ede0ebe```](https://github.com/ktmeaton/NCBImeta/commit/ede0ebe) example dir, more unicode testing
* [```13b657b```](https://github.com/ktmeaton/NCBImeta/commit/13b657b) Bug fixing more unicode
* [```de6e9bb```](https://github.com/ktmeaton/NCBImeta/commit/de6e9bb) Update NCBImeta.py
* [```c405178```](https://github.com/ktmeaton/NCBImeta/commit/c405178) Update README.md
* [```7c13655```](https://github.com/ktmeaton/NCBImeta/commit/7c13655) Update README.md
* [```aa71c05```](https://github.com/ktmeaton/NCBImeta/commit/aa71c05) Update README.md
* [```4df6f8f```](https://github.com/ktmeaton/NCBImeta/commit/4df6f8f) Update README.md
* [```46041ed```](https://github.com/ktmeaton/NCBImeta/commit/46041ed) Update README.md
* [```fdc0029```](https://github.com/ktmeaton/NCBImeta/commit/fdc0029) Update README.md
* [```9ffcf42```](https://github.com/ktmeaton/NCBImeta/commit/9ffcf42) Update README.md
* [```cf734b6```](https://github.com/ktmeaton/NCBImeta/commit/cf734b6) Update README.md
* [```e2937a4```](https://github.com/ktmeaton/NCBImeta/commit/e2937a4) Update README.md
* [```95228b1```](https://github.com/ktmeaton/NCBImeta/commit/95228b1) Update README.md
* [```3921724```](https://github.com/ktmeaton/NCBImeta/commit/3921724) Update README.md
* [```e74a81a```](https://github.com/ktmeaton/NCBImeta/commit/e74a81a) Update README.md
* [```3fc5925```](https://github.com/ktmeaton/NCBImeta/commit/3fc5925) Update README.md
* [```88dd3d8```](https://github.com/ktmeaton/NCBImeta/commit/88dd3d8) Update README.md
* [```c406454```](https://github.com/ktmeaton/NCBImeta/commit/c406454) Update README.md
* [```611501b```](https://github.com/ktmeaton/NCBImeta/commit/611501b) Update README.md
* [```06193bf```](https://github.com/ktmeaton/NCBImeta/commit/06193bf) Update README.md
* [```65594c5```](https://github.com/ktmeaton/NCBImeta/commit/65594c5) Update README.md
* [```a33ef41```](https://github.com/ktmeaton/NCBImeta/commit/a33ef41) Update README.md
* [```ee70c20```](https://github.com/ktmeaton/NCBImeta/commit/ee70c20) Update README.md
* [```c0ba54c```](https://github.com/ktmeaton/NCBImeta/commit/c0ba54c) Update README.md
* [```b956374```](https://github.com/ktmeaton/NCBImeta/commit/b956374) Update README.md
* [```f8bd00d```](https://github.com/ktmeaton/NCBImeta/commit/f8bd00d) Update README.md
* [```5d504ec```](https://github.com/ktmeaton/NCBImeta/commit/5d504ec) Update README.md
* [```a4a7a42```](https://github.com/ktmeaton/NCBImeta/commit/a4a7a42) Update README.md
* [```576f7df```](https://github.com/ktmeaton/NCBImeta/commit/576f7df) Update README.md
* [```6c3c5f6```](https://github.com/ktmeaton/NCBImeta/commit/6c3c5f6) Update README.md
* [```4716032```](https://github.com/ktmeaton/NCBImeta/commit/4716032) Update README.md
* [```5360f6d```](https://github.com/ktmeaton/NCBImeta/commit/5360f6d) Badge Testing
* [```c268d12```](https://github.com/ktmeaton/NCBImeta/commit/c268d12) Update geocode.R
* [```21dcccd```](https://github.com/ktmeaton/NCBImeta/commit/21dcccd) Update CHANGELOG.md
* [```995b702```](https://github.com/ktmeaton/NCBImeta/commit/995b702) Update CHANGELOG.md
* [```7cdae9c```](https://github.com/ktmeaton/NCBImeta/commit/7cdae9c) Update CHANGELOG.md
* [```e2d363d```](https://github.com/ktmeaton/NCBImeta/commit/e2d363d) Merge branch 'master' of https://github.com/ktmeaton/NCBImeta
* [```2b24140```](https://github.com/ktmeaton/NCBImeta/commit/2b24140) Remove R tmp
* [```d98765f```](https://github.com/ktmeaton/NCBImeta/commit/d98765f) Delete .Rhistory
* [```97df6dc```](https://github.com/ktmeaton/NCBImeta/commit/97df6dc) Ignore .pyc
* [```b759f15```](https://github.com/ktmeaton/NCBImeta/commit/b759f15) Cleanup
* [```ab545b4```](https://github.com/ktmeaton/NCBImeta/commit/ab545b4) cleanup .pyc
* [```7a13369```](https://github.com/ktmeaton/NCBImeta/commit/7a13369) License and setup
* [```7118d10```](https://github.com/ktmeaton/NCBImeta/commit/7118d10) Set theme jekyll-theme-slate
* [```b946a24```](https://github.com/ktmeaton/NCBImeta/commit/b946a24) Delete _config.yml
* [```17af0ac```](https://github.com/ktmeaton/NCBImeta/commit/17af0ac) Update _config.yml
* [```c9c7331```](https://github.com/ktmeaton/NCBImeta/commit/c9c7331) Set theme jekyll-theme-architect
* [```eb87a53```](https://github.com/ktmeaton/NCBImeta/commit/eb87a53) Update
* [```5fa31d5```](https://github.com/ktmeaton/NCBImeta/commit/5fa31d5) Update _config.yml
* [```09b4b07```](https://github.com/ktmeaton/NCBImeta/commit/09b4b07) Update _config.yml
* [```2471e05```](https://github.com/ktmeaton/NCBImeta/commit/2471e05) Update _config.yml
* [```63caeaa```](https://github.com/ktmeaton/NCBImeta/commit/63caeaa) Update _config.yml
* [```50b0045```](https://github.com/ktmeaton/NCBImeta/commit/50b0045) Initial config
* [```3dd8ba9```](https://github.com/ktmeaton/NCBImeta/commit/3dd8ba9) GH Pages
* [```0d15ea9```](https://github.com/ktmeaton/NCBImeta/commit/0d15ea9) Set theme jekyll-theme-architect
* [```7aa3cf5```](https://github.com/ktmeaton/NCBImeta/commit/7aa3cf5) Set theme jekyll-theme-architect
* [```02115ff```](https://github.com/ktmeaton/NCBImeta/commit/02115ff) Python 2/3 compatability fix for unicode
* [```e3faf17```](https://github.com/ktmeaton/NCBImeta/commit/e3faf17) Update README.md
* [```42d82a1```](https://github.com/ktmeaton/NCBImeta/commit/42d82a1) Update requirements
* [```aed30d8```](https://github.com/ktmeaton/NCBImeta/commit/aed30d8) v0.3.2 Begins
* [```0eb820c```](https://github.com/ktmeaton/NCBImeta/commit/0eb820c) Merge branch 'v0.3.1'
* [```7a23c24```](https://github.com/ktmeaton/NCBImeta/commit/7a23c24) Ignore pycache directories
* [```d570112```](https://github.com/ktmeaton/NCBImeta/commit/d570112) Ignore
* [```4e90e7f```](https://github.com/ktmeaton/NCBImeta/commit/4e90e7f) Cleanup
* [```a606819```](https://github.com/ktmeaton/NCBImeta/commit/a606819) Cleanup
* [```c541c9b```](https://github.com/ktmeaton/NCBImeta/commit/c541c9b) READMEs, bug fix on error raising
* [```6bddcde```](https://github.com/ktmeaton/NCBImeta/commit/6bddcde) remove kits
* [```9d21d46```](https://github.com/ktmeaton/NCBImeta/commit/9d21d46) v0.3.1 begins
* [```1c97b6b```](https://github.com/ktmeaton/NCBImeta/commit/1c97b6b) Merge branch 'v0.3.0'
* [```d78b88b```](https://github.com/ktmeaton/NCBImeta/commit/d78b88b) merge with master
* [```ed26911```](https://github.com/ktmeaton/NCBImeta/commit/ed26911) Final commit of v0.3.0
* [```df8cedb```](https://github.com/ktmeaton/NCBImeta/commit/df8cedb) Fixed node-attribute confliction
* [```c57bafa```](https://github.com/ktmeaton/NCBImeta/commit/c57bafa) Cleanup old source files
* [```b069d4f```](https://github.com/ktmeaton/NCBImeta/commit/b069d4f) Proper unicode re-encoding
* [```5f673b7```](https://github.com/ktmeaton/NCBImeta/commit/5f673b7) Function Assembly, Biosample, BioProject, SRA
* [```e49648e```](https://github.com/ktmeaton/NCBImeta/commit/e49648e) Functional state SRA
* [```ace7c45```](https://github.com/ktmeaton/NCBImeta/commit/ace7c45) Temp save
* [```5693e3e```](https://github.com/ktmeaton/NCBImeta/commit/5693e3e) Functional automation of Assembly and BioSample
* [```7790c65```](https://github.com/ktmeaton/NCBImeta/commit/7790c65) Remove some files
* [```700bded```](https://github.com/ktmeaton/NCBImeta/commit/700bded) Rename and reorganize
* [```3bf8e94```](https://github.com/ktmeaton/NCBImeta/commit/3bf8e94) Dynamic play
* [```f39e9c6```](https://github.com/ktmeaton/NCBImeta/commit/f39e9c6) Annotate update
* [```2ff7b25```](https://github.com/ktmeaton/NCBImeta/commit/2ff7b25) Annotations
* [```8cdc64c```](https://github.com/ktmeaton/NCBImeta/commit/8cdc64c) Duplicate removal
* [```7edfbc1```](https://github.com/ktmeaton/NCBImeta/commit/7edfbc1) XML frustration
* [```26b5018```](https://github.com/ktmeaton/NCBImeta/commit/26b5018) On the fly, database creation
* [```27eb8a3```](https://github.com/ktmeaton/NCBImeta/commit/27eb8a3) v0.3.0 rewrite main script
* [```efa56ce```](https://github.com/ktmeaton/NCBImeta/commit/efa56ce) v0.3.0 Development Beings - Automated
* [```00b6f1b```](https://github.com/ktmeaton/NCBImeta/commit/00b6f1b) v0.3.0 Development Begins - Automated
* [```7e77ac2```](https://github.com/ktmeaton/NCBImeta/commit/7e77ac2) Merge branch 'v0.2.1'
* [```9383f6c```](https://github.com/ktmeaton/NCBImeta/commit/9383f6c) More annotations
* [```97f6419```](https://github.com/ktmeaton/NCBImeta/commit/97f6419) More annotation
* [```4965472```](https://github.com/ktmeaton/NCBImeta/commit/4965472) Massive annotation
* [```a09d401```](https://github.com/ktmeaton/NCBImeta/commit/a09d401) CURATED annotation files
* [```5df0f48```](https://github.com/ktmeaton/NCBImeta/commit/5df0f48) bioproject addition to annotate
* [```4b3a27f```](https://github.com/ktmeaton/NCBImeta/commit/4b3a27f) Check multiple matches in annotation, removed BioSample columns
* [```82661a4```](https://github.com/ktmeaton/NCBImeta/commit/82661a4) Consensus file after plotting
* [```710d6b5```](https://github.com/ktmeaton/NCBImeta/commit/710d6b5) First pass at plotting geographic distribution, Pandemic, PubUse
* [```d851f70```](https://github.com/ktmeaton/NCBImeta/commit/d851f70) Plotting
* [```fc2a071```](https://github.com/ktmeaton/NCBImeta/commit/fc2a071) Plotting
* [```1248978```](https://github.com/ktmeaton/NCBImeta/commit/1248978) Consensus lat,long,pandemic,established
* [```2716c6b```](https://github.com/ktmeaton/NCBImeta/commit/2716c6b) Geocoding switched to Doogal
* [```d8161ae```](https://github.com/ktmeaton/NCBImeta/commit/d8161ae) geocoding
* [```9a59fd2```](https://github.com/ktmeaton/NCBImeta/commit/9a59fd2) Yp database for plotting
* [```e5036c0```](https://github.com/ktmeaton/NCBImeta/commit/e5036c0) More annotations
* [```69a3865```](https://github.com/ktmeaton/NCBImeta/commit/69a3865) Annotation files
* [```1f09561```](https://github.com/ktmeaton/NCBImeta/commit/1f09561) Dare I say.. functional v0.2?
* [```2d4399d```](https://github.com/ktmeaton/NCBImeta/commit/2d4399d) Heuristisc approach to merge
* [```fcce755```](https://github.com/ktmeaton/NCBImeta/commit/fcce755) Biosample improve
* [```4450cd1```](https://github.com/ktmeaton/NCBImeta/commit/4450cd1) README Title Fix
* [```665b47c```](https://github.com/ktmeaton/NCBImeta/commit/665b47c) v0.2.1 README Update
* [```64467ae```](https://github.com/ktmeaton/NCBImeta/commit/64467ae) Assembly, SRA, Bioproject functional
* [```d33fc24```](https://github.com/ktmeaton/NCBImeta/commit/d33fc24) Update CHANGELOG.md
* [```38be909```](https://github.com/ktmeaton/NCBImeta/commit/38be909) v0.2.0 In Development Begins, Rename NCBInfect
* [```b010314```](https://github.com/ktmeaton/NCBImeta/commit/b010314) v0.2.0 Begins Development
* [```9dfe4d7```](https://github.com/ktmeaton/NCBImeta/commit/9dfe4d7) Merge branch 'v0.1.2'
* [```594b4bc```](https://github.com/ktmeaton/NCBImeta/commit/594b4bc) Sanity commit, Jan 2018
* [```0d6b931```](https://github.com/ktmeaton/NCBImeta/commit/0d6b931) Update CHANGELOG.md
* [```ae8d3b3```](https://github.com/ktmeaton/NCBImeta/commit/ae8d3b3) Update CHANGELOG.md
* [```ee48571```](https://github.com/ktmeaton/NCBImeta/commit/ee48571) Update CHANGELOG.md
* [```0e3e2f3```](https://github.com/ktmeaton/NCBImeta/commit/0e3e2f3) v0.1.2 changes
* [```973d416```](https://github.com/ktmeaton/NCBImeta/commit/973d416) Update README.md
* [```b7326d9```](https://github.com/ktmeaton/NCBImeta/commit/b7326d9) Update README.md
* [```6e1a700```](https://github.com/ktmeaton/NCBImeta/commit/6e1a700) Create CHANGELOG.md
* [```0805301```](https://github.com/ktmeaton/NCBImeta/commit/0805301) Moved annotate file
* [```54f194f```](https://github.com/ktmeaton/NCBImeta/commit/54f194f) Merge branch 'master' of https://github.com/ktmeaton/GenomeCollector
* [```b8af949```](https://github.com/ktmeaton/NCBImeta/commit/b8af949) Big Changes
* [```43c6869```](https://github.com/ktmeaton/NCBImeta/commit/43c6869) Annotate File
* [```3e0e006```](https://github.com/ktmeaton/NCBImeta/commit/3e0e006) python3 version fix
* [```e8e4041```](https://github.com/ktmeaton/NCBImeta/commit/e8e4041) python3 print statement fix
* [```a80db3f```](https://github.com/ktmeaton/NCBImeta/commit/a80db3f) Added extra fields in SRA table
* [```5453c8d```](https://github.com/ktmeaton/NCBImeta/commit/5453c8d) Version 2.0 fully functional
* [```2a1c06e```](https://github.com/ktmeaton/NCBImeta/commit/2a1c06e) SRA functionality as separate module
* [```7280508```](https://github.com/ktmeaton/NCBImeta/commit/7280508) Partially working version of SRA Table
* [```5c8f3e9```](https://github.com/ktmeaton/NCBImeta/commit/5c8f3e9) Added SRA function
* [```85654bc```](https://github.com/ktmeaton/NCBImeta/commit/85654bc) Added more to do in README.md
* [```2e859b1```](https://github.com/ktmeaton/NCBImeta/commit/2e859b1) Added genomeannotator log file writing
* [```c6ecaf0```](https://github.com/ktmeaton/NCBImeta/commit/c6ecaf0) Added the GenomeAnnotator tool
* [```5de9f53```](https://github.com/ktmeaton/NCBImeta/commit/5de9f53) Removed tmp file writing for record_dict testing
* [```fdbe076```](https://github.com/ktmeaton/NCBImeta/commit/fdbe076) Changed biosample Entrez parsing to validate=false
* [```bac6eb7```](https://github.com/ktmeaton/NCBImeta/commit/bac6eb7) Field AsmReleaseDate replaced with SubmissionDate
* [```642c2a2```](https://github.com/ktmeaton/NCBImeta/commit/642c2a2) Updated README.md to version 2.0
* [```506b976```](https://github.com/ktmeaton/NCBImeta/commit/506b976) Merge pull request #1 from ktmeaton/1.1
* [```fdc45fb```](https://github.com/ktmeaton/NCBImeta/commit/fdc45fb) Update README.md
* [```cde62fa```](https://github.com/ktmeaton/NCBImeta/commit/cde62fa) Update README.md
* [```a521f23```](https://github.com/ktmeaton/NCBImeta/commit/a521f23) Update README.md
* [```fc5bea4```](https://github.com/ktmeaton/NCBImeta/commit/fc5bea4) Update genomecollector.py
* [```c433ef2```](https://github.com/ktmeaton/NCBImeta/commit/c433ef2) Update genomecollector.py
* [```7736c0e```](https://github.com/ktmeaton/NCBImeta/commit/7736c0e) Add files via upload
* [```dfff42c```](https://github.com/ktmeaton/NCBImeta/commit/dfff42c) Delete genomeconcatenate.py
* [```ac643db```](https://github.com/ktmeaton/NCBImeta/commit/ac643db) Add files via upload
