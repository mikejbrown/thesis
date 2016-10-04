Symmetry improvement techniques for non-perturbative quantum field theory
=========================================================================


About
-----

This is the source distribution for Michael J. Brown's PhD thesis entitled "Symmetry improvement techniques for non-perturbative quantum field theory".

Except where otherwise noted, all parts of this thesis are Copyright © 2016 by Michael J. Brown.


License
-------

Except where otherwise noted, the following licenses apply:

- All non-code assets (LyX/LaTeX source files, accompanying bibtex files, images etc.) are licensed under a [Creative Commons Attribution-ShareAlike 4.0 International License](http://creativecommons.org/licenses/by-sa/4.0/).
  See the accompanying license file for full terms.

  ![CC-BY-SA 4.0 badge](cc-by-sa.png)

- All code (Python scripts, Mathematica notebooks, shell scripts etc.) are licensed under the [Academic Free License (“AFL”) v. 3.0](https://opensource.org/licenses/afl-3.0).
  See the accompanying license file for full terms.

- The James Cook University logo (jcucrestcolour.jpg  jculogo.pdf  jculogo.svg) is the trademark of [James Cook University](https://www.jcu.edu.au/).
  The views expressed in this thesis do not represent any official position of James Cook University.

  ![JCU logo](jculogo.svg)


Version Info
------------

This thesis uses a modified form of [semantic versioning](http://semver.org/) suitable for scientific documents.
The version number is of the form v_MAJOR_._MINOR_._PATCH_ where MAJOR, MINOR and PATCH are non-negative integers.

v1.0.0 **will always** point to the final approved version of the thesis with examiner corrections.

The PATCH version is bumped for any changes that
  * Alter the rendered `pdf`, and
  * Affect only formatting, spelling or grammar,
  * **Without** changing the meaning of any passage or affecting the scientific conclusions of the thesis.

The MINOR version is bumped for any changes that
  * Alter the rendered `pdf`, and
  * Alter the meaning or structure of a passage, but
  * **Do not** affect or change any of the scientific conclusions of the thesis.

The MAJOR version is bumped for any changes that
  * Alter the rendered `pdf`, and
  * **Do** affect or change any of the scientific conclusions of the thesis.
  * An exception to this rule (the _only_ exception) is any changes made _before_ the thesis is sent for examination. The MAJOR version number for all pre-examination versions of the thesis shall always be zero.

Changes which do not alter the rendered `pdf` do not bump the version number.

Occasionally pre-release tags of the form v0.1.0-draft.1, v1.0.0-rc1 or the like may be used if this proves useful for something.
Such pre-release tags are lower in precedence than the corresponding release tags, e.g. v0.1.0 and v1.0.0 respectively for those examples.


Building
--------

This thesis is built using [LyX](http://www.lyx.org/) / LaTeX. Supporting code / scripts are written in Python and Mathematica. It has been tested to work with:

- LyX Version 2.2.0 with MiKTeX 2.9
- Python 3.5.1 :: Anaconda 4.0.0 (64-bit)
- Mathematica 10.2.0 for Microsoft Windows (64-bit) (July 28, 2015)
- Git Bash: GNU bash, version 4.3.42(5)-release (x86_64-pc-msys)
- Microsoft Windows 7 (64-bit)

To compile the `pdf` version of the thesis in batch mode run
```sh
$ ./build_thesis.sh
```
This creates a file named something like `thesis-03102016-193203.pdf` with the date stamp in the format DDMMYYYY-HHMMSS.

To make a pure `pdflatex` version of the thesis run
```sh
$ ./build_arxiv_dist.sh
```
This creates a file `thesis-bundle.tar.gz` containing the `pdflatex` source and all dependencies in a for suitable for, e.g., uploading to the physics arXiv.

The `code/` directory contains the programs used to reproduce the scientific content of the thesis from scratch.
This includes all of the figure files.
The `thesis_*_plots.py` programs generate the plots from data stored in the `hartree-results-03102014.npy` file be default.
The `sym_imp_scalars_hartree.py` program runs the computation to produce the data set and checks it against the stored data to guarantee reproducibility.

The scripts must be run inside the `code/` directory.
Doing
```sh
$ cd code
$ python sym_imp_scalars_hartree.py
```
**should** (after a wait) produce the output
```
Wrapping <function hartree_integral at 0x00000000032E9B70> with PrecomputeInterpolate
Took 1.063158 seconds.
Step 0 of 305. temp/Tc = 0.000000
Step 1 of 305. temp/Tc = 0.007662
 < cut out ~300 lines >
 Step 304 of 305. temp/Tc = 1.900825
 *****  Results  *****
 Higgs model: vev = 93 mH = 500 N = 4 lambda = 86.715 CRIT_TEMP = 131.52
 Total time: 147.384346
 *****  Reproducibility Check  *****
 Loading results file: 'hartree-results-03102014.npy'.
 Results are the same shape: True
 Results are nans in all the same places: True
 RMS Residual = 1.276132 / 12522.904213 == 0.000102
```
Thinking of the data set (with nans filtered) as a vector in a high dimensional space, `RMS Residual` is (the Euclidean norm of the error vector / the Euclidean norm of recomputed results).
