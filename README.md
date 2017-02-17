Symmetry improvement techniques for non-perturbative quantum field theory
=========================================================================


About
-----

This is the source distribution for Michael J. Brown's PhD thesis entitled "Symmetry improvement techniques for non-perturbative quantum field theory".

Except where otherwise noted, all parts of this thesis are Copyright © 2016 by Michael J. Brown.

Version History
---------------

- [v1.0.0 (Github link)](https://github.com/mikejbrown/thesis/releases/tag/v1.0.0): The final approved version with all examiner corrections. [Download the rendered pdf here](https://github.com/mikejbrown/thesis/releases/download/v1.0.0/michael.j.brown-thesis-v1.0.0.pdf). [![DOI](https://zenodo.org/badge/69083959.svg)](https://zenodo.org/badge/latestdoi/69083959)


Abstract
--------

Quantum field theory is the mathematical language of nature at the deepest levels currently known. Viewed in this light, much of the last century of theoretical physics can be seen as quantum field theory calculations performed in a variety of approximations. However, despite having a very successful standard model of physics and a phenomenally useful perturbation theory for extracting predictions from it, persistent problems remain. Many of these are due to the inability to calculate in practice what is calculable in principle. The difficulty of these problems comes from the importance of a large number of degrees of freedom, or a strong coupling between degrees of freedom, or both. Only in the case of very simple systems or systems close to thermal equilibrium can their properties be determined with confidence. The remaining cases require a resummation, or re-organisation, of perturbation theory. However, ad hoc resummations are problematic because perturbation series are asymptotic in nature and the re-organisation of the terms of an asymptotic series is mathematically questionable. Systematic resummation schemes are required to guarantee consistency with the original non-perturbative theory.

n-particle irreducible effective actions (nPIEAs; n=1,2,3,...) are a proven method of resummation which can be thought of as generalisations of mean field theory which are (a) elegant, (b) general, (c) in principle exact, and (d) have been promoted for their applicability to non-equilibrium situations. These properties make the nPIEAs compelling candidates for filling the gap between theory and experiment anywhere the problem is the inability to compute the behaviour of strongly interacting degrees of freedom or degrees of freedom which are out of thermal equilibrium. Unfortunately, nPIEAs are known to violate the symmetries of a field theory when they are truncated to a form that can be solved in practice. This can lead to qualitatively wrong physical predictions as an artefact of the method of calculation. If one takes the nPIEA predictions at face value one may reject a physically correct model on the basis of an incorrect prediction. Therefore it is critical to gain a better understanding of the symmetry problem in order for nPIEAs to be useful in the future.

This is where this thesis aims to make a contribution. This thesis examines the theory of global symmetries in the nPIEA formalism for quantum field theory with the goal of improving the symmetry properties of solutions obtained from practical approximation schemes. Following pedagogical reviews of quantum field theory and nPIEAs, symmetries are introduced by considering a scalar \phi^{4} field theory with O(N) symmetry and spontaneous symmetry breaking. The symmetries are embodied by Ward identities which relate the various correlation functions of the theory. When truncated, the nPIEAs for n>1 violate these identities with phenomenological consequences which are reviewed. The symmetry improvement method (SI) of [Pilaftsis and Teresi (Nuclear Physics B, 874(2):594-619, 2013)](http://dx.doi.org/10.1016/j.nuclphysb.2013.06.004) addresses this by imposing the 1PI Ward identities directly onto the solutions of nPI equations of motion through the use of Lagrange multipliers. In novel work, SI is extended to the 3PIEA. The SI-3PIEA is then renormalised in the Hartree-Fock and two loop truncations in four dimensions and in the three loop truncation in three dimensions. The Hartree-Fock equations of motion are solved and compared to the unimproved and SI-2PI Hartree-Fock solutions. In both the SI-2PI and SI-3PI methods Goldstone's theorem is satisfied. However, the SI-2PI solution correctly predicts that the phase transition is second order while the SI-3PI solution predicts a weak (weaker than the unimproved 2PIEA) first order transition. Further checks of the formalism are performed: it is shown that the SI-3PIEA obeys the Coleman-Mermin-Wagner theorem and is consistent with unitarity to O(\hbar), but only if the full three loop truncation is kept.

Following this a novel method of soft symmetry improvement (SSI) is introduced. SSI differs from SI in that the Ward identities are imposed softly in the sense of least squared error. A new stiffness parameter controls the strength of the constraint. The motivation for this method is that the singular nature of the constraint in SI leads to pathological behaviour such as the non-existence of solutions in certain truncations and the breakdown of the theory out of equilibrium. In order to regulate the IR behaviour, a box of finite volume is used and the limit of infinite volume taken carefully. Three limiting regimes are found. Two are equivalent to the unimproved 2PIEA and SI-2PIEA respectively. The third is a novel limit which has pathological behaviour including a strongly first order phase transition and a regime where solutions cease to exist for a range of temperatures below the critical temperature. As the stiffness increases this range increases and at a finite value of the stiffness the solution ceases to exist even at zero temperature. This loss of solution is confirmed both analytically and numerically.

Finally, the linear response of a system in equilibrium subject to external perturbations is studied for the first time in the SI-2PIEA formalism. It is shown that, beyond equilibrium, the scheme is inconsistent in that the equations of motion are generically over-constrained. Thus, the SI-2PIEA predicts qualitatively incorrect physics out of equilibrium. This result can be understood as a result of the decoupling of the propagator fluctuation from the field fluctuation in the 2PIEA formalism and the failure of higher order Ward identities. This argument generalises so that, as a result, SI-nPIEA are invalid out of equilibrium, at least within the linear response approximation in generic truncations.

License
-------

Except where otherwise noted, the following licenses apply:

- All non-code assets (LyX/LaTeX source files, accompanying bibtex files, images etc.) are licensed under a [Creative Commons Attribution-ShareAlike 4.0 International License](http://creativecommons.org/licenses/by-sa/4.0/).
  See the accompanying [license file](LICENSE-CC-BY-SA) for full terms.

  ![CC-BY-SA 4.0 badge](cc-by-sa.png)

- All code (Python scripts, Mathematica notebooks, shell scripts etc.) are licensed under the [Academic Free License (“AFL”) v. 3.0](https://opensource.org/licenses/afl-3.0).
  See the accompanying [license file](LICENSE-AFL) for full terms.

- The James Cook University logo (jcucrestcolour.jpg  jculogo.pdf  jculogo.svg) is the trademark of [James Cook University](https://www.jcu.edu.au/).
  The views expressed in this thesis do not represent any official position of James Cook University.

  ![JCU logo](jculogo.svg)


Version Number Info
-------------------

This thesis uses a modified form of [semantic versioning](http://semver.org/) suitable for scientific documents.
The version number is of the form v_MAJOR_._MINOR_._PATCH_ where MAJOR, MINOR and PATCH are non-negative integers.

v1.0.0 **will always** point to the final approved version of the thesis with examiner corrections.

The PATCH version is bumped for any changes that
  * Alter the rendered pdf (apart from version number bumping), and
  * Affect only formatting, spelling or grammar,
  * **Without** changing the meaning of any passage or affecting the scientific conclusions of the thesis.

The MINOR version is bumped for any changes that
  * Alter the rendered pdf (apart from version number bumping), and
  * Alter the meaning or structure of a passage, but
  * **Do not** affect or change any of the scientific conclusions of the thesis.

The MAJOR version is bumped for any changes that
  * Alter the rendered pdf (apart from version number bumping), and
  * **Do** affect or change any of the scientific conclusions of the thesis.
  * An exception to this rule (the _only_ exception) is any changes made _before_ the thesis is sent for examination. The MAJOR version number for all pre-examination versions of the thesis shall always be zero.

Changes which do not alter the rendered pdf (apart from version number bumping) do not bump the version number.

Occasionally pre-release tags of the form v0.1.0-draft.1, v1.0.0-rc1 or the like may be used if this proves useful for something.
Such pre-release tags are lower in precedence than the corresponding release tags, e.g. v0.1.0 and v1.0.0 respectively for those examples.


Building
--------

This thesis is built using [LyX](http://www.lyx.org/) / LaTeX. Supporting code / scripts are written in shell, Python and Mathematica. It has been tested to work with:

- LyX Version 2.2.0 with MiKTeX 2.9
- Python 3.5.1 :: Anaconda 4.0.0 (64-bit)
- Mathematica 10.2.0 for Microsoft Windows (64-bit) (July 28, 2015)
- Git Bash: GNU bash, version 4.3.42(5)-release (x86_64-pc-msys)
- Microsoft Windows 7 (64-bit)

To compile the pdf version of the thesis in batch mode run
```sh
$ ./build_thesis.sh
```
This creates a file named something like `thesis-03102016-193203.pdf` with the date stamp in the format DDMMYYYY-HHMMSS.

To make a pure `pdflatex` version of the thesis run
```sh
$ ./build_arxiv_dist.sh
```
This creates a file `thesis-bundle.tar.gz` containing the `pdflatex` source and all dependencies in a form suitable for, e.g., uploading to the physics arXiv.

The `code/` directory contains the programs used to reproduce the scientific content of the thesis from scratch.
This includes all of the figure files.
The `thesis_*_plots.py` programs generate the plots from data stored in the `hartree-results-03102014.npy` file by default.
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
Thinking of the data set (with nans filtered) as a vector in a high dimensional space, `RMS Residual` is the Euclidean norm of the error vector divided by the Euclidean norm of recomputed results.
