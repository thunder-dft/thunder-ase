[![Upload Python Package](https://github.com/thunder-dft/thunder-ase/actions/workflows/python-publish.yml/badge.svg)](https://github.com/thunder-dft/thunder-ase/actions/workflows/python-publish.yml)

# Thunder ASE

### Description
The interface of ASE for FIREBALL.  

[FIREBALL](https://sites.google.com/site/fireballofficialsite/) is a local-orbital ab-initio tight binding implementation of  molecular dynamics. The method allows for the simulation and calculation of very large supercells of thousands of atoms or very long MD  simulations with ease.

[ASE](https://wiki.fysik.dtu.dk/ase/index.html) is a popular python package for atomic structure modeling. ASE provide interfaces with many QM and MM software. This package is the interface for FIREBALL. 

### Installation

`pip install thunder-ase`

#### Requirements

* ase

## Roadmap

* v0.1: run basic fireball calculation and read basic result from it.
  * v0.1.1: jupyter-notebook for examples
  * v0.1.2: multiple atoms for one calculator
  * v0.1.3: read fireball input to construct calculator
* v0.2: band structure calculation.
  * v0.2.1: DOS calculation.
* v0.3.1: fit basis to gaussian, write orbitals to mwfn.
  * v0.3.2: interface to multiwfn.
* v0.4: MD calculation.
* v0.5: Fdata management.
* v0.6: interactive running.
