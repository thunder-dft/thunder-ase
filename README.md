[![Upload Python Package](https://github.com/thunder-dft/thunder-ase/actions/workflows/python-publish.yml/badge.svg)](https://github.com/thunder-dft/thunder-ase/actions/workflows/python-publish.yml)

# Thunder ASE

### Description
The interface of ASE for FIREBALL.

[FIREBALL](https://fireball-dft.org) is an open source and high efficiency DFT software based on local-orbital ab-initio method, allows for the simulation very large supercells of thousands of atoms. [ASE](https://wiki.fysik.dtu.dk/ase/index.html) is a popular python package for atomic structure modeling. ASE provide interfaces with many QM and MM software. This `thunder-ase` package is not only the ASE interface for FIREBALL, but also provides many other useful functions. 

![thunder-ase-framework](./Docs/img/thunder-ase-framework.jpeg)

### Installation

`pip install -U thunder-ase`

#### Requirements

* ase

## Roadmap

* v0.1: run basic fireball calculation and read basic result from it.
  * jupyter-notebook for examples
  * multiple atoms for one calculator
  * read fireball input to construct calculator
* v0.2: band structure calculation.
  * DOS calculation.
* v0.3: fit basis to gaussian basis set, write orbitals to mwfn.
  * interface to multiwfn.
* v0.4: Read charge and population.
* v0.5: interactive running.
* v0.6: Fdata management.
