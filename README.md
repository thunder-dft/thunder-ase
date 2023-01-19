# fireball_ase_calculator

### Description
The interface of ASE for FIREBALL.  

[FIREBALL](https://sites.google.com/site/fireballofficialsite/) is a local-orbital ab-initio tight binding implementation of  molecular dynamics. The method allows for the simulation and calculation of very large supercells of thousands of atoms or very long MD  simulations with ease.

[ASE](https://wiki.fysik.dtu.dk/ase/index.html) is a popular python package for atomic structure modeling. ASE provide interfaces with many QM and MM software. This package is the interface for FIREBALL. 

### Installation

`pip install git+https://gitee.com/pjren/fireball_ase_calculator.git`

#### Requirements
* ase

## Roadmap

* v0.1: run basic fireball calculation and read basic result from it.
  * v0.1.1: jupyter-notebook for tests
* v0.2: band structure calculation.
* v0.3: MD calculation.
* v0.4: update Fdata management.
* v0.5: interactive running.