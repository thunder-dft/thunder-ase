# fireball_ase_calculator

#### Description
The interface of ASE for FIREBALL.  

[FIREBALL](https://sites.google.com/site/fireballofficialsite/) is a local-orbital ab-initio tight binding implementation of  molecular dynamics. The method allows for the simulation and calculation of very large supercells of thousands of atoms or very long MD  simulations with ease.

[ASE](https://wiki.fysik.dtu.dk/ase/index.html) is a popular python package for atomic structure modeling. ASE provide interfaces with many QM and MM software. This package is the interface for FIREBALL. 



### FIREBALL calculation

#### Excute command:

* **fireball.x**: full version of fireball, calculate both energy and force.
* **lightning.x**: minimum version, only calculate energy.

#### Input files

* `structure.inp`: main input file for calculation, include the atomic structure filename, io options.
* `001.inp` and `001.kpoints`: atomic structure and kpoints files. If no `.kpoints` file provided, gamma-only calculation will be applied.
* `Fdata` directory: the precomputed two center and three center interaction.
* `Fdata.inp`: The input file when creating `Fdata`.

#### structure.inp



#### 001.inp



#### 001.kpoints