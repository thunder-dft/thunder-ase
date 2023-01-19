### FIREBALL calculation

#### Excute command:

* **fireball.x**: full version of fireball, calculate both energy and force.
* **lightning.x**: minimum version, only calculate energy.

#### Input files

* `structure.inp`: main input file for calculation, include the atomic structure filename, io options.
* `001.inp` and `001.KPOINTS`: atomic structure and kpoints files. If no `.KPOINTS` file provided, gamma-only calculation will be applied.
* `Fdata` directory: the precomputed two center and three center interaction.
* `Fdata.inp`: The input file when creating `Fdata`.

#### structure.inp



#### 001.inp



#### 001.KPOINTS