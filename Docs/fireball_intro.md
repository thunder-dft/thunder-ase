### FIREBALL calculation

#### Excute command:

* **fireball.x**: full version of fireball, calculate both energy and force.
* **lightning.x**: minimum version, only calculate energy.

#### Input files

* `structure.inp`: main input file for calculation, include the atomic structure filename, io options.
* `001.inp` and `001.KPOINTS`: atomic structure and kpoints files. If no `.KPOINTS` file provided, gamma-only calculation will be applied.
* `Fdata` directory: the precomputed two center and three center interaction.
* `Fdata.inp`: The input file when creating `Fdata`.

### structure.inp

#### Basic Parameters of FIREBALL

* `qstate` : The charge state of the system. Default is 0.
* `efermi_T`: Smearing temperature at Fermi level. Default is 100.0, which is about 0.0086 eV. For metal system, larger value is recommended.
* `max_scf_iterations_set`: SCF steps maximum. Default is 50.
* `scf_tolerance_set`: SCF tolerance of charge. Default is 1.0E-6. Smaller value is recommended for higher precision, e.g. 1.0E-8.
* `beta_set`: Broyden’s mix factor of charges. Default is 0.08.
* `Ecut_set`: To control mesh grid density. Default is 200.0.
* `iwriteout_charges`: Write out charges. Default is  0 (False).

#### Advanced Parameters of FIREBALL

* `iconstraint_rcm`: Whether shifts molecule to center of mass (COM) before the simulation. Default is 1 (True).
* `ifix_neighbors`: Fix the neighbor list of the system to initial structure. Default is  0 (False).
* `ifix_CHARGES`: Fix the charges of the system to initial structure. Default is  0 (False). This is useful sometimes for very unreasonable initial structures.
* `iwriteout_me_sandh`: Write out overlap (S) and Hamiltonian (H) matrix. Default is 0 (False).
* `iwriteout_cdcoeffs`: Write out orbital coefficients. Default is 0 (False).
* `iwriteout_density`: Write out density matrix. Default is  0 (False).
* `iwriteout_energies`: Write out all energy terms. Default is  0 (False).
* `iwriteout_populations`: Measure the localization by the entropic quantity W. Default is  0 (False). W is the number of accessible atoms for the given electronic energy and which describes the spatial extent of the electronic state. See more in Refs. H. Wang and J. P. Lewis, J. Phys.: Condens. Matter 18, 421–434 (2005) and  H. Wang and J. P. Lewis, J. Phys.: Condens. Matter 17, L209–L213 (2005).
* `iwriteout_forces`: Write out all force terms. Default is  0 (False).
* `iwriteout_neighbors`: Write out all neighbors. Default is  0 (False).
* `iwriteout_abs`: Write out  absorption spectra. Default is  0 (False).
* `iwriteout_ewf`: Write out  wavefuncion. Default is  0 (False).
* `iwriteout_dos`: Write out  density of states. Default is  0 (False).

#### Deprecated Parameters of FIREBALL

The following parameters are deprecated due to the usage of thunder-ase, in the future, they will be removed from FIREBALL. Most parameters are for MD simulations.

* `nstepi`: Initial step number. Default is 1.
* `nstepf`: Final step number. Default is 1. `nstepf - nstepi` is the required MD steps.
* `iconstraint_vcm`: Constraining the velocities about the COM, it means whether fix  the COM during simulation. Default is 1 (True).
* `iconstraint_L`: Angular momentum constraint of the whole system. Default is  0 (False).
* `iconstraint_KE`: Kinetic energy constraint of the whole system. Default is  1 (True).
* `T_initial`: Initial temperature of MD simulation. Default is 300.0. 
* `T_final`: Final temperature of MD simulation. Default is 0.. 
* `iquench`: Quench method. Default is 0. 
  * `iquench = n`: quench the velocities on every n’th step; 
  * `iquench = -1`: quench whenever the instantaneous temperature `T_instantaneous` is lower than the instantaneous temperature on the previous step `T_previous`;
  * `iquench = -2`: annealing simulation.  `T_want` is the desired annealing temperature and `taurelax` is for specifies how rapidly of quench;
  * `iquench = -3`: Coordinate power quench, quench velocities one coordinate at a time.
* `T_want`: Annealing temperature to reach.  Default is 300.0 K.
* `taurelax`: Quench rate. Default is 5.0 fs.
* `dt`: Time step of MD simulation. Default is 0.25 fs.
* `iensemble`: Thermodynamic ensemble. Default is 0.
  * `iensemble = 1`: Constant temperature ensemble
* `rho_surface_min` and `rho_surface_max`: Electron density minimum and maximum for XCrySDen Structure File. Default are 0.5E-3 and 0.1. 







#### 001.inp



#### 001.KPOINTS