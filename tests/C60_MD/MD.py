from ase.build import molecule
from ase import units
from ase.md.verlet import VelocityVerlet as NVE
from ase.md.langevin import Langevin
from ase.md.nvtberendsen import NVTBerendsen
from thunder_ase.fireball import Fireball

atoms = molecule('C60')
Fdata_path = '/home/ren/Programs/Github/thunder-ase/data/Fdata-McWEDA-0.15-3SN.Hs3.75.Cs4.00p4.45.Os3.35p3.80-3SNP.Fes5.30p5.30d4.80'

kwargs = {'iconstraint_l': 1,
          'iwriteout_charges': 1,  # Writing out the charges.
          'taurelax': 5.0,
          'efermi_T': 200.0,
          'ifix_CHARGES': 0,
          'max_scf_iterations_set': 100,
          'scf_tolerance_set': 0.00000001,
          'beta_set': 0.04,
          }

calc = Fireball(command='/home/ren/bin/fireball.3.x',
                Fdata_path=Fdata_path,
                **kwargs)
atoms.set_calculator(calc)

dyn = NVE(atoms, timestep=1.0 * units.fs, trajectory='md.traj', logfile='md.log')
dyn.run(1000)  # take 1000 steps

print("Done!")
