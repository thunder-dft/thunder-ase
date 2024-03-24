import ase.io
from ase.optimize.bfgs import BFGS
from thunder_ase.fireball import Fireball
from ase.build import molecule
from ase.md.verlet import VelocityVerlet as NVE
from ase import units

atoms = ase.io.read('/home/ren/data/spectroscopy_dev/STM/crystal_data/C/POSCAR')

Fdata_path = '/home/ren/Fdata/Fdata-Horsfield-0.10-9SN.Hs4.10-9DN.Cs4.35p4.80.Ns3.95p4.40.Os3.35p3.80'
kpt_size = (8, 8, 1)
#max_step = 10
kwargs = {
    'kpt_size': kpt_size,
    'kpt_reduced': True,
    'kpt_gamma': True,
    'ipi': 1,  # use i-pi socket to save reading Fdata time
#    'nstepf': max_step + 1,  # max step should set as this
}

fireball = Fireball(command='fireball-ase.9.x', Fdata_path=Fdata_path, **kwargs)
#dyn = BFGS(atoms, trajectory='test_dyn_run.traj')
#fireball.dynamics(dyn, fmax=0.1)

max_step = 10
dyn = NVE(atoms, timestep=1.0 * units.fs, trajectory='md-nve.traj', logfile='md-nve.log')
fireball.dynamics(dyn, steps=max_step)

#socket_run(atoms, dyn=dyn, calculator=fireball, fmax=0.1, max_step=max_step)
print('Done')
