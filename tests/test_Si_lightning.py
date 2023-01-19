from fireball_calculator.fireball import Fireball
import ase
from ase.io.trajectory import Trajectory
import numpy as np
from ase.io import read
from ase.units import kJ
from ase.eos import EquationOfState


cell = np.array([[2.715000, 2.715000, 0.000000],
                 [2.715000, 0.000000, 2.715000],
                 [0.000000, 2.715000, 2.715000]])
positions = np.array([[0.0000000, 0.0000000, 0.0000000],
                      [1.3575000, 1.3575000, 1.3575000]])
kwargs = {'kpt_size': [2, 2, 2],
          'iwriteout_ME_SandH': 0,
          'iwriteout_density': 0,
          'iwriteout_cdcoeffs': 0,
          'iwriteout_charges': 1,
          'iwriteout_energies': 0,
          'iwriteout_populations': 0,
          'iwriteout_forces': 1,
          'iwriteout_neighbors': 0,
          'iwriteout_dos': 0,
          'iwriteout_abs': 0,
          'iwriteout_ewf': 0,
          'nstepi': 1,
          'nstepf': 100,
          'iquench': -3,
          'T_initial': 0.0,
          'T_final': 0.0,
          'T_want': 0.0,
          'taurelax': 5.0,
          'efermi_T': 200.0,
          'dt': 1.00,
          'iensemble': 0,
          'iconstraint_rcm': 1,
          'iconstraint_vcm': 1,
          'iconstraint_L': 0,
          'iconstraint_KE': 1,
          'ifix_neighbors': 0,
          'ifix_CHARGES': 0,
          'max_scf_iterations_set': 100,
          'scf_tolerance_set': 0.00000001,
          'beta_set': 0.04,
          'Ecut_set': 200.0,
          'rho_surface_min': 0.0005,
          'rho_surface_max': 0.01000000,
          }

cell_factors = np.linspace(0.8, 1.2, 8)
traj = Trajectory('Si.traj', 'w')
for cf in cell_factors:
    atoms = ase.Atoms(numbers=[14, 14],
                      cell=cell,
                      pbc=True,
                      positions=positions,
                      )
    atoms.set_cell(cell=cell*cf, scale_atoms=True)
    calc = Fireball(command='lightning.3.x', **kwargs)
    atoms.set_calculator(calc)
    e0 = atoms.get_potential_energy()
    traj.write(atoms)
    print("The energy for cell factor {:.3f} is {:.3f}".format(cf, e0))

configs = read('Si.traj', index=':')
volumes = [si.get_volume() for si in configs]
energies = [si.get_potential_energy() for si in configs]
eos = EquationOfState(volumes, energies)
v0, e0, B = eos.fit()
print(B / kJ * 1.0e24, 'GPa')
eos.plot('Si-eos.png')

print("Done!!")
