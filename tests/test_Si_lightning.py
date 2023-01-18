from fireball_calculator.fireball import GenerateFireballInput, Fireball

import ase
import numpy as np

atoms = ase.Atoms(numbers=[14, 14],
                  cell=[
                      [2.715000, 2.715000, 0.000000],
                      [2.715000, 0.000000, 2.715000],
                      [0.000000, 2.715000, 2.715000]],
                  positions=[
                      [0.0000000, 0.0000000, 0.0000000],
                      [1.3575000, 1.3575000, 1.3575000]],
                  pbc=True,
                  )
# Write input
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
writer = GenerateFireballInput(atoms, **kwargs)

writer.write_options()
writer.write_atoms(pbc=atoms.pbc)
writer.write_kpts()

print("Done")
