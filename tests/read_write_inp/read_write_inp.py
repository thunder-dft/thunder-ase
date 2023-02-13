# ### Clean old results
import os
import sys

for f in [
    '001.inp',
    'Fdata.inp',
    'structures.inp',
    '001.KPOINTS',
]:
    if os.path.isfile(f) or os.path.islink(f):
        os.remove(f)

# add the package to import path
sys.path.append("../../")

from thunder_ase.fireball import GenerateFireballInput, Fireball
import ase
import numpy as np

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

atoms = ase.Atoms('N2', pbc=True, cell=np.eye(3) * 10)

# ### test write of GenerateFireballInput
writer = GenerateFireballInput(atoms, **kwargs)
writer.write_options()
writer.write_atoms(pbc=atoms.pbc)
writer.write_kpts()

# This will generate 3 files: '001.inp', 'structures.inp', '001.KPOINTS'

# ### test read of GenerateFireballInput
reader = GenerateFireballInput()
reader.read_kpts('001.KPOINTS')
print(reader.get_kpoints())  # results in a list of kpts and weights

reader.read_options()
print(reader.output_params)
print(reader.options_params)
print(reader.xsfoptions_params)
