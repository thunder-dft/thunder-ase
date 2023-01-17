from fireball_calculator.fireball import GenerateFireballInput, Fireball
import ase
import numpy as np

kwargs = {}

atoms = ase.Atoms('N2', pbc=True, cell=np.eye(3)*10)

writer = GenerateFireballInput(atoms, **kwargs)

writer.write_options()
writer.write_atoms(pbc=atoms.pbc)
writer.write_kpts()
