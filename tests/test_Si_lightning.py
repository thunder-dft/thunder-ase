from fireball_calculator.fireball import GenerateFireballInput, Fireball

import ase
import numpy as np

atoms = ase.Atoms(numbers=[14, 14],
                  cell=[
                      [2.715000, 2.715000, 0.000000],
                      [2.715000, 0.000000, 2.715000],
                      [0.000000, 2.715000, 2.715000]],
                  positions=[
                      [3.141593, 0.367879, 1.414214],
                      [4.499093, 1.725379, 2.771714]],
                  pbc=True,
                  )
# Write input
kwargs = {'kpt_size': [3, 3, 3]}
writer = GenerateFireballInput(atoms, **kwargs)

writer.write_options()
writer.write_atoms(pbc=atoms.pbc)
writer.write_kpts()

print("Done")
