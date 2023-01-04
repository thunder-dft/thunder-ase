from fireball_calculator.fireball import GenerateFireballInput, Fireball
import ase


kwargs = {}

atoms = ase.Atoms('N2', pbc=False)

writer = GenerateFireballInput(atoms, **kwargs)


writer.write_options()
writer.write_atoms(pbc=atoms.pbc)
writer.write_kpts()