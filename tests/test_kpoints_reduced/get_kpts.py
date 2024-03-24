from thunder_ase.fireball import get_kpts
import ase.io

atoms = ase.io.read('POSCAR')
kpts = get_kpts(atoms, size=(8, 8, 8), offset=None, reduced=True, gamma=False)

print(len(kpts))

print('Done')
