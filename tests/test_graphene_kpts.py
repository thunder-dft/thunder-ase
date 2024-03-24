import ase.io

from thunder_ase.fireball import Fireball

atoms = ase.io.read('/home/ren/data/spectroscopy_dev/STM/crystal_data/C/POSCAR')

Fdata_path = '/home/ren/Fdata/Fdata-Horsfield-0.10-9SN.Hs4.10-9DN.Cs4.35p4.80.Ns3.95p4.40.Os3.35p3.80'
kpt_size = (8, 8, 1)
kwargs = {
    'kpt_size': kpt_size,
    'kpt_reduced': True,
    'kpt_gamma': True,
}

fireball = Fireball(command='fireball-ase.9.x', Fdata_path=Fdata_path, **kwargs)

atoms.calc = fireball
print(atoms.get_potential_energy())

kwargs = {
    'kpt_size': kpt_size,
    'kpt_reduced': False,
    'kpt_gamma': True,
}

fireball = Fireball(command='fireball-ase.9.x', Fdata_path=Fdata_path, **kwargs)

atoms.calc = fireball
print(atoms.get_potential_energy())

print('Done')
