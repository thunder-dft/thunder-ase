import shutil

import pytest
from ase.build import molecule
import os
from thunder_ase.fireball import Fireball


class TestFireball:
    def setup_method(self):
        # create tmp directory and get in
        if os.path.isdir('tmp'):
            shutil.rmtree('tmp')
        os.mkdir('tmp')
        os.chdir('tmp')

    def teardown_method(self):
        # clean all result files
        os.chdir('..')
        shutil.rmtree('tmp')

    def test_molecule_sp(self):
        atoms = molecule('C6H6')
        Fdata_path = '/home/ren/Fdata/Fdata-Horsfield-0.10-9SN.Hs4.10-9DN.Cs4.35p4.80.Ns3.95p4.40.Os3.35p3.80'
        kwargs = {}
        fireball = Fireball(command='fireball-ase.9.x', Fdata_path=Fdata_path, **kwargs)
        atoms.calc = fireball
        E0 = atoms.get_potential_energy()

        assert E0 == -1017.283695

    def test_crystal_sp(self):
        from ase.build import bulk
        diamond = bulk('C')
        Fdata_path = '/home/ren/Fdata/Fdata-Horsfield-0.10-9SN.Hs4.10-9DN.Cs4.35p4.80.Ns3.95p4.40.Os3.35p3.80'
        kwargs = {'kpt_size': [6, 6, 6],}
        fireball = Fireball(command='fireball-ase.9.x', Fdata_path=Fdata_path, **kwargs)
        diamond.calc = fireball
        E0 = diamond.get_potential_energy()

        assert E0 == -308.687861

    def test_minimize(self):
        atoms = molecule('C6H6')
        Fdata_path = '/home/ren/Fdata/Fdata-Horsfield-0.10-9SN.Hs4.10-9DN.Cs4.35p4.80.Ns3.95p4.40.Os3.35p3.80'
        kwargs = {'ipi': 1}
        fireball = Fireball(command='fireball-ase.9.x', Fdata_path=Fdata_path, **kwargs)
        atoms.calc = fireball
        fireball.minimize(atoms, method='MDMin', fmax=0.1, trajectory='minimize.traj',
                          logfile='minimize.log')
        Eopt = atoms.get_potential_energy()

        assert Eopt == -1017.7185886746582

