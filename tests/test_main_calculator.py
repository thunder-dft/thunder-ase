import shutil

import pytest
from ase import units, Atoms
from ase.build import molecule
import os
from thunder_ase.fireball import Fireball
import subprocess


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
        kwargs = {'kpt_size': [6, 6, 6]}
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

    def test_md(self):
        """
        Use diamond as test
        :return:
        """
        from ase.cluster import Icosahedron
        from ase.md.verlet import VelocityVerlet as NVE

        max_step = 10  # run 10 steps
        diamond = Icosahedron('C', noshells=3, latticeconstant=2.4)
        diamond.positions[:, 0] = diamond.positions[:, 0] + 10
        Fdata_path = '/home/ren/Fdata/Fdata-Horsfield-0.10-9SN.Hs4.10-9DN.Cs4.35p4.80.Ns3.95p4.40.Os3.35p3.80'
        kwargs = {'ipi': 1}
        fireball = Fireball(command='fireball-ase.9.x', Fdata_path=Fdata_path, **kwargs)
        dyn = NVE(diamond, timestep=1.0 * units.fs, trajectory='md-nve.traj', logfile='md-nve.log')
        fireball.dynamics(dyn, steps=max_step)
        E0 = diamond.get_potential_energy()

        assert E0 == -308.687861

    def test_md_C60(self):
        """
        C60 as test
        Time: 209.54s
        :return:
        """
        from ase.build import molecule
        C60 = molecule('C60')
        from ase.md.verlet import VelocityVerlet as NVE
        max_step = 10  # run 10 steps
        Fdata_path = '/home/ren/Fdata/Fdata-Horsfield-0.10-9SN.Hs4.10-9DN.Cs4.35p4.80.Ns3.95p4.40.Os3.35p3.80'
        kwargs = {'ipi': 1}
        dyn = NVE(C60, timestep=1.0 * units.fs, trajectory='md-nve.traj', logfile='md-nve.log')
        fireball = Fireball(command='fireball-ase.9.x', Fdata_path=Fdata_path, **kwargs)
        fireball.dynamics(dyn, steps=max_step)

    def test_mwfn(self):
        # Example from Preeya using Multiwfn calculate dipole of N2
        atoms = Atoms('N2', [(-0.420, 1.215, 2.018), (-1.561, 1.215, 2.018)])
        Fdata_path = '/home/ren/Fdata/Fdata-Horsfield-0.10-9SN.Hs4.10-9DN.Cs4.35p4.80.Ns3.95p4.40.Os3.35p3.80'
        kwargs = {'iwriteout_cdcoeffs':1}
        fireball = Fireball(command='fireball-ase.9.x', Fdata_path=Fdata_path, **kwargs)
        atoms.calc = fireball
        E0 = atoms.get_potential_energy()
        atoms.calc.write_mwfn()
        # Use Multiwfn to calc dipole
        multiwfn_input = [
            '300',
            '5',
            '0',
            'q',
        ]
        input_name = 'dipole.inp'
        with open(input_name, 'w') as f:
            f.write('\n'.join(multiwfn_input))

        os.system(f'Multiwfn 001.mwfn < {input_name} > mwfn.log')
        command = ["awk", "/Magnitude of dipole moment/ {print $7,$8}", "mwfn.log"]
        result = subprocess.run(command, capture_output=True, text=True)
        assert result.returncode == 0
        assert pytest.approx(float(result.stdout.split()[0]),abs=0.001)==0