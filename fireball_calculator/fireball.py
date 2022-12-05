"""
The interface of ASE for FIREBALL.

* create input files for FIREBALL
* run fireball calculations
* read output


TODO:
* 确认大小写是否敏感
* Kpoints 采样，放到哪一部分？

"""
import os
from typing import Dict, Any

import ase
import numpy as np
from ase.calculators.calculator import Calculator, FileIOCalculator
from ase.dft.kpoints import monkhorst_pack


class GenerateFireballInput:

    def __int__(self, atoms=None):
        if atoms is None or len(atoms) == 0:
            raise ValueError
        self.atoms = atoms
        if type(self.atoms) in (ase.Atoms,):
            self.atoms_lst = [atoms]
        elif type(self.atoms) in (tuple, list):
            self.atoms_lst = atoms
        else:
            raise NotImplementedError

        self.sname_lst = ["{:03d}".format(idx+1) for idx in range(len(self.atoms_lst))]

        self.options_params = {
            'nstepi': 1,
            'nstepf': 1,
            'iquench': 0,
            'T_initial': 300.0,
            'T_final': 0.0,
            'T_want': 300.0,
            'taurelx': 5.0,
            'efermi_T': 100.0,
            'dt': 0.25,  # fs
            'iensembel': 0,
            'iconstraint_rcm': 1,
            'iconstraint_vcm': 1,
            'iconstraint_L': 0,
            'iconstraint_KE': 1,
            'ifix_neighbors': 0,
            'ifix_CHARGES': 1,
            'max_scf_iterations_set': 50,
            'scf_tolerance_set': 0.00000001,
            'beta_set': 0.08,
            'Ecut_set': 200.0,
        }

        self.output_params = {
            'iwriteout_ME_SandH': 0,
            'iwriteout_density': 0,
            'iwriteout_cdcoeffs': 0,
            'iwriteout_charges': 0,
            'iwriteout_energies': 0,
            'iwriteout_populations': 0,
            'iwriteout_forces': 0,
            'iwriteout_neighbors': 0,
            'iwriteout_dos': 0,
            'iwriteout_abs': 0,
            'iwriteout_ewf': 0,
        }

        self.xsfoptions_params = {
            'rho_surface_min': 0.0005,
            'rho_surface_max': 0.1,
        }

        self.all_params = self.options_params | self.output_params | self.xsfoptions_params

    def write_structure(self):

        with open('structure.inp', 'w') as f:
            f.write("{}\n".format(len(self.atoms_lst)))
            for sname in self.sname_lst:
                f.write("{}.inp\n".format(sname))

            f.write("! Write out options - leave this line as comment line or null\n")

            f.write("&OUTPUT\n")
            for k, v in self.output_params.items():
                f.write("{} = {}\n".format(k, v))
            f.write("&END\n")

            f.write("&OPTIONS\n")
            for k, v in self.options_params.items():
                f.write("{} = {}\n".format(k, v))
            f.write("&END\n")

            f.write("&XSFOPTIONS\n")
            for k, v in self.xsfoptions_params.items():
                f.write("{} = {}\n".format(k, v))
            f.write("&END\n")

    def write_atoms(self, pbc=None):

        for sname, iatoms in zip(self.sname_lst,self.atoms_lst):
            with open(sname+".inp", 'w') as f:
                if pbc is None:
                    ipbc = 1 if np.any(iatoms.pbc) else 0
                else:
                    ipbc = 1 if pbc else 0
                f.write("{:3d}{:12d}\n".format(len(iatoms), ipbc))

                if ipbc == 1:
                    cell = iatoms.cell[:]
                else:
                    cell = np.eye(3) * 999.0
                for row in cell:
                    f.write("{:11.6f} {:11.6f} {:11.6f}\n".format(*row))

                for num, xyz in zip(iatoms.numbers, iatoms.positions):
                    x, y, z = xyz
                    f.write("{:3d} {:11.6f} {:11.6f} {:11.6f}\n".format(num, x, y, z))

    def write_kpts(self, kpts=None, offset=None):
        if offset is None:
            offset = [0., 0., 0.]
        if kpts is None:
            kpts = [1, 1, 1]
        offset = np.asarray(offset)
        kpts = np.asarray(kpts, dtype=int)

        if len(kpts.shape) == 1:
            kpts = [kpts for _ in self.atoms_lst]
        else:
            assert kpts.shape == (len(self.atoms_lst), 3)

        if len(offset.shape) == 1:
            offset = [offset for _ in self.atoms_lst]
        else:
            assert offset.shape == (len(self.atoms_lst), 3)

        for ikpt, ioffset, sname in zip(kpts, offset, self.sname_lst):
            kpoints = monkhorst_pack(ikpt) + ioffset
            weight = 1.0
            with open(sname+'.kpoints', 'w') as f:
                f.write("{}\n".format(len(kpoints)))
                for k in kpoints:
                    f.write("{:8.6f} {:8.6f} {:8.6f} {:8.6f}\n".format(*k, weight))




class Fireball(GenerateFireballInput, Calculator):
    name = 'fireball'
    ase_objtype = 'fireball_calculator'  # For JSON storage

    # Environment commands
    env_commands = None

    implemented_properties = [
        'energy', 'forces',
    ]

    # Can be used later to set defaults
    default_parameters: Dict[str, Any] = {}

    def __init__(self, atoms=None,
                 Fdata_path='Fdata',
                 command='fireball.x',
                 check_version=False, **kwargs):

        self.__name__ = 'fireball'

        self.Fdata_path = Fdata_path
        self.command = command
        # if Fdata not in current directory, make a symbolic link
        if not os.path.isdir('./Fdata'):
            if not os.path.isdir(Fdata_path):
                # check the existence of Fdata directory
                raise FileNotFoundError
            else:
                os.symlink(Fdata_path, '.', target_is_directory=True)
        if not os.path.isfile("Fdata.inp"):
            Fdata_inp_path = os.join(os.path.dirname(Fdata_path), "Fdata.inp")
            if not os.path.isfile(os.join(os.path.dirname(Fdata_path), "Fdata.inp")):
                # check the existence of Fdata.inp
                raise FileNotFoundError
            else:
                os.symlink(Fdata_inp_path, '.', target_is_directory=False)

        GenerateFireballInput.__init__(self)
        # initialize the ase.calculators.general calculator
        Calculator.__init__(self, atoms=atoms)
