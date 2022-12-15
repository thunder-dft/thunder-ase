"""
The interface of ASE for FIREBALL.

* create input files for FIREBALL
* run fireball calculations
* read output


TODO:
* input 大小敏感： 需要去掉敏感性

"""
import os
import subprocess
from typing import Dict, Any
from collections import defaultdict
import ase
import numpy as np
from ase.calculators.calculator import Calculator, FileIOCalculator, CalculatorError, CalculationFailed, all_changes
from ase.dft.kpoints import monkhorst_pack
from ase.geometry import get_distances
import ase.spacegroup


def get_kpts(atoms, size=None, offset=None, reduced=True, **kwargs):
    """
    Get reduced kpoints.
    Reference:
    * https://wiki.fysik.dtu.dk/ase/ase/dft/kpoints.html
    * https://wiki.fysik.dtu.dk/ase/ase/spacegroup/spacegroup.html
    :param reduced:
    :param offset:
    :param atoms:
    :param size:
    :return: np.array with shape (N, 4), coordinates and weights

    TODO: Need check the algorithm from gpaw/kpt_refine.py
    TODO: compare with vasp in several systems
    """
    if 'symprec' in kwargs:
        symprec = kwargs['symprec']
    else:
        symprec = 1e-05  # 注意这是笛卡尔空间的精度
    # generate MP kpoints
    kpoints = monkhorst_pack(size) + np.asarray(offset)
    # get spacegroup from atoms  ase.spacegroup.get_spacegroup
    sg = ase.spacegroup.get_spacegroup(atoms, symprec=symprec)
    # get equivalent sites sg.equivalent_sites
    irreducible_dct = defaultdict(set)
    reducible_set = set()
    for idx, kpt in enumerate(kpoints):
        if idx not in reducible_set:
            sites, kinds = sg.equivalent_sites([kpt])
            # 计算 kpoints 与 sites 之间的距离，注意这是倒易空间
            dist = get_distances(kpoints, sites, cell=np.eye(3), pbc=True)
            # 如果最小的距离 < symprec ，得到对应的 kpoints 的 index
            kpts_green = set(np.argwhere(dist < symprec/atoms.cell.cellpar.min(), axis=0).tolist()[0])
            if len(kpts_green) > 0:
                irreducible_dct[idx] += kpts_green
                reducible_set += kpts_green
                reducible_set += {idx}
    # get kpt weights
    result = []
    for idx, v in irreducible_dct.items():
        kx, ky, kz = kpoints[idx]
        result.append([kx, ky, kz, len(v)])

    return result


class GenerateFireballInput:

    def __int__(self, atoms=None, **kwargs):
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
            'iensemble': 0,
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

    def write_options(self):

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

    def write_kpts(self, size=None, offset=None, reduced=True, **kwargs):
        if offset is None:
            offset = [0., 0., 0.]
        if size is None:
            size = [1, 1, 1]
        offsets = np.asarray(offset)
        sizes = np.asarray(size, dtype=int)

        if len(sizes.shape) == 1:
            sizes = [sizes for _ in self.atoms_lst]
        else:
            assert sizes.shape == (len(self.atoms_lst), 3)

        if len(offsets.shape) == 1:
            offsets = [offsets for _ in self.atoms_lst]
        else:
            assert offsets.shape == (len(self.atoms_lst), 3)

        for isize, ioffset, sname in zip(sizes, offsets, self.sname_lst):
            kpoints = get_kpts(self.atoms, size=isize, offset=ioffset, reduced=reduced, **kwargs)
            with open(sname+'.kpoints', 'w') as f:
                f.write("{}\n".format(len(kpoints)))
                for k in kpoints:
                    f.write("{:8.6f} {:8.6f} {:8.6f} {:8.6f}\n".format(*k))

    def write_input(self):
        self.write_options()
        self.write_atoms(pbc=self.atoms.pbc)
        self.write_kpts()  # TODO: 确定参数


class Fireball(GenerateFireballInput, Calculator):
    name = 'fireball'
    ase_objtype = 'fireball_calculator'  # For JSON storage

    # Environment commands
    env_commands = None

    implemented_properties = [
        'energy',
    ]

    # Can be used later to set defaults
    default_parameters: Dict[str, Any] = {}

    def __init__(self, atoms=None,
                 Fdata_path='Fdata',
                 command='fireball.x',
                 **kwargs):

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
            Fdata_inp_path = os.path.join(os.path.dirname(Fdata_path), "Fdata.inp")
            if not os.path.isfile(os.path.join(os.path.dirname(Fdata_path), "Fdata.inp")):
                # check the existence of Fdata.inp
                raise FileNotFoundError
            else:
                os.symlink(Fdata_inp_path, '.', target_is_directory=False)

        GenerateFireballInput.__init__(self)
        # initialize the ase.calculators.general calculator
        Calculator.__init__(self, atoms=atoms, **kwargs)

    def calculate(self,
                  atoms=None,
                  properties=('energy', ),
                  system_changes=tuple(all_changes)):
        """Do a VASP calculation in the specified directory.

        This will generate the necessary VASP input files, and then
        execute VASP. After execution, the energy, forces. etc. are read
        from the VASP output files.
        """

        if atoms is not None:
            self.atoms = atoms.copy()

        self.write_input()

        errorcode = self._run(command=self.command,
                              directory=self.directory)

        if errorcode:
            raise CalculationFailed(
                '{} in {} returned an error: {:d}'.format(
                    self.name, self.directory, errorcode))

        # TODO: Read results from calculation
        # self.update_atoms(atoms)
        # self.read_results()

    def _run(self, command=None, directory=None):
        """Method to explicitly execute VASP"""
        if command is None:
            command = self.command
        if directory is None:
            directory = self.directory
        errorcode = subprocess.call(command,
                                    shell=True,
                                    cwd=directory)
        return errorcode
