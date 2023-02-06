"""
The interface of ASE for FIREBALL.
* create input files for FIREBALL
* run fireball calculations
* read output
"""
import os
import subprocess
from typing import Dict, Any
import ase
import numpy as np
from ase.calculators.calculator import Calculator, CalculationFailed, all_changes
from ase.dft.kpoints import monkhorst_pack, kpoint_convert
from ase.geometry import get_distances
import ase.spacegroup
from ase.io import jsonio


def get_kpts(atoms, size=None, offset=None, reduced=False, **kwargs):
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
    # generate MP kpoints
    kpoints = monkhorst_pack(size) + np.asarray(offset)
    kpoints = kpoint_convert(cell_cv=atoms.cell, skpts_kc=kpoints)  # convert from scaled and cartesian coordinates
    nkpt = len(kpoints)
    if not reduced or nkpt <= 1:
        kpoints_weight = []
        for k in kpoints:
            kx, ky, kz = k
            kpoints_weight.append([kx, ky, kz, 1.0/nkpt])
        return kpoints_weight

    if 'symprec' in kwargs:
        symprec = kwargs['symprec']
    else:
        symprec = 1e-05  # 注意这是笛卡尔空间的精度
    try:
        # get spacegroup from atoms ase.spacegroup.get_spacegroup
        sg = ase.spacegroup.get_spacegroup(atoms, symprec=symprec)
    except RuntimeError:
        print("Didn't find symmetry. No reducing is needed.")
        kpoints_weight = []
        for k in kpoints:
            kx, ky, kz = k
            kpoints_weight.append([kx, ky, kz, 1.0/nkpt])
        return kpoints_weight

    kpts_extend = []  # all the symmetry points for each kpts
    for kpt in kpoints:
        sites, kinds = sg.equivalent_sites([kpt])
        kpts_extend.append(sites)

    irreducible_dct = dict()
    kpoints_lst = kpoints.tolist()
    while len(kpoints_lst) > 0:
        kpt = kpoints_lst.pop(0)
        kpt_extend = kpts_extend.pop(0)
        key = tuple(kpt)
        irreducible_dct[key] = [kpt]
        # calc the distance between kpt_extend and kpts_extend
        if len(kpoints_lst) > 0:
            reduced_idx = []
            for idx, other_extend in enumerate(kpts_extend[:]):
                dist_array, dist = get_distances(kpt_extend, other_extend, cell=atoms.cell, pbc=True)
                if np.any(dist < symprec):
                    reduced_idx.append(idx)
            irreducible_dct[key] += [kpoints_lst[idx] for idx in reduced_idx]
            # update kpoints_lst and kpts_extend by removing the reduced_idx
            kpoints_lst = [v for idx, v in enumerate(kpoints_lst) if idx not in reduced_idx]
            kpts_extend = [v for idx, v in enumerate(kpts_extend) if idx not in reduced_idx]

    kpoints_weight = [[k[0], k[1], k[2], len(v) / nkpt] for k, v in irreducible_dct.items()]
    return kpoints_weight


options_params = {
    'nstepi': {'type': (int,), 'name': 'nstepi', 'default': 1},
    'nstepf': {'type': (int,), 'name': 'nstepf', 'default': 1},
    'iquench': {'type': (int,), 'name': 'iquench', 'default': 0},
    't_initial': {'type': (int, float), 'name': 'T_initial', 'default': 300.0},
    't_final': {'type': (int, float), 'name': 'T_final', 'default': 0.0},
    't_want': {'type': (int, float), 'name': 'T_want', 'default': 300.0},
    'taurelax': {'type': (int, float), 'name': 'taurelax', 'default': 5.0},
    'efermi_t': {'type': (int, float), 'name': 'efermi_T', 'default': 100.0},
    'dt': {'type': (int, float), 'name': 'dt', 'default': 0.25},  # fs
    'iensemble': {'type': (int,), 'name': 'iensemble', 'default': 0},
    'iconstraint_rcm': {'type': (int,), 'name': 'iconstraint_rcm', 'default': 1},
    'iconstraint_vcm': {'type': (int,), 'name': 'iconstraint_vcm', 'default': 1},
    'iconstraint_l': {'type': (int,), 'name': 'iconstraint_L', 'default': 0},
    'iconstraint_ke': {'type': (int,), 'name': 'iconstraint_KE', 'default': 1},
    'ifix_neighbors': {'type': (int,), 'name': 'ifix_neighbors', 'default': 0},
    'ifix_charges': {'type': (int,), 'name': 'ifix_CHARGES', 'default': 1},
    'max_scf_iterations_set': {'type': (int,), 'name': 'max_scf_iterations_set', 'default': 50},
    'scf_tolerance_set': {'type': (int, float), 'name': 'scf_tolerance_set', 'default': 0.00000001},
    'beta_set': {'type': (int, float), 'name': 'beta_set', 'default': 0.08},
    'ecut_set': {'type': (int, float), 'name': 'Ecut_set', 'default': 200.0},
}

output_params = {
    'iwriteout_me_sandh': {'type': (int,), 'name': 'iwriteout_ME_SandH', 'default': 0},
    'iwriteout_density': {'type': (int,), 'name': 'iwriteout_density', 'default': 0},
    'iwriteout_cdcoeffs': {'type': (int,), 'name': 'iwriteout_cdcoeffs', 'default': 0},
    'iwriteout_charges': {'type': (int,), 'name': 'iwriteout_charges', 'default': 0},
    'iwriteout_energies': {'type': (int,), 'name': 'iwriteout_energies', 'default': 0},
    'iwriteout_populations': {'type': (int,), 'name': 'iwriteout_populations', 'default': 0},
    'iwriteout_forces': {'type': (int,), 'name': 'iwriteout_forces', 'default': 0},
    'iwriteout_neighbors': {'type': (int,), 'name': 'iwriteout_neighbors', 'default': 0},
    'iwriteout_dos': {'type': (int,), 'name': 'iwriteout_dos', 'default': 0},
    'iwriteout_abs': {'type': (int,), 'name': 'iwriteout_abs', 'default': 0},
    'iwriteout_ewf': {'type': (int,), 'name': 'iwriteout_ewf', 'default': 0},
}

xsfoptions_params = {
    'rho_surface_min': {'type': (int, float), 'name': 'rho_surface_min', 'default': 0.0005},
    'rho_surface_max': {'type': (int, float), 'name': 'rho_surface_max', 'default': 0.1},
}

calc_params = {
    'kpt_size': {'type': (list, np.array), 'name': 'kpt_size', 'default': [1, 1, 1]},
    'kpt_offset': {'type': (list, np.array), 'name': 'kpt_offset', 'default': [0., 0., 0.]},
    'kpt_interval': {'type': (list, np.array, float, int), 'name': 'kpt_interval', 'default': None},
}

fireball_params = options_params | output_params | xsfoptions_params | calc_params


def write_params(dct, f):
    for k, v in dct.items():
        kname = fireball_params[k]['name']
        f.write("{} = {}\n".format(kname, v))


class GenerateFireballInput:
    def __init__(self, atoms=None, **kwargs):
        self.atoms = atoms
        self.sname = '001'  # TODO: name for input file

        self.output_params = {}
        self.options_params = {}
        self.xsfoptions_params = {}

        self.kpt_size = None
        self.kpt_interval = None
        self.kpt_offset = None

        self.check_input(kwargs)

    def check_input(self, kwargs):
        for key, v in kwargs.items():
            k = key.lower()
            if k not in fireball_params:
                print("The option {} not supported!".format(k))
                raise KeyError
            if type(v) not in fireball_params[k]['type']:
                print("The type of {} should be {}".format(k, ' or '.join(fireball_params[k]['type'])))
                raise TypeError
            if k in output_params:
                self.output_params[k] = v
            elif k in options_params:
                self.options_params[k] = v
            elif k in xsfoptions_params:
                self.xsfoptions_params[k] = v
            elif k == 'kpt_size':
                self.kpt_size = v
            elif k == 'kpt_offset':
                self.kpt_offset = v
            elif k == 'kpt_interval':
                self.kpt_interval = v
        return

    def write_options(self):

        with open('structures.inp', 'w') as f:
            f.write("1\n")
            f.write("{}.inp\n".format(self.sname))

            f.write("! Write out options - leave this line as comment line or null\n")

            f.write("&OUTPUT\n")
            write_params(self.output_params, f)
            f.write("&END\n")

            f.write("&OPTIONS\n")
            write_params(self.options_params, f)
            f.write("&END\n")

            f.write("&XSFOPTIONS\n")
            write_params(self.xsfoptions_params, f)
            f.write("&END\n")

    def write_atoms(self, pbc=None):

        with open(self.sname + ".inp", 'w') as f:
            if pbc is None:
                ipbc = 0 if np.any(self.atoms.pbc) else 1
            else:
                ipbc = 0 if np.any(pbc) else 1
            f.write("{:3d}{:12d}\n".format(len(self.atoms), ipbc))
            if ipbc == 0:
                cell = self.atoms.cell[:]
            else:
                cell = np.eye(3) * 999.0
            for row in cell:
                f.write("{:11.6f} {:11.6f} {:11.6f}\n".format(*row))
            for num, xyz in zip(self.atoms.numbers, self.atoms.positions):
                x, y, z = xyz
                f.write("{:3d} {:11.6f} {:11.6f} {:11.6f}\n".format(num, x, y, z))

    def write_kpts(self, reduced=True, **kwargs):
        """
        :param reduced:
        :param kwargs:
        :return:
        """
        if not np.all(self.atoms.pbc):
            return

        if self.kpt_offset is None:
            offset = [0., 0., 0.]
        else:
            offset = self.kpt_offset
        if self.kpt_size is None:
            if self.kpt_interval is None:
                size = [1, 1, 1]
            else:
                size = np.ceil(1 / self.atoms.cell.cellpar()[0:3] / self.kpt_interval)
        else:
            size = self.kpt_size

        kpoints = get_kpts(self.atoms, size=size, offset=offset, reduced=reduced, **kwargs)
        with open(self.sname + '.KPOINTS', 'w') as f:
            f.write("{}\n".format(len(kpoints)))
            for k in kpoints:
                f.write("{:8.6f} {:8.6f} {:8.6f} {:8.6f}\n".format(*k))

    def write_input(self):
        self.write_options()
        self.write_atoms(pbc=self.atoms.pbc)
        self.write_kpts()


class Fireball(GenerateFireballInput, Calculator):
    name = 'fireball'
    ase_objtype = 'fireball_calculator'  # For JSON storage

    # Environment commands
    env_commands = None

    implemented_properties = [
        'energy',
        'forces',
        'fermi',
    ]

    # Can be used later to set defaults
    default_parameters: Dict[str, Any] = {}

    def __init__(self, atoms=None,
                 Fdata_path='Fdata',
                 command='fireball.x',
                 **kwargs):
        self.__name__ = 'fireball'
        self._atoms = None
        GenerateFireballInput.__init__(self, **kwargs)  # TODO: kwargs 使用set 函数进行管理，参考vasp
        self.Fdata_path = Fdata_path
        self.command = command
        # if Fdata not in current directory, make a symbolic link
        if not os.path.isdir('Fdata'):
            if not os.path.isdir(Fdata_path):
                # check the existence of Fdata directory
                raise FileNotFoundError
            else:
                Fdata_dir_path = os.path.join(os.path.dirname(Fdata_path), "Fdata")
                os.symlink(Fdata_dir_path, 'Fdata', target_is_directory=True)
        if not os.path.isfile("Fdata.inp"):
            Fdata_inp_path = os.path.join(os.path.dirname(Fdata_path), "Fdata.inp")
            if not os.path.isfile(Fdata_inp_path):
                # check the existence of Fdata.inp
                raise FileNotFoundError
            else:
                os.symlink(Fdata_inp_path, "Fdata.inp", target_is_directory=False)

        Calculator.__init__(self, atoms=atoms, **kwargs)

    @property
    def atoms(self):
        return self._atoms

    @atoms.setter
    def atoms(self, atoms):
        if atoms is None:
            self._atoms = None
        else:
            self._atoms = atoms.copy()

    def calculate(self,
                  atoms=None,
                  properties=('energy',),
                  system_changes=tuple(all_changes)):
        """Do a fireball calculation in the specified directory.

        This will generate the necessary fireball input files, and then
        execute fireball. After execution, the energy, forces. etc. are read
        from the fireball output files.
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

        # self.update_atoms(atoms)
        self.read_results()

    def _run(self, command=None, directory=None):
        """Method to explicitly execute Fireball"""
        if command is None:
            command = self.command
        if directory is None:
            directory = self.directory
        errorcode = subprocess.call(command,
                                    shell=True,
                                    cwd=directory)
        return errorcode

    def read_results(self):
        output = self.sname + ".log.json"
        result = jsonio.read_json(output)
        self.results = result['fireball'][-1]
