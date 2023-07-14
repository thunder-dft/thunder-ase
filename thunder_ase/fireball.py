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
from ase.units import Hartree
import numpy as np
from ase.calculators.calculator import Calculator, CalculationFailed, all_changes
from ase.dft.kpoints import monkhorst_pack, kpoint_convert
from ase.geometry import get_distances
import ase.spacegroup
from ase.io import jsonio
from ase.dft.kpoints import BandPath
from thunder_ase.utils.mwfn import MWFN_FORMAT, MWFN_DEFAULT, MWFN_TEMPLATE, \
    CELL_TEMPLATE, format_data, read_cdcoeffs, reorder_cdcoeffs
from thunder_ase.utils.basis_set import read_info, read_gaussian
from thunder_ase.utils.shell_dict import SHELL_PRIMARY_NUMS, SHELL_PRIMITIVE


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
        return [[kx, ky, kz, 1.0 / nkpt] for kx, ky, kz in kpoints]

    if 'symprec' in kwargs:
        symprec = kwargs['symprec']
    else:
        symprec = 1e-05  # for Cartesian coordinate in real space
    try:
        # get spacegroup from atoms ase.spacegroup.get_spacegroup
        sg = ase.spacegroup.get_spacegroup(atoms, symprec=symprec)
    except RuntimeError:
        print("Didn't find symmetry. No reducing is needed.")
        return [[kx, ky, kz, 1.0 / nkpt] for kx, ky, kz in kpoints]

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
    'qstate': {'type': (int, float), 'name': 'qstate', 'default': 0},  # Extra electron. +1 mean adding 1 electron.
    'iquench': {'type': (int,), 'name': 'iquench', 'default': 0},
    # optimization algorithm, 0: MD, -1: hard quench at KE achieve maximum, -3: power (force * velocity) quench,
    # positive number N: quench every N step
    't_initial': {'type': (int, float), 'name': 'T_initial', 'default': 300.0},
    't_final': {'type': (int, float), 'name': 'T_final', 'default': 0.0},
    't_want': {'type': (int, float), 'name': 'T_want', 'default': 300.0},
    'taurelax': {'type': (int, float), 'name': 'taurelax', 'default': 5.0},
    'efermi_t': {'type': (int, float), 'name': 'efermi_T', 'default': 100.0},
    'dt': {'type': (int, float), 'name': 'dt', 'default': 0.25},  # fs
    'iensemble': {'type': (int,), 'name': 'iensemble', 'default': 0},
    'iconstraint_rcm': {'type': (int,), 'name': 'iconstraint_rcm', 'default': 1},  # shift molecule to COM or not
    'iconstraint_vcm': {'type': (int,), 'name': 'iconstraint_vcm', 'default': 1},  # whether keep COM fixed during md
    'iconstraint_l': {'type': (int,), 'name': 'iconstraint_L', 'default': 0},
    'iconstraint_ke': {'type': (int,), 'name': 'iconstraint_KE', 'default': 1},
    'ifix_neighbors': {'type': (int,), 'name': 'ifix_neighbors', 'default': 0},
    'ifix_charges': {'type': (int,), 'name': 'ifix_CHARGES', 'default': 1},
    'max_scf_iterations_set': {'type': (int,), 'name': 'max_scf_iterations_set', 'default': 50},
    'scf_tolerance_set': {'type': (int, float), 'name': 'scf_tolerance_set', 'default': 0.000001},
    'beta_set': {'type': (int, float), 'name': 'beta_set', 'default': 0.08},  # mix factor
    'ecut_set': {'type': (int, float), 'name': 'Ecut_set', 'default': 200.0},  # control mesh grid number
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
    'kpt_size': {'type': (list, np.array), 'name': 'kpt_size', 'default': None},
    'kpt_offset': {'type': (list, np.array), 'name': 'kpt_offset', 'default': None},
    'kpt_interval': {'type': (list, np.array, float, int), 'name': 'kpt_interval', 'default': None},
    'kpt_path': {'type': (BandPath, str, list, np.array), 'name': 'kpt_path', 'default': None},
    'nkpt': {'type': (int,), 'name': 'nkpt', 'default': None},  # n kpoints on path. Used if kpt_path is string.
    'kpt_reduced': {'type': (bool,), 'name': 'kpt_reduced', 'default': False},
}

fireball_params = options_params | output_params | xsfoptions_params | calc_params


def get_params_from_string(s):
    k, v = s.split('=')
    k = k.strip().lower()
    v = v.strip()
    if k not in fireball_params:
        raise KeyError
    vtype = fireball_params[k]['type']
    if float in vtype:
        v = float(v)
    else:
        v = vtype[0](v)
    return k, v


def write_params(dct, f):
    for k, v in dct.items():
        kname = fireball_params[k]['name']
        f.write("{} = {}\n".format(kname, v))


class GenerateFireballInput:
    def __init__(self, atoms=None, **kwargs):
        if atoms is not None:
            self._atoms = atoms.copy()
        else:
            self._atoms = None
        self.sname = '001'
        self.output_params = {}
        self.options_params = {}
        self.xsfoptions_params = {}

        self.kpt_size = fireball_params['kpt_size']['default']
        self.kpt_interval = fireball_params['kpt_interval']['default']
        self.kpt_offset = fireball_params['kpt_offset']['default']
        self.kpt_path = fireball_params['kpt_path']['default']
        self.nkpt = fireball_params['nkpt']['default']
        self.kpt_reduced = fireball_params['kpt_reduced']['default']
        self._kpoints = None

        self.check_kwargs(kwargs)

    @property
    def atoms(self):
        return self._atoms

    @atoms.setter
    def atoms(self, atoms):
        self._atoms = atoms.copy()

    def get_kpoints(self, reduced=True, **kwargs):
        if self._kpoints is not None:
            return self._kpoints

        if not np.all(self.atoms.pbc):
            return None

        if self.kpt_path is not None:
            if type(self.kpt_path) in (str,):
                lat = self.atoms.cell.get_bravais_lattice()
                path = lat.bandpath(self.kpt_path, npoints=self.nkpt)
                self.kpt_path = path
                kpts = path.kpts
            elif type(self.kpt_path) in (np.array, list, tuple):
                kpts = self.kpt_path
            elif type(self.kpt_path) in (BandPath,):  # TODO: kpt_path must be BandPath
                kpts = self.kpt_path.kpts
            else:
                raise TypeError

            kpts = kpoint_convert(cell_cv=self.atoms.cell, skpts_kc=kpts)  # convert to Cartesian coordinate
            self.nkpt = len(kpts)
            self._kpoints = [[k[0], k[1], k[2], 1.0 / self.nkpt] for k in kpts]
            return self._kpoints

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

        self._kpoints = get_kpts(self.atoms, size=size, offset=offset, reduced=reduced, **kwargs)
        self.nkpt = len(self._kpoints)
        return self._kpoints

    def get_k_point_weights(self):
        return np.asarray(self.get_kpoints())[:, -1]

    def set_kpoints(self, kpoints):
        self._kpoints = kpoints

    def write_Fdata_inp(self, atoms=None, Fdata_path=None):
        if atoms is None:
            atoms = self.atoms
        if Fdata_path is None:
            Fdata_path = 'Fdata'
        species = sorted([i for i in set(atoms.numbers)])

        # TODO: if Fdata_path longer than 128, make a symbolic link.
        with open('Fdata.inp', 'w') as f:
            f.write("{}\n".format(len(species)))
            for s in species:
                f.write("{}\n".format(s))
            f.write("'{}'".format(Fdata_path))

    def check_kwargs(self, kwargs):
        for key, v in kwargs.items():
            k = key.lower()
            if k not in fireball_params:
                print("The option {} not supported!".format(k))
                raise KeyError
            if type(v) not in fireball_params[k]['type']:
                print("The type of {} should be {}, but type {} is provided!".format(k, ' or '.join(
                    [i.__name__ for i in fireball_params[k]['type']]), type(v).__name__))
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
            elif k == 'kpt_path':
                self.kpt_path = v
            elif k == 'nkpt':
                self.nkpt = v
            elif k == 'kpt_reduced':
                self.kpt_reduced = v
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

    def write_atoms(self, atoms=None, pbc=None):
        """
        TODO: for non-PBC system, set cell to [0.000] * 9
        :param atoms:
        :param pbc:
        :return:
        """
        if atoms is not None:
            self.atoms = atoms.copy()

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
        with open(self.sname + '.KPOINTS', 'w') as f:
            f.write("{}\n".format(len(self.get_kpoints(reduced=reduced, **kwargs))))
            for k in self.get_kpoints(reduced=reduced, **kwargs):
                f.write("{:8.6f} {:8.6f} {:8.6f} {:8.6f}\n".format(*k))

    def write_input(self, atoms=None, Fdata_path=None):
        if atoms is not None:
            self.atoms = atoms.copy()
        self.write_Fdata_inp(atoms=self.atoms, Fdata_path=Fdata_path)
        self.write_options()
        self.write_atoms(pbc=self.atoms.pbc)
        if np.any(self.atoms.pbc):
            self.write_kpts(reduced=self.kpt_reduced)

    def read_options(self, input_file='structures.inp', read_atoms=False):
        # read structures.inp, get names for atoms and kpoints
        with open(input_file, 'r') as f:
            lines = f.readlines()

        natoms_list = int(lines[0].strip())
        sname_list = []
        for iline in lines[1:natoms_list + 1]:
            sname_list.append(iline.strip().strip('.inp'))

        loption, loutput, lxsfoptions = False, False, False
        for line in lines[natoms_list + 1:]:
            content = line.strip()
            if '!' in content:
                content = content.split('!')[0]
            if len(content) == 0:
                continue

            if '&OPTIONS' in content:
                loption = True
                continue
            if '&OUTPUT' in content:
                loutput = True
                continue
            if "&XSFOPTIONS" in content:
                lxsfoptions = True
                continue

            if loption:
                if '&END' in content:
                    loption = False
                    continue
                k, v = get_params_from_string(content)
                self.options_params[k] = v

            elif loutput:
                if '&END' in content:
                    loutput = False
                    continue
                k, v = get_params_from_string(content)
                self.output_params[k] = v
            elif lxsfoptions:
                if '&END' in content:
                    lxsfoptions = False
                    continue
                k, v = get_params_from_string(content)
                self.xsfoptions_params[k] = v

        if read_atoms:
            return sname_list

    def read_kpts(self, input_file=None):
        with open(input_file, 'r') as f:
            lines = f.readlines()
        nkpts = int(lines[0].strip())
        kpts = [list(map(float, line.strip().split())) for line in lines[1:] if len(line.strip()) > 0]
        if nkpts != len(kpts):
            raise ValueError("The kpoints number is not inconsistent.")
        self.set_kpoints(kpts)


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
        self._eigenvalues = None
        self._shell_info = None
        GenerateFireballInput.__init__(self, **kwargs)  # TODO: set kwargs use set() function, see vasp calculator
        Calculator.__init__(self, atoms=atoms, **kwargs)
        if not os.path.isdir(Fdata_path):
            # check the existence of Fdata directory
            print("Error: Can't find Fdata! Pls Check the setting!")
            raise FileNotFoundError
        self.Fdata_path = Fdata_path
        self.command = command

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

        self.write_input(Fdata_path=self.Fdata_path)

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
        output = self.sname + ".json"
        result = jsonio.read_json(output)
        if 'charges' in result['fireball'][-1]:
            shell_charges = result['fireball'][-1]['charges']
            result['fireball'][-1]['shell_charges'] = shell_charges
            del(result['fireball'][-1]['charges'])

        self.results = result['fireball'][-1]

    def get_charges(self, atoms=None):
        if self.results is None:
            try:
                self.read_results()
            except AttributeError:
                return None
        else:
            if 'charges' in self.results:
                return self.results['charges']

        shell_charges = self.results['shell_charges']
        sum_charges = [sum(isc) for isc in shell_charges]
        if atoms is None:
            atoms = self.atoms
        ref_charges = [sum(self.shell_info[s]['occupation']) for s in atoms.symbols]
        charges = [ref - sc for sc, ref in zip(sum_charges, ref_charges)]
        self.results['charges'] = charges
        return charges

    def get_forces(self, atoms=None):
        forces = self.get_property('forces', atoms)
        return np.asarray(forces)

    def get_fermi_level(self):
        return self.get_property('fermi')

    def get_ibz_k_points(self):
        if self.kpt_path is None:
            raise ValueError
        kpts = np.asarray(self.get_kpoints())[:, :3]
        return kpts

    def get_eigenvalues(self, kpt=0, spin=0):
        if self._eigenvalues is not None:
            return self._eigenvalues[spin][kpt]

        if not self.results:
            self.calculate()
        # Read eigen file
        filename = self.sname + '.eigen'
        if not os.path.isfile(filename):
            raise FileExistsError

        with open(filename, 'r') as f:
            lines = f.readlines()

        eigenvalues = []
        eigen = []
        kpts_count = [str(k + 1) for k in range(self.nkpt)]
        for line in lines:
            content = [i.strip() for i in line.strip().split() if len(i.strip()) != 0]
            if len(content) == 0:
                pass
            if len(kpts_count) != 0:
                if content[0] == kpts_count[0]:
                    if len(eigen) != 0:
                        eigenvalues.append(eigen)
                    eigen = []
                    kpts_count.pop(0)
                else:
                    eigen += [float(i) for i in content]
            else:
                eigen += [float(i) for i in content]
        eigenvalues.append(eigen)

        self._eigenvalues = [eigenvalues]  # only one spin

        return self._eigenvalues[spin][kpt]

    def get_number_of_spins(self):
        """
        Current Fireball doesn't support spin polarization.
        :return:
        """
        return 1

    def band_structure(self, atoms=None, reference=None):
        """Create band-structure object for plotting."""
        if atoms is not None:
            self.atoms = atoms.copy()

        from ase.spectrum.band_structure import get_band_structure
        if self.kpt_path is None:
            raise ValueError

        if type(self.kpt_path) in (BandPath,):
            path = self.kpt_path
        else:
            path = None
        return get_band_structure(atoms=self.atoms, calc=self, path=path, reference=reference)

    def read_options(self, input_file='structures.inp', read_atoms=False):
        sname_list = super().read_options(input_file=input_file, read_atoms=read_atoms)
        if not read_atoms or len(sname_list) == 0:
            return None

        if len(sname_list) == 1:
            self.atoms = read_inp(input_file='{}.inp'.format(sname_list[0]))
            kpts_file = '{}.KPOINTS'.format(sname_list[0])
            if os.path.isfile(kpts_file):
                self.read_kpts(kpts_file)
            return None
            # multiple atoms: if it contains multiple atoms, return MultiFireball instead of Fireball
        print("Warning: This file contains multiple atoms, return MultiFireball instead of Fireball")
        atoms_list = []
        for sname in sname_list:
            # read atoms
            atoms = read_inp(input_file='{}.inp'.format(sname))
            # read kpts
            kpts_file = '{}.KPOINTS'.format(sname)
            if os.path.isfile(kpts_file):
                self.read_kpts(kpts_file)
            atoms.set_calculator(self)
            atoms_list.append(atoms)
        multicalc = MultiFireball(atoms_list=atoms_list, sname_list=sname_list, calc=self)
        return multicalc

    @property
    def shell_info(self):
        if self._shell_info is None:
            # read Fdata/info.dat
            self._shell_info = read_info(os.path.join(self.Fdata_path, 'info.dat'))

            # TODO: this is a dirty way to get exicted_label.
            #  The excited label can be obtained by set(ishell). However, still dirty.
            for symbol, value in self._shell_info.items():
                ishell = value['shells']
                len_ground = len(SHELL_PRIMARY_NUMS[value['number']])
                assert len(set(ishell)) <= len_ground  # if not, something wrong during generating Fdata
                excited_label = [0] * len_ground + [1] * (len(ishell) - len_ground)
                self._shell_info[symbol]['excited'] = excited_label
        return self._shell_info

    def get_valence_charge(self, i):
        isymbol = self.atoms.symbols[i]
        shell_info = self.shell_info[isymbol]
        return sum(shell_info['occupation'])

    def write_mwfn(self, kpoint=0, filename=None):
        mwfn_dict = MWFN_DEFAULT.copy()
        # Initialize default data format
        for k, v in mwfn_dict.items():
            if v is not None:
                mwfn_dict[k] = MWFN_FORMAT[k].format(v)

        # atom information
        mwfn_dict['ncenter'] = MWFN_FORMAT['ncenter'].format(len(self.atoms))
        atoms_coord = [MWFN_FORMAT['atoms_coord'].format(idx + 1,
                                                         iatom.symbol,
                                                         iatom.number,
                                                         self.get_valence_charge(idx),
                                                         *iatom.position)
                       for idx, iatom in enumerate(self.atoms)]
        mwfn_dict['atoms_coord'] = '\n'.join(atoms_coord)

        # cell info
        if np.any(self.atoms.pbc):
            cell_info = {
                'ndim': MWFN_FORMAT['ndim'].format(3),
                'cellv1': (MWFN_FORMAT['cellv1'] * 3).format(*self.atoms.cell[0]),
                'cellv2': (MWFN_FORMAT['cellv2'] * 3).format(*self.atoms.cell[1]),
                'cellv3': (MWFN_FORMAT['cellv3'] * 3).format(*self.atoms.cell[2]),
            }
            mwfn_dict['cell_info'] = CELL_TEMPLATE.substitute(cell_info)
        else:
            mwfn_dict['cell_info'] = ''

        # electron for alpha and beta
        tot_elec = sum([self.get_valence_charge(i) for i in range(len(self.atoms))]) - sum(self.get_charges())
        naelec = 0.5 * tot_elec
        nbelec = 0.5 * tot_elec
        mwfn_dict['naelec'] = MWFN_FORMAT['naelec'].format(naelec)
        mwfn_dict['nbelec'] = MWFN_FORMAT['nbelec'].format(nbelec)
        # total energy
        mwfn_dict['e_tot'] = MWFN_FORMAT['e_tot'].format(self.get_potential_energy() / Hartree)
        # Basis set info
        shell_types = []
        shell_centers = []
        shell_contraction_degress = []
        primitive_exponents = []
        contraction_coefficients = []
        for idx, symbol in enumerate(self.atoms.symbols):
            ishells = self.shell_info[symbol]['shells']
            shell_centers += [idx + 1] * len(ishells)
            shell_types += ishells
            excited_label = self.shell_info[symbol]['excited']
            for shell, is_excited in zip(ishells, excited_label):
                gbs = read_gaussian(self.shell_info[symbol]['number'], shell, is_excited, Fdata_path=self.Fdata_path)
                shell_contraction_degress.append(gbs['degree'])
                primitive_exponents += gbs['alpha'].tolist()
                contraction_coefficients += gbs['coefficient'].tolist()
        nbasis = sum([shell * 2 + 1 for shell in shell_types])
        nshell = len(shell_types)
        nprimshell = len(primitive_exponents)
        nprims = sum([SHELL_PRIMITIVE[st] * degree
                      for st, degree in zip(shell_types, shell_contraction_degress)])
        shell_types = [((i < 2) * 2 - 1) * i for i in shell_types]  # if i >= 2, use negative value

        mwfn_dict['nbasis'] = MWFN_FORMAT['nbasis'].format(nbasis)
        mwfn_dict['nindbasis'] = MWFN_FORMAT['nindbasis'].format(nbasis)
        mwfn_dict['nprims'] = MWFN_FORMAT['nprims'].format(nprims)
        mwfn_dict['nshell'] = MWFN_FORMAT['nshell'].format(nshell)
        mwfn_dict['nprimshell'] = MWFN_FORMAT['nprimshell'].format(nprimshell)

        mwfn_dict['shell_types'] = format_data('shell_types', shell_types)
        mwfn_dict['shell_centers'] = format_data('shell_centers', shell_centers)
        mwfn_dict['shell_contraction_degress'] = format_data('shell_contraction_degress', shell_contraction_degress)
        mwfn_dict['primitive_exponents'] = format_data('primitive_exponents', primitive_exponents)
        mwfn_dict['contraction_coefficients'] = format_data('contraction_coefficients', contraction_coefficients)

        cdcoeffs = read_cdcoeffs(self.sname + '.cdcoeffs-mwfn')[kpoint]
        orbital_coeffs = []
        for orbital in cdcoeffs:
            info = ''.join(orbital['info'])
            coeff = orbital['coeff']
            coeff = reorder_cdcoeffs(coeff, shell_types)
            orbital_coeffs += [info + format_data('MO_coefficients', coeff)]
        mwfn_dict['orbital_coeffs'] = '\n\n'.join(orbital_coeffs)

        self.mwfn_dict = mwfn_dict
        content = MWFN_TEMPLATE.substitute(mwfn_dict)
        # TODO: read orbital info, append to the content

        # write out
        if filename is None:
            filename = self.sname + '.mwfn'
        with open(filename, 'w') as f:
            f.write(content)


class MultiFireball:
    name = 'multi_fireball'

    def __init__(self, atoms_list=None, calc=None, sname_list=None):
        if calc is not None:
            self.calc = calc
        elif atoms_list is not None:
            if atoms_list[0].calc is not None:
                self.calc = atoms_list[0].calc
        else:
            self.calc = None

        _ = [atoms.set_calc(self.calc) for atoms in atoms_list if atoms.calc is None]

        self.atoms_list = atoms_list

        if sname_list is None:
            self.sname_list = ["{:03d}".format(idx + 1) for idx in range(len(self.atoms_list))]

    def write_Fdata_inp(self, Fdata_path=None):
        if Fdata_path is None:
            Fdata_path = 'Fdata'
        species = sorted(list(set([i for atoms in self.atoms_list for i in set(atoms.numbers)])))
        with open('Fdata.inp', 'w') as f:
            f.write("{}\n".format(len(species)))
            for s in species:
                f.write("{}\n".format(s))
            f.write("'{}'".format(Fdata_path))

    def write_options(self):
        with open('structures.inp', 'w') as f:
            f.write("{:d}\n".format(len(self.atoms_list)))
            for sname in self.sname_list:
                f.write("{}.inp\n".format(sname))

            f.write("! Write out options - leave this line as comment line or null\n")

            f.write("&OUTPUT\n")
            write_params(self.calc.output_params, f)
            f.write("&END\n")

            f.write("&OPTIONS\n")
            write_params(self.calc.options_params, f)
            f.write("&END\n")

            f.write("&XSFOPTIONS\n")
            write_params(self.calc.xsfoptions_params, f)
            f.write("&END\n")

    def write_input(self):
        Fdata_path = self.calc.Fdata_path
        self.write_Fdata_inp(Fdata_path=Fdata_path)
        self.write_options()
        for atoms, sname in zip(self.atoms_list, self.sname_list):
            atoms.calc.sname = sname
            atoms.calc.write_atoms(atoms=atoms, pbc=atoms.pbc)
            atoms.calc.write_kpts()

    def calculate(self):
        errorcode = self._run(command=self.calc.command,
                              directory=self.calc.directory)
        if errorcode:
            raise CalculationFailed(
                '{} in {} returned an error: {:d}'.format(
                    self.name, self.calc.directory, errorcode))
        for atoms in self.atoms_list:
            atoms.calc.read_results()

    def _run(self, command=None, directory=None):
        """Method to explicitly execute Fireball"""
        if command is None:
            command = self.calc.command
        if directory is None:
            directory = self.calc.directory
        errorcode = subprocess.call(command, shell=True, cwd=directory)
        return errorcode


def read_inp(input_file):
    # read fireball input file, return ase.Atoms
    with open(input_file, 'r') as f:
        lines = f.readlines()

    natoms, is_cluster = lines[0].split()
    pbc = not bool(int(is_cluster))
    cell = [[float(i) for i in lines[iline + 1].split()] for iline in range(3)]
    numbers = []
    positions = []
    for line in lines[4:]:
        content = line.strip().split()
        if len(content) == 4:
            num, x, y, z = content
            numbers.append(int(num))
            positions.append(list(map(float, [x, y, z])))

    atoms = ase.Atoms(numbers=numbers,
                      positions=positions,
                      cell=cell,
                      pbc=pbc)
    return atoms
