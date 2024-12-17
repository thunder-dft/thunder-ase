import numpy as np
import spglib
from ase.dft.kpoints import kpoint_convert
from thunder_ase.utils.basis_set import fit_gaussian_command, fit_gaussian


def atoms2spg(atoms):
    """
    Convert ASE Atoms object to spglib input cell format
    :param atoms: ASE Atoms object
    :return: spglib input cell format
    """
    lattice = np.array(atoms.get_cell().T, dtype='double', order='C')
    positions = np.array(atoms.get_scaled_positions(), dtype='double', order='C')
    numbers = np.array(atoms.get_atomic_numbers(), dtype='intc')
    cell = (lattice, positions, numbers)
    return cell


def get_kpts(atoms, size=None, offset=None, reduced=True, **kwargs):
    """
    Get irreducible kpoints.
    Reference:
    * https://spglib.readthedocs.io/en/latest/python-spglib.html#methods-kpoints
    :param reduced:
    :param offset:
    :param atoms:
    :param size:
    :return: np.array with shape (N, 4), coordinates and weights
    """
    if 'symprec' in kwargs:
        symprec = kwargs['symprec']
    else:
        symprec = 1e-05  # for Cartesian coordinate in real space
    if 'gamma' in kwargs:
        gamma = kwargs['gamma']
    else:
        gamma = True
    if gamma:
        is_shift = [0, 0, 0]
        if offset is None:
            offset = [0., 0., 0.]
    else:
        is_shift = [1, 1, 1]
        if offset is None:
            offset = [0.5, 0.5, 0.5]

    mapping, grid = spglib.get_ir_reciprocal_mesh(size, atoms2spg(atoms), is_shift=is_shift, symprec=symprec)
    N = len(grid)
    if reduced:
        ir_idx = np.unique(mapping)
        ir_grid = (grid[ir_idx] + offset) / np.asarray(size)
    else:
        ir_idx = range(N)
        ir_grid = (grid + offset) / np.asarray(size)
    kpoints = kpoint_convert(cell_cv=atoms.cell, skpts_kc=ir_grid)  # convert from scaled and cartesian coordinates
    weights = [sum(mapping == i) / N for i in ir_idx]

    kpoints_weight = [[k[0], k[1], k[2], w] for k, w in zip(kpoints, weights)]
    return kpoints_weight
