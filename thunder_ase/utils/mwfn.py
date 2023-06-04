from string import Template

import numpy as np

CELL_TEMPLATE = Template(
"""Ndim=$ndim
cellv1= $cellv1
cellv2= $cellv2
cellv3= $cellv3
"""
)

MWFN_TEMPLATE = Template(
    """# Generate by Thunder-ASE
Wfntype=$wfntype
Charge=$charge
Naelec=$naelec
Nbelec=$nbelec
E_tot=$e_tot
VT_ratio=$vt_ratio
$cell_info
# Atom information
Ncenter=$ncenter
$$Centers
$atoms_coord

# Basis function information
Nbasis=    $nbasis
Nindbasis= $nindbasis
Nprims=    $nprims
Nshell=    $nshell
Nprimshell=$nprimshell
$$Shell types
$shell_types
$$Shell centers
$shell_centers
$$Shell contraction degrees
$shell_contraction_degress
$$Primitive exponents
$primitive_exponents
$$Contraction coefficients
$contraction_coefficients

# Orbital information (nindbasis orbitals)

$orbital_coeffs
"""
)

MWFN_DEFAULT = {
    'wfntype': 0,
    'charge': 0.0,
    'naelec': None,
    'nbelec': None,
    'e_tot': 0.0,
    'vt_ratio': 0.0,
    'ndim': None,
    'cellv1': None,
    'cellv2': None,
    'cellv3': None,
    'ncenter': None,
    'atoms_coord': None,
    'nbasis': None,
    'nindbasis': None,
    'nprims': None,
    'nshell': None,
    'nprimshell': None,
    'shell_types': None,
    'shell_centers': None,
    'shell_contraction_degress': None,
    'primitive_exponents': None,
    'contraction_coefficients': None,
}

MWFN_FORMAT = {
    'wfntype': '{:4d}',
    'charge': '{:15.6f}',
    'naelec': '{:15.6f}',
    'nbelec': '{:15.6f}',
    'e_tot': '{:16.8E}',
    'vt_ratio': '{:12.8f}',
    'ndim': '{:4d}',
    'cellv1': '{:12.8f}',
    'cellv2': '{:12.8f}',
    'cellv3': '{:12.8f}',
    'ncenter': '{:8d}',
    'atoms_coord': '{:6d} {:<2s}{:4d}{:6.1f}{:16.8f}{:16.8f}{:16.8f}',
    'nbasis': '{:8d}',
    'nindbasis': '{:8d}',
    'nprims': '{:8d}',
    'nshell': '{:8d}',
    'nprimshell': '{:8d}',
    'shell_types': '{:3d}',
    'shell_centers': '{:8d}',
    'shell_contraction_degress': '{:4d}',
    'primitive_exponents': '{:16.8E}',
    'contraction_coefficients': '{:16.8E}',
    'MO_coefficients': '{:16.8E}',
}

KEY_LINE_LEN = {
    'shell_types': 25,
    'shell_centers': 10,
    'shell_contraction_degress': 20,
    'primitive_exponents': 5,
    'contraction_coefficients': 5,
    'MO_coefficients': 5,
}


def format_data(key, data):
    maxlen = KEY_LINE_LEN[key]
    res = len(data) % maxlen
    nline = int(len(data) / maxlen)
    result = [(MWFN_FORMAT[key] * maxlen).format(*data[i*maxlen:(i+1)*maxlen])
              for i in range(nline)]
    result.append((MWFN_FORMAT[key] * res).format(*data[-res:]))
    result = [i for i in result if len(i.strip()) > 0]
    return '\n'.join(result)


def read_cdcoeffs(filename):
    with open(filename, 'r') as f:
        lines = f.readlines()

    block_start, block_end = [], []
    for il, line in enumerate(lines):
        if "Kpoint=" in line:
            block_start.append(il+1)
            if il > 0:
                block_end.append(il-1)
    block_end.append(len(lines))

    coeffs_kpoints = []
    for start, end in zip(block_start, block_end):
        coeff_orbitals = []
        coeff_orbital = {
            'info': [],
            'coeff': []
        }
        for line in lines[start:end]:
            content = line.strip()
            if len(content) == 0:
                if len(coeff_orbital['info']) > 0:
                    coeff_orbitals.append(coeff_orbital)
                coeff_orbital = {
                    'info': [],
                    'coeff': []
                }
                continue
            if 'Index' in content:
                start_orbital = True
            elif '$Coeff' in content:
                coeff_orbital['info'].append(line)
                start_orbital = False
                continue
            if start_orbital:
                coeff_orbital['info'].append(line)
            else:
                coeff_orbital['coeff'] += [float(i) for i in content.split()]

        coeffs_kpoints.append(coeff_orbitals)
    return coeffs_kpoints


def reorder_cdcoeffs(coeff_list, shell_types):
    """
    Re-order the orbital from (Y, Z, X) and (D-2, D-1, D0, D+1, D+2) to
    (X, Y, Z) and (D0, D+1, D-1, D+2, D-2)

    :param coeff_list:
    :param shell_types:
    :return:
    """
    new_coeff_list = []
    idx_curr = 0
    for st in shell_types:
        if st == 0:  # s orbital
            new_coeff_list += [coeff_list[idx_curr]]
        elif st == 1:  # p orbital
            new_coeff_list += [
                coeff_list[idx_curr+2],
                coeff_list[idx_curr],
                coeff_list[idx_curr+1]
            ]
        elif st == -2:  # pure d orbital. The scale factor, see fireball/b.FUNCTIONS/clm.f90 and Ylm.f90
            new_coeff_list += [
                coeff_list[idx_curr+2],  # D0
                coeff_list[idx_curr+3],  # D+1
                coeff_list[idx_curr+1],  # D-1
                coeff_list[idx_curr+4],  # D+2
                coeff_list[idx_curr+0],  # D-2
            ]
        else:
            NotImplementedError

        idx_curr += np.abs(st) * 2 + 1

    return new_coeff_list
