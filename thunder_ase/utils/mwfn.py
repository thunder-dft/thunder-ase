from string import Template

MWFN_TEMPLATE = Template(
    """# Generate by Thunder-ASE
Wfntype=$wfntype
Charge=$charge
Naelec=$naelec
Nbelec=$nbelec
E_tot=$e_tot
VT_ratio=$vt_ratio
Ndim=$ndim
cellv1= $cellv1
cellv2= $cellv2
cellv3= $cellv3

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
$v
$$Primitive exponents
$primitive_exponents
$$Contraction coefficients
$contraction_coefficients
    
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
    'wfntype': '{:4i}',
    'charge': '{:15.6f}',
    'naelec': '{:15.6f}',
    'nbelec': '{:15.6f}',
    'e_tot': '{:16.8E}',
    'vt_ratio': '{:12.8f}',
    'ndim': '{:4i}',
    'cellv1': '{:12.8f}{:12.8f}{:12.8f}',
    'cellv2': '{:12.8f}{:12.8f}{:12.8f}',
    'cellv3': '{:12.8f}{:12.8f}{:12.8f}',
    'ncenter': '{:8i}',
    'atoms_coord': '{:6i} {:<2s}{:4i}{:6.1f}{:16.8f}{:16.8f}{:16.8f}',
    'nbasis': '{:8i}',
    'nindbasis': '{:8i}',
    'nprims': '{:8i}',
    'nshell': '{:8i}',
    'nprimshell': '{:8i}',
    'shell_types': '{:3i}',
    'shell_centers': '{:8i}',
    'shell_contraction_degress': '{:4i}',
    'primitive_exponents': '{:16.8E}',
    'contraction_coefficients': '{:16.8E}',
}

KEY_LINE_LEN = {
    'shell_types': 25,
    'shell_centers': 10,
    'shell_contraction_degress': 20,
    'primitive_exponents': 5,
    'contraction_coefficients': 5,
}


def format_data(key, data):
    maxlen = KEY_LINE_LEN[key]
    res = len(data) % maxlen
    nline = int(len(data) / maxlen)
    result = ''
    for i in range(nline):
        result += (MWFN_FORMAT[key] * maxlen).format(data[i*maxlen:(i+1)*maxlen]) + '\n'
    result += (MWFN_FORMAT[key] * res).format(data[-res:])
    return result
