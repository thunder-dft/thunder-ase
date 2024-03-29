SHELL_NAME = {
    0: 's',
    1: 'p',
    2: 'd',
    3: 'f',
}

SHELL_NUM = {
    's': 0,
    'p': 1,
    'd': 2,
    'f': 3,
}

SHELL_PRIMARY_NUMS = {
    1: [1, ],  # H
    2: [1, ],  # He
    3: [2, 2, ],  # Li
    4: [2, 2, ],  # Be
    5: [2, 2, ],  # B
    6: [2, 2, ],  # C
    7: [2, 2, ],  # N
    8: [2, 2, ],  # O
    9: [2, 2, ],  # F
    10: [2, 2, ],  # Ne
    11: [3, 3, 3, ],  # Na
    12: [3, 3, 3, ],  # Mg
    13: [3, 3, 3, ],  # Al
    14: [3, 3, 3, ],  # Si
    15: [3, 3, 3, ],  # P
    16: [3, 3, 3, ],  # S
    17: [3, 3, 3, ],  # Cl
    18: [3, 3, 3, ],  # Ar
    19: [4, 4, 3, ],  # K
    20: [4, 4, 3, ],  # Ca
    21: [4, 4, 3, ],  # Sc
    22: [4, 4, 3, ],  # Ti
    23: [4, 4, 3, ],  # V
    24: [4, 4, 3, ],  # Cr
    25: [4, 4, 3, ],  # Mn
    26: [4, 4, 3, ],  # Fe
    27: [4, 4, 3, ],  # Co
    28: [4, 4, 3, ],  # Ni
    29: [4, 4, 3, ],  # Co
    30: [4, 4, 3, ],  # Zn
    31: [4, 4, 3, ],  # Ga
    32: [4, 4, 3, ],  # Ge
    33: [4, 4, 3, ],  # As
    34: [4, 4, 3, ],  # Se
    35: [4, 4, 3, ],  # Br
    36: [4, 4, 3, ],  # Kr
    37: [5, 5, 4, ],  # Rb
    38: [5, 5, 4, ],  # Sr
    39: [5, 5, 4, ],  # Y
    40: [5, 5, 4, ],  # Zr
    41: [5, 5, 4, ],  # Nb
    42: [5, 5, 4, ],  # Mo
    43: [5, 5, 4, ],  # Tc
    44: [5, 5, 4, ],  # Ru
    45: [5, 5, 4, ],  # Rh
    46: [5, 5, 4, ],  # Pd
    47: [5, 5, 4, ],  # Ag
    48: [5, 5, 4, ],  # Cd
    49: [5, 5, 4, ],  # In
    50: [5, 5, 4, ],  # Sn
    51: [5, 5, 4, ],  # Sb
    52: [5, 5, 4, ],  # Te
    53: [5, 5, 4, ],  # I
    54: [5, 5, 4, ],  # Xe
    55: [6, 6, 5, 4, ],  # Cs
    56: [6, 6, 5, 4, ],  # Ba
    57: [6, 6, 5, 4, ],  # La
    58: [6, 6, 5, 4, ],  # Ce
    59: [6, 6, 5, 4, ],  # Pr
    60: [6, 6, 5, 4, ],  # Nd
    61: [6, 6, 5, 4, ],  # Pm
    62: [6, 6, 5, 4, ],  # Sm
    63: [6, 6, 5, 4, ],  # Eu
    64: [6, 6, 5, 4, ],  # Gd
    65: [6, 6, 5, 4, ],  # Tb
    66: [6, 6, 5, 4, ],  # Dy
    67: [6, 6, 5, 4, ],  # Ho
    68: [6, 6, 5, 4, ],  # Er
    69: [6, 6, 5, 4, ],  # Tm
    70: [6, 6, 5, 4, ],  # Yb
    71: [6, 6, 5, 4, ],  # Lu
    72: [6, 6, 5, 4, ],  # Hf
    73: [6, 6, 5, 4, ],  # Ta
    74: [6, 6, 5, 4, ],  # W
    75: [6, 6, 5, 4, ],  # Re
    76: [6, 6, 5, 4, ],  # Os
    77: [6, 6, 5, 4, ],  # Ir
    78: [6, 6, 5, 4, ],  # Pt
    79: [6, 6, 5, 4, ],  # Au
    80: [6, 6, 5, 4, ],  # Hg
    81: [6, 6, 5, 4, ],  # Tl
    82: [6, 6, 5, 4, ],  # Pb
    83: [6, 6, 5, 4, ],  # Bi
    84: [6, 6, 5, 4, ],  # Po
    85: [6, 6, 5, 4, ],  # At
    86: [6, 6, 5, 4, ],  # Rn
}

EXCITED_SHELL_PRIMARY_NUMS = {k: [i + 1 for i in v] for k, v in SHELL_PRIMARY_NUMS.items()}


SHELL_PRIMITIVE = {
    0: 1,  # s
    1: 3,  # p
    2: 6,  # d
    3: 10,  # f
}
