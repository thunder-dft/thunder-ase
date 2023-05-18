def ordinal(n: int):
    """
    Convert an integer to "1st", "2nd", etc,
    :param n: int
    :return: str
    """
    if 11 <= (n % 100) <= 13:
        suffix = 'th'
    else:
        suffix = ['th', 'st', 'nd', 'rd', 'th'][min(n % 10, 4)]
    return str(n) + suffix


from thunder_ase.utils.basis_set import fit_gaussian
