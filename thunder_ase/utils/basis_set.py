import argparse
import os.path

import numpy as np
from scipy.optimize import minimize
import matplotlib.pyplot as plt
from thunder_ase.utils.shell_dict import SHELL_PRIMARY_NUMS, SHELL_NUM, SHELL_NAME
from thunder_ase.utils.ANO_DK3_GBS import ANO_DK3_GBS as GBS
from ase.data import chemical_symbols


# read wf file
def read_wf(filename):
    with open(filename, 'r') as f:
        lines = f.readlines()
    data = [list(map(float, line.strip().split())) for line in lines if len(line.strip().split()) == 2]
    return np.asarray(data)


def read_info(filename='Fdata/info.dat'):
    """
    :param filename:
    :return: dict, {'H': {
        'number': 1,
        'nshell': 1,
        'shells': [0, 1],
        'occupation': [1.0, 0.0],
    }}
    """
    with open(filename, 'r') as f:
        lines = f.readlines()

    result_dict = {}

    status = 0  # 0: start block, 1: in block, 2: need to end block
    for il, line in enumerate(lines):
        content = line.strip()
        if len(content) == 0:
            continue
        if content == '=' * 70:
            if status == 0:
                status = 1
            elif status == 2:
                status = 0
                continue

        if status == 1:
            result_dict[lines[il + 2].strip().split()[0]] = {
                'number': int(lines[il + 3].strip().split()[0]),
                'nshell': int(lines[il + 5].strip().split()[0]),
                'shells': [int(i) for i in lines[il + 6].strip().split()],
                'occupation': [float(i) for i in lines[il + 7].strip().split()],
            }
            status = 2
        else:
            continue

    return result_dict


# gaussian function
def gaussian(r, l=0, A=np.array([1.0]), a=np.array([0.1])):
    if l == 0:
        f = np.sum([Ai * np.exp(-ai * (r ** 2))
                    for Ai, ai in zip(A, a)], axis=0)
    else:
        f = np.sum([Ai * (r ** l) * np.exp(-ai * (r ** 2))
                    for Ai, ai in zip(A, a)], axis=0)
    return f


# loss function
def loss_function(x, *args):
    len_x = int(len(x) / 2)
    alpha = 10 ** x[0:len_x]
    A = np.asarray(x[len_x:])
    l, r, Y = args
    if r.shape != Y.shape:
        raise ValueError('Shapes for args[1] and args[2] are not equal!')
    # loss = np.abs(Y - gaussian(r, l, A, alpha)).sum()
    loss = np.mean((gaussian(r, l, A, alpha) - Y)**2)
    return loss


def loss_jac(x, *args):
    len_x = int(len(x) / 2)
    alpha = 10 ** x[0:len_x]
    A = np.asarray(x[len_x:])
    l, r, Y = args
    if r.shape != Y.shape:
        raise ValueError('Shapes for args[1] and args[2] are not equal!')
    if l == 0:
        d_alpha = np.asarray([Ai * (-1) * (r ** 2) * np.exp(-ai * (r ** 2)) * ai * np.log(10)
                              for Ai, ai in zip(A, alpha)])
        d_A = np.asarray([np.exp(-ai * (r ** 2)) for ai in alpha])
    else:
        d_alpha = np.asarray([Ai * (r ** l) * (-1) * (r ** 2) * np.exp(-ai * (r ** 2)) * ai * np.log(10)
                              for Ai, ai in zip(A, alpha)])
        d_A = np.asarray([(r ** l) * np.exp(-ai * (r ** 2)) for ai in alpha])
    der = np.concatenate([d_alpha, d_A])
    jac = 2 * np.mean((gaussian(r, l, A, alpha) - Y) * der, axis=1)
    return jac


def fit_wf(data, x0, l=0, bnds=None):
    R, Y = np.asarray(data)
    nz = int(len(x0) / 2)
    if bnds is not None:
        bnds = [bnds[0]] * nz + [bnds[1]] * nz
    res = minimize(loss_function, x0, bounds=bnds, args=(l, R, Y))
    alpha = res.x[0:nz]  # exponential parameters,
    Ae = res.x[nz:]  # coefficients
    Y_fit = gaussian(R, l, Ae, alpha)
    error = (np.linalg.norm(Y - Y_fit)) / len(Y)
    return [Ae, alpha, Y_fit, error]


def fit_wf_from_random(data, l=0, tol=None, Nzeta=None, bnds=None):
    """

    :param data: shape = [N, 2]
    :param l: principal quantum number
    :param tol: Error tolerance
    :param Nzeta: number of gaussian
    :param bnds: boundary for fitting
    :return:
    """
    if Nzeta is None:
        if tol is None:
            tol = 1.0E-5
        Nzeta0 = 3  # initial zeta number
        Nzeta = 20  # so the maximum Nzeta is 20
    else:
        Nzeta0 = Nzeta

    success = False
    Ae, ae, Y_fit, error = None, None, None, None
    for nz in range(Nzeta0, Nzeta + 1):
        A0 = np.random.rand(nz)  # initial value for A
        a0 = np.random.rand(nz)  # initial value for alpha0
        x0 = np.concatenate([A0, a0])  # initial para vector
        Ae, ae, Y_fit, error = fit_wf(data, x0, l=l, bnds=bnds)
        if tol is not None:
            if error < tol:
                success = True
                break

    if not success:
        if tol is None:
            print("Fitting error is {} for {} gaussians.".format(error, Nzeta))
        else:
            print("Warning: Fitting error {} didn't meet the tolerance {} for {} gaussians.".format(error, tol, Nzeta))

    return [Ae, ae, Y_fit]


def plot_fitting(data, l, Ae, ae, output=None):
    R, Y = np.asarray(data)
    plt.plot(R, Y, '*')
    plt.plot(R, gaussian(R, l, Ae, ae))
    Y_min = min(Y)
    Y_max = max(Y)
    Y_range = Y_max - Y_min
    ymin = Y_min - Y_range * 0.2
    ymax = Y_max + Y_range * 0.2
    plt.ylim([ymin, ymax])
    for Ai, ai in zip(Ae, ae):
        plt.plot(R, gaussian(R, l, [Ai], [ai]))
    if output is not None:
        plt.savefig(output)
    else:
        plt.show()


def get_initial_gbs_guess(element_num, primary_num, shell_name):
    # loop GBS, find minimum element contain excited_pn and shell_name
    pn = str(primary_num)
    for n_element in chemical_symbols[element_num:]:
        if n_element not in GBS:
            continue
        if pn not in GBS[n_element]:
            continue
        if shell_name in GBS[n_element][pn]:
            gauss = GBS[n_element][pn][shell_name]
            if n_element != chemical_symbols[element_num]:
                print("Initial guess for {} {}{} is from {}!".format(
                    chemical_symbols[element_num], pn, shell_name, n_element))
            return np.concatenate(np.asarray(gauss).T)

    print("No initial guess for {} {}{}!".format(chemical_symbols[element_num], pn, shell_name))
    print("Maybe too large element number!")
    raise NotImplementedError


def get_primary_number(number, shell, is_excited):
    # TODO: this is a dirty way to get primary number. It depends on same shell configuration.
    # The best way is write primary number in info.dat or wf.dat @JamesLewis .
    n = SHELL_PRIMARY_NUMS[number][shell]
    if is_excited == '1':
        n += 1
    return n


# run this after begin.x
def fit_gaussian(prog='fit_basis_set',
                 description='Fit Fireball basis set to Gaussian-type basis set.'):
    parser = argparse.ArgumentParser(
        prog=prog,
        description=description, )
    parser.add_argument('input_name', nargs='+', help='Fireball wave function file.')
    parser.add_argument('-p', '--plot', action='store_true')
    args = parser.parse_args()

    for input_name in args.input_name:
        name_list = input_name.split('.')  # format: '001.wf-s0.dat', element_number, shell, is_excited
        element_number = int(name_list[0])
        shell_name, is_excited = name_list[1][-2:]
        shell = SHELL_NUM[shell_name]
        if not os.path.exists(input_name):
            print("{} doesn't exist!".format(input_name))
            raise FileNotFoundError
        wf_data = read_wf(input_name)

        # get initial guess
        primary_number = get_primary_number(element_number, shell, is_excited)
        ini_guess = get_initial_gbs_guess(element_number, primary_number, shell_name)
        # fitting by Gaussian
        Ae, ae, Y_fit, error = fit_wf(wf_data.T, ini_guess, l=shell)
        if args.plot:
            output_fig = "{:03d}.wf-{}{}-fitting.png".format(element_number, shell_name, is_excited)
            plot_fitting(wf_data.T, shell, Ae, ae, output=output_fig)

        output = "{:03d}.wf-{}{}.gbs".format(element_number, shell_name, is_excited)
        with open(output, 'w') as f:
            for A, a in zip(Ae, ae):
                f.write("{: .8E}  {: .8E}\n".format(A, a))


def read_gaussian(element_number, shell, is_excited, Fdata_path=None):
    shell_name = SHELL_NAME[shell]
    filename = os.path.join(Fdata_path, "{:03d}.wf-{}{}.gbs".format(element_number, shell_name, is_excited))
    with open(filename, 'r') as f:
        lines = f.readlines()
    array = np.asarray([[float(i) for i in line.split()] for line in lines])
    return {'alpha': array[:, 0], 'coefficient': array[:, 1], 'degree': len(array)}
