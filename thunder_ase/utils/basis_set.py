import argparse
import os.path

import numpy as np
from scipy.optimize import minimize, basinhopping
import matplotlib.pyplot as plt

from thunder_ase.utils import ordinal
from thunder_ase.utils.shell_dict import SHELL_NUM, SHELL_NAME
from ase.data import chemical_symbols
from ase.units import Bohr


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


def norm_c(a, shell=0):
    # normalization constant for pure Gaussian basis functions
    # see https://iodata.readthedocs.io/en/latest/basis.html
    norm = (2 * a / np.pi) ** 0.75  # s orbital
    if shell == 1:  # p orbital
        norm = norm * 2 * np.sqrt(a)
    elif shell == 2:  # d orbital
        norm = norm * 4 * a / np.sqrt(3)
    elif shell > 2:
        raise NotImplementedError
    return norm


def dnorm_c(a, shell=0):
    # derivative of constant part of normalization constant
    dnorm = (2 * a / np.pi) ** (-0.25) * (2 / np.pi)  # s orbital
    if shell == 1:  # p orbital
        dnorm = dnorm * 2 * np.sqrt(a) + norm_c(a, 0) * 2 / np.sqrt(a)
    elif shell == 2:  # d orbital
        dnorm = dnorm * 4 * a / np.sqrt(3) + norm_c(a, 0) * 4 / np.sqrt(3)
    elif shell > 2:
        raise NotImplementedError
    return dnorm


def gaussian(r, l=0, A=np.array([1.0]), a=np.array([0.1])):
    # gaussian function
    f = np.sum([Ai * (r ** l) * np.exp(-ai * (r ** 2)) * norm_c(ai, l)
                for Ai, ai in zip(A, a)], axis=0)
    return f


# loss function
def loss_function(x, *args):
    """
    TODO: use sum( (4 * pi * r**2 * (wf**2 - wf0**2))**2 ) as loss function
    :param x:
    :param args:
    :return:
    """
    len_x = int(len(x) / 2)
    alpha = 10 ** x[0:len_x]
    A = np.asarray(x[len_x:])
    l, r, Y = args
    if r.shape != Y.shape:
        raise ValueError('Shapes for args[1] and args[2] are not equal!')

    loss = np.mean((gaussian(r, l, A, alpha) - Y) ** 2)
    return loss


def loss_jac(x, *args):
    len_x = int(len(x) / 2)
    alpha = 10 ** x[0:len_x]
    A = np.asarray(x[len_x:])
    l, r, Y = args
    if r.shape != Y.shape:
        raise ValueError('Shapes for args[1] and args[2] are not equal!')
    d_alpha = np.asarray([Ai * (r ** l) * np.exp(-ai * (r ** 2)) * ai * np.log(10) *
                          (dnorm_c(ai, l) - r ** 2 * norm_c(ai, l)) for Ai, ai in zip(A, alpha)])
    d_A = np.asarray([(r ** l) * np.exp(-ai * (r ** 2)) * norm_c(ai, l) for ai in alpha])
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
    error = res.fun
    return [Ae, alpha, error]


def fit_wf_from_random(data, l=0, tol=1e-5, Nzeta0=3, Nzeta_max=10, Niter=3, bnds=None):
    """

    :param Niter:
    :param data: shape = [N, 2]
    :param l: principal quantum number
    :param tol: Error tolerance
    :param Nzeta: number of gaussian
    :param bnds: boundary for fitting
    :return:
    """
    R, Y = np.asarray(data)
    error = np.inf
    result_res = None

    print("Fitting parameters: shell = {}, Nzeta0 = {}, Nzeta_max = {}, Ntry = {}, tol = {}"
          .format(l, Nzeta0, Nzeta_max, Niter, tol))

    for nz in range(Nzeta0, Nzeta_max + 1):
        if bnds is None:
            bnds_lb = [-1.5] * nz + [-1.0] * nz
            bnds_ub = [4.0] * nz + [1.0] * nz
            bounds = list(zip(bnds_lb, bnds_ub))
        else:
            bounds = bnds
        for itry in range(Niter):
            print("::: The {} try for Nzeta = {} ...".format(ordinal(itry + 1), nz))
            init_alpha = np.random.random(nz) * 3.5 - 1.0  # -1.0 ~ 2.5
            init_coeff = np.random.random(nz) * 0.5 + 0.5  # 0.5 ~ 1.0
            init_guess = np.concatenate([init_alpha, init_coeff])
            res = basinhopping(loss_function, init_guess, T=0.0001, niter=100, disp=False, stepsize=0.5,
                               niter_success=50,
                               interval=20,
                               stepwise_factor=0.8,
                               minimizer_kwargs={'args': (l, R, Y),
                                                 'bounds': bounds,
                                                 'jac': loss_jac,
                                                 }
                               )
            if res.fun < error:
                result_res = res
                error = res.fun
                if error < tol:
                    alpha = 10 ** result_res.x[0:nz]  # exponential parameters,
                    Ae = result_res.x[nz:]  # coefficients
                    print("Success!: Fitting error {} meet the tolerance {} for {} gaussians."
                          .format(error, tol, nz))
                    return [Ae, alpha, error]
        print("Fitting error {} didn't meet the tolerance {} for {} gaussians after {} try."
              .format(error, tol, nz, Niter))

    alpha = 10 ** result_res.x[0:Nzeta_max]  # exponential parameters,
    Ae = result_res.x[Nzeta_max:]  # coefficients
    print("Warning: Fitting error {} didn't meet the tolerance {} for {} gaussians."
          .format(error, tol, Nzeta_max))

    return [Ae, alpha, error]


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


def expand_data(wf_data):
    """
    Expand the cutoff to current_cutoff + 3.0 angstrom.
    This is to get better fitting for the gaussian tail.
    :param wf_data:
    :return:
    """
    r, y = wf_data
    r_max = r.max()
    dr = r_max / len(r)
    r_plus = 3.0
    r_extend = np.arange(r_max + dr, r_max + r_plus, dr)
    new_r = np.concatenate([r, r_extend])
    new_y = np.concatenate([y, np.zeros(len(r_extend))])

    return np.asarray([new_r, new_y])


# run this after begin.x
def fit_gaussian(prog='fit_gaussians',
                 description='Fit Fireball basis set to Gaussian-type basis set.'):
    parser = argparse.ArgumentParser(
        prog=prog,
        description=description, )
    parser.add_argument('input_name', nargs='+', help='Fireball wave function file.')
    parser.add_argument('-p', '--plot', action='store_true')
    parser.add_argument('-t', '--tolerance', type=float, dest='tolerance', default=1e-4)
    parser.add_argument('-Nz', '--Nzeta_max', type=int, dest='nzeta_max', default=10)
    parser.add_argument('-Nz0', '--Nzeta0', type=int, dest='nzeta0', default=4)
    parser.add_argument('-Nt', '--Ntry', type=int, dest='ntry', default=3)
    args = parser.parse_args()

    for input_name in args.input_name:
        name_list = input_name.split('.')  # format: '001.wf-s0.dat', element_number, shell, is_excited
        element_number = int(name_list[0])
        shell_name, is_excited = name_list[1][-2:]
        shell = SHELL_NUM[shell_name]
        if not os.path.exists(input_name):
            print("{} doesn't exist!".format(input_name))
            raise FileNotFoundError
        wf_data = read_wf(input_name).T

        coeff_angular = Bohr**1.5
        # normalize wf_data for different l: sqrt((2*l+1) / (4*pi))
        if shell == 0:
            coeff_angular = coeff_angular * np.sqrt(1.0 / (4.0 * np.pi))
        elif shell == 1:
            coeff_angular = coeff_angular * np.sqrt(3.0 / (4.0 * np.pi))
        elif shell == 2:
            coeff_angular = coeff_angular * np.sqrt(5.0 / (4.0 * np.pi))
        else:
            NotImplementedError

        wf_data[1] = wf_data[1] * coeff_angular
        wf_data = expand_data(wf_data)
        # fitting from random
        Ae, ae, error = fit_wf_from_random(data=wf_data, l=shell,
                                           tol=args.tolerance,
                                           Nzeta_max=args.nzeta_max,
                                           Nzeta0=args.nzeta0,
                                           Niter=args.ntry)

        if args.plot:
            output_fig = "{:03d}.wf-{}{}-fitting.png".format(element_number, shell_name, is_excited)
            plot_fitting(wf_data, shell, Ae, ae, output=output_fig)

        output = "{:03d}.wf-{}{}.gbs".format(element_number, shell_name, is_excited)
        with open(output, 'w') as f:
            for A, a in zip(Ae, ae):
                f.write("{: .8E}  {: .8E}\n".format(a, A))  # order is alpha, co_effi


def read_gaussian(element_number, shell, is_excited, Fdata_path=None):
    shell_name = SHELL_NAME[shell]
    filename = os.path.join(Fdata_path, 'basis', "{:03d}.wf-{}{}.gbs".format(element_number, shell_name, is_excited))
    with open(filename, 'r') as f:
        lines = f.readlines()
    array = np.asarray([[float(i) for i in line.split()] for line in lines])
    return {'alpha': array[:, 0], 'coefficient': array[:, 1], 'degree': len(array)}
