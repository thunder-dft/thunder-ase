import argparse
import os.path

import numpy as np
from numpy.ma.core import around
from scipy.optimize import minimize, basinhopping
from scipy.integrate import quad  # gaussian quadrature.
import matplotlib.pyplot as plt
from thunder_ase.utils.shell_dict import SHELL_NUM, SHELL_NAME
from ase.units import Bohr


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
    """
    Gaussian function
    :param r: radial distance
    :param l: angular momentum
    :param A: coefficients
    :param a: exponential parameters
    :return: 
    """
    f = np.sum([Ai * (r ** l) * np.exp(-ai * (r ** 2)) * norm_c(ai, l)
                for Ai, ai in zip(A, a)], axis=0)
    return f

def sphere_rho(r, l, A, a):
    # function for integral of fitted density over the whole space
    f = 4 * np.pi * r**2 * (gaussian(r, l, A, a))**2
    return f

# loss function
def loss_function(x, *args):
    """
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


    # Fitting type: Multiwfn manual 3.300.2 as reference
    #rho0 = Y ** 2
    #rho = gaussian(r, l, A, alpha) ** 2
    # type 1: absolute error of density: sum((rho - rho0)**2)
    #loss = np.sum((rho - rho0) ** 2)

    # type 2: relative error of density: sum(((rho - rho0)/rho)**2)
    #loss = np.sum(((rho - rho0) / rho) ** 2)

    # type 3: error of density radial distribution function: sum((4 * pi * r**2 * (rho - rho0))**2)
    #loss = np.sum((4 * np.pi * r**2 * (rho - rho0))**2)

    # type 4: mean error of wave function: mean((wf-wf0)**2)  # not in multiwfn, my version
    loss = np.mean((gaussian(r, l, A, alpha) - Y) ** 2)

    # "Commonly, fitting type (3) is preferred over others, because the density fitted in this way can
    # reproduce actual density over the entire range. Fitted density corresponding to fitting type (1) can
    # only well represent the region very close to nucleus, since electron density in this region is
    # significantly larger than other regions. The density fitted by type (3) usually represents valence and
    # tail regions well, but has a poor description in the region very close to nucleus."

    return loss


def loss_jac(x, *args):
    len_x = int(len(x) / 2)
    alpha = 10 ** x[0:len_x]
    A = np.asarray(x[len_x:])
    l, r, Y = args
    rho0 = Y ** 2
    rho = gaussian(r, l, A, alpha) ** 2

    if r.shape != Y.shape:
        raise ValueError('Shapes for args[1] and args[2] are not equal!')
    d_alpha = np.asarray([Ai * (r ** l) * np.exp(-ai * (r ** 2)) * ai * np.log(10) *
                          (dnorm_c(ai, l) - r ** 2 * norm_c(ai, l)) for Ai, ai in zip(A, alpha)])
    d_A = np.asarray([(r ** l) * np.exp(-ai * (r ** 2)) * norm_c(ai, l) for ai in alpha])
    der = np.concatenate([d_alpha, d_A])
    # type 3
    #k = 64 * np.pi**2 * r**4
    #jac = np.sum(k * (rho - rho0) * gaussian(r, l, A, alpha) * der, axis=1)

    # type 4
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

    rho0 = Y ** 2
    dr = R.max() / len(R)
    Nele = np.sum(4 * np.pi * R**2 * rho0 * dr)

    # "Integral of fitted density over the whole space may deviate from actual number of electrons
    # (Nelec), clearly this breaks physical meaning of fitted density and makes it useless in many scenarios.
    # In order to solve this problem, the fitted coefficients should be scaled by a factor"
    Ae = Ae * Nele / quad(sphere_rho, a=0, b=10, args=(l, Ae, alpha))[0]

    return [Ae, alpha, error]


def fit_wf_from_random(data, l=0, tol=1e-5, Nzeta0=3, Nzeta_max=10, Niter=3, bnds=None):
    """

    :param Niter:
    :param data: shape = [N, 2]
    :param l: principal quantum number
    :param tol: Error tolerance
    :param Nzeta0: number of initial gaussian
    :param Nzeta_max: maximum number of gaussian
    :param bnds: boundary for fitting
    :return:
    """
    R, Y = np.asarray(data)
    error = np.inf
    result_res = None
    rho0 = Y ** 2
    dr = R.max() / len(R)
    Nele = np.sum(4 * np.pi * R ** 2 * rho0 * dr)
    print(f"Number of electrons: {Nele} ~ {around(Nele, 0)}")
    Nele = around(Nele, 0)
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

                    # "Integral of fitted density over the whole space may deviate from actual number of electrons
                    # (Nelec), clearly this breaks physical meaning of fitted density and makes it useless in many scenarios.
                    # In order to solve this problem, the fitted coefficients should be scaled by a factor"
                    Ae = Ae * Nele / quad(sphere_rho, a=0, b=10, args=(l, Ae, alpha))[0]
                    return [Ae, alpha, error]
        print("Fitting error {} didn't meet the tolerance {} for {} gaussians after {} try."
              .format(error, tol, nz, Niter))

    alpha = 10 ** result_res.x[0:Nzeta_max]  # exponential parameters,
    Ae = result_res.x[Nzeta_max:]  # coefficients
    print("Warning: Fitting error {} didn't meet the tolerance {} for {} gaussians."
          .format(error, tol, Nzeta_max))
    Ae = Ae * Nele / quad(sphere_rho, a=0, b=10, args=(l, Ae, alpha))[0]
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
        plt.plot(R, gaussian(R, l, np.array([Ai]), np.array([ai])))
    if output is not None:
        plt.savefig(output)
    else:
        plt.show()


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


def angular_coeff(shell):
    coeff_angular = Bohr ** 1.5
    # normalize wf_data for different l: sqrt((2*l+1) / (4*pi))
    if shell == 0:
        coeff_angular = coeff_angular * np.sqrt(1.0 / (4.0 * np.pi))
    elif shell == 1:
        coeff_angular = coeff_angular * np.sqrt(3.0 / (4.0 * np.pi))
    elif shell == 2:
        coeff_angular = coeff_angular * np.sqrt(5.0 / (4.0 * np.pi))
    else:
        raise NotImplementedError
    return coeff_angular


def get_full_wf(input_name):
    name_list = input_name.split('.')  # format: '001.wf-s0.dat', element_number, shell, is_excited
    shell_name, is_excited = name_list[1][-2:]
    shell = SHELL_NUM[shell_name]
    if not os.path.exists(input_name):
        print("{} doesn't exist!".format(input_name))
        raise FileNotFoundError
    wf_data = read_wf(input_name).T
    wf_data[1] = wf_data[1] * angular_coeff(shell)
    wf_data = expand_data(wf_data)
    return wf_data


# run this after begin.x
def fit_gaussian(parser=None, **kwargs):
    """
    :param parser: argparse.ArgumentParser
    :param kwargs:
        input_name: input file name, can be a list or str.
        nzeta_max: maximum number of Nzeta, default 10
        nzeta0: initial Nzeta, default 3
        ntry: number of trying with random initial guest, default 4
        tolerance: fitting error tolerance, default 0.00001
    :return:
    """
    if parser is not None:
        args = parser.parse_args()
    else:
        default = {
            'tolerance': 1e-5,
            'nzeta_max': 10,
            'nzeta0': 3,
            'ntry': 4,
        }
        default.update(kwargs)

        # convert dict to a AttributeDict object
        class AttributeDict(dict):
            __getattr__ = dict.__getitem__
            __setattr__ = dict.__setitem__
            __delattr__ = dict.__delitem__
        args = AttributeDict(default)
        if type(args.input_name) == str:
            args.input_name = [args.input_name]

    for input_name in args.input_name:
        wf_data = get_full_wf(input_name)
        name_list = input_name.split('.')  # format: '001.wf-s0.dat', element_number, shell, is_excited
        element_number = int(name_list[0])
        shell_name, is_excited = name_list[1][-2:]
        shell = SHELL_NUM[shell_name]
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


# run this after begin.x
def fit_gaussian_command(prog='fit-gaussians',
                         description='Fit Fireball basis set to Gaussian-type basis set.'):
    parser = argparse.ArgumentParser(
        prog=prog,
        description=description, )
    parser.add_argument('input_name', nargs='+', help='Fireball wave function file.')
    parser.add_argument('-p', '--plot', action='store_true')
    parser.add_argument('-t', '--tolerance', type=float, dest='tolerance', default=1e-5)
    parser.add_argument('-Nz', '--Nzeta_max', type=int, dest='nzeta_max', default=10)
    parser.add_argument('-Nz0', '--Nzeta0', type=int, dest='nzeta0', default=3)
    parser.add_argument('-Nt', '--Ntry', type=int, dest='ntry', default=4)
    fit_gaussian(parser=parser)


def read_gaussian(element_number, shell, is_excited, Fdata_path=None):
    shell_name = SHELL_NAME[shell]
    filename = os.path.join(Fdata_path, 'basis', "{:03d}.wf-{}{}.gbs".format(element_number, shell_name, is_excited))
    with open(filename, 'r') as f:
        lines = f.readlines()
    array = np.asarray([[float(i) for i in line.split()] for line in lines])
    return {'alpha': array[:, 0], 'coefficient': array[:, 1], 'degree': len(array)}
