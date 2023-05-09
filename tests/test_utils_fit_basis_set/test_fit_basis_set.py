from thunder_ase.utils.basis_set import read_wf, loss_function, plot_fitting, loss_jac
import numpy as np
from scipy.optimize import basinhopping, minimize


wf_file = '/home/ren/Programs/Github/thunder-ase/data/Si_thunder_begin/014.wf-p1.dat'

data = read_wf(wf_file).T
R, Y = np.asarray(data)
nz = 10
l = 1

init_alpha = np.random.random(nz) * 3.5 - 1.0  # -1.0 ~ 2.5
init_coeff = np.random.random(nz) * 0.5 + 0.5  # 0.5 ~ 1.0
init_guess = np.concatenate([init_alpha, init_coeff])

bnds_lb = [-1.5] * nz + [-1.5] * nz
bnds_ub = [4.0] * nz + [1.5] * nz
bnds = list(zip(bnds_lb, bnds_ub))

res = basinhopping(loss_function, init_guess, T=0.0001, niter=100, disp=True, stepsize=0.5,
                   niter_success=50,
                   interval=20,
                   stepwise_factor=0.8,
                   minimizer_kwargs={'args': (l, R, Y),
                                     'bounds': bnds,
                                     'jac': loss_jac,
                                     }
                   )

alpha = 10 ** res.x[0:nz]  # exponential parameters,
Ae = res.x[nz:]  # coefficients

plot_fitting(data, l, Ae, alpha)

print("Done!")
