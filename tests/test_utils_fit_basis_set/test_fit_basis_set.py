from thunder_ase.utils.basis_set import read_wf, loss_function, plot_fitting, loss_jac, fit_wf_from_random
import numpy as np
from scipy.optimize import basinhopping, minimize


wf_file = '/home/ren/Programs/Github/thunder-ase/tests/test_utils_fit_basis_set/007.wf-p1.dat'

data = read_wf(wf_file).T
nzeta_max = 15
ntry = 3
tolerance = 1e-5
shell = 1

Ae, ae, error = fit_wf_from_random(data=data, l=shell,
                                          tol=tolerance, Nzeta_max=nzeta_max, Niter=ntry, Nzeta0=15)
plot_fitting(data, shell, Ae, ae)

print("Done!")
