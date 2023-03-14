import numpy as np
from scipy.optimize import minimize
import matplotlib.pyplot as plt


# read wf file
def read_wf(filename):
    with open(filename, 'r') as f:
        lines = f.readlines()
    data = [list(map(float, line.strip().split())) for line in lines if len(line.strip().split()) == 2]
    return np.asarray(data)


# gaussian function
def gaussian(r, n=1, A=np.array([1]), a=np.array([1])):
    f = np.sum([Ai * r ** (n - 1) * np.exp(-ai * r ** 2)
                for Ai, ai in zip(A, a)], axis=0)
    return f


# loss function
def loss_function(x, *args):
    len_x = int(len(x) / 2)
    A = x[0:len_x]
    a = x[len_x:]
    n, r, Y = args
    if r.shape != Y.shape:
        raise ValueError('Shapes for args[1] and args[2] are not equal!')
    loss = np.abs(Y - gaussian(r, n, A, a)).sum()
    return loss


if "__name__" == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        prog='FitBasisSet',
        description='Fit Fireball basis set to Gaussian-type basis set.',)

# main
# filename = '../tests/test_utils_fit_basis_set/001.wf-s0.dat'
filename = '../tests/test_utils_fit_basis_set/007.wf-p1.dat'

n = 1
A0 = [-1.0, 1.0, 2.0, 3.0]
a0 = [1.1, 1.0, 2.0, 3.0]
# bnds = [(0, None)] * len(A0) + [(0, None)] * len(a0)
bnds = None

data = read_wf(filename)
R, Y = data.T
x0 = np.asarray(A0 + a0)

res = minimize(loss_function, x0, bounds=bnds, args=(n, R, Y))

print(res.message)

# plot result
len_x = int(len(x0) / 2)
Ae = res.x[0:len_x]
ae = res.x[len_x:]
plt.plot(R, Y, '-')
plt.plot(R, gaussian(R, n, Ae, ae))
for Ai, ai in zip(Ae, ae):
    plt.plot(R, gaussian(R, n, [Ai], [ai]))
plt.show()
print("Done")
