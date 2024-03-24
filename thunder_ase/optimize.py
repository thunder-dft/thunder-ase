"""
Due to the FIREBALL use Harris functional, the fmax as convergence criteria is too strict.
Herein, we rewrite some of the optimizer by using force rms as criteria.
Furthermore, we rewrite the step() function for MDMin, BFGS, LBFGS, and FIRE, where the rotation and
displacement of the whole molecule were removed.

NOTE: We encourage to use MDMin, BFGS, LBFGS, FIRE as optimizer to instead of original ase optimizer.
"""
import numpy as np
from ase.optimize import MDMin as _MDMin
from ase.optimize import LBFGS as _LBFGS
from ase.optimize import BFGS as _BFGS
from ase.optimize import FIRE as _FIRE
from ase.md.velocitydistribution import ZeroRotation, Stationary
from ase.optimize.optimize import Optimizer
from numpy.linalg import eigh


def rms_converged(self, forces=None):
    """Did the optimization converge?"""
    if forces is None:
        forces = self.atoms.get_forces()
    if hasattr(self.atoms, "get_curvature"):
        return (forces ** 2).sum(
            axis=1
        ).max() < self.fmax ** 2 and self.atoms.get_curvature() < 0.0

    fmax2 = (forces ** 2).sum(axis=1).max()
    rms2 = (forces ** 2).sum(axis=1).mean()
    print("fmax: {:.5f}, force rms: {:.5f}".format(np.sqrt(fmax2), np.sqrt(rms2)))
    return rms2 < self.fmax ** 2


def zero_rotation(atoms, v):
    atoms_tmp = atoms.copy()
    atoms_tmp.set_velocities(v)
    Stationary(atoms_tmp)
    ZeroRotation(atoms_tmp)
    return atoms_tmp.get_velocities()


class MDMin(_MDMin):
    defaults = {**Optimizer.defaults, 'dt': 0.1}

    def step(self, f=None):
        atoms = self.atoms

        if f is None:
            f = atoms.get_forces()

        if self.v is None:
            self.v = np.zeros((len(atoms), 3))
        else:
            self.v += 0.5 * self.dt * f
            # Correct velocities:
            vf = np.vdot(self.v, f)
            if vf < 0.0:
                self.v[:] = 0.0
            else:
                self.v[:] = f * vf / np.vdot(f, f)

        self.v += 0.5 * self.dt * f
        # set displacement and angular moment to zero
        self.v = zero_rotation(atoms, self.v)

        r = atoms.get_positions()
        atoms.set_positions(r + self.dt * self.v)
        self.dump((self.v, self.dt))


class BFGS(_BFGS):
    def step(self, f=None):
        atoms = self.atoms

        if f is None:
            f = atoms.get_forces()

        r = atoms.get_positions()
        f = f.reshape(-1)
        self.update(r.flat, f, self.r0, self.f0)
        omega, V = eigh(self.H)

        dr = np.dot(V, np.dot(f, V) / np.fabs(omega)).reshape((-1, 3))
        steplengths = (dr ** 2).sum(1) ** 0.5
        dr = self.determine_step(dr, steplengths)
        # set displacement and angular moment to zero
        dr = zero_rotation(atoms, dr)
        atoms.set_positions(r + dr)
        self.r0 = r.flat.copy()
        self.f0 = f.copy()
        self.dump((self.H, self.r0, self.f0, self.maxstep))


class LBFGS(_LBFGS):
    def step(self, f=None):
        """Take a single step

        Use the given forces, update the history and calculate the next step --
        then take it"""

        if f is None:
            f = self.atoms.get_forces()

        r = self.atoms.get_positions()

        self.update(r, f, self.r0, self.f0)

        s = self.s
        y = self.y
        rho = self.rho
        H0 = self.H0

        loopmax = np.min([self.memory, self.iteration])
        a = np.empty((loopmax,), dtype=np.float64)

        # ## The algorithm itself:
        q = -f.reshape(-1)
        for i in range(loopmax - 1, -1, -1):
            a[i] = rho[i] * np.dot(s[i], q)
            q -= a[i] * y[i]
        z = H0 * q

        for i in range(loopmax):
            b = rho[i] * np.dot(y[i], z)
            z += s[i] * (a[i] - b)

        self.p = - z.reshape((-1, 3))

        g = -f
        if self.use_line_search is True:
            e = self.func(r)
            self.line_search(r, g, e)
            dr = (self.alpha_k * self.p).reshape(len(self.atoms), -1)
        else:
            self.force_calls += 1
            self.function_calls += 1
            dr = self.determine_step(self.p) * self.damping

        # set displacement and angular moment to zero
        dr = zero_rotation(self.atoms, dr)

        self.atoms.set_positions(r + dr)

        self.iteration += 1
        self.r0 = r
        self.f0 = -g
        self.dump((self.iteration, self.s, self.y,
                   self.rho, self.r0, self.f0, self.e0, self.task))


class FIRE(_FIRE):
    def step(self, f=None):
        atoms = self.atoms

        if f is None:
            f = atoms.get_forces()

        if self.v is None:
            self.v = np.zeros((len(atoms), 3))
            if self.downhill_check:
                self.e_last = atoms.get_potential_energy(
                    force_consistent=self.force_consistent)
                self.r_last = atoms.get_positions().copy()
                self.v_last = self.v.copy()
        else:
            is_uphill = False
            if self.downhill_check:
                e = atoms.get_potential_energy(
                    force_consistent=self.force_consistent)
                # Check if the energy actually decreased
                if e > self.e_last:
                    # If not, reset to old positions...
                    if self.position_reset_callback is not None:
                        self.position_reset_callback(atoms, self.r_last, e,
                                                     self.e_last)
                    atoms.set_positions(self.r_last)
                    is_uphill = True
                self.e_last = atoms.get_potential_energy(
                    force_consistent=self.force_consistent)
                self.r_last = atoms.get_positions().copy()
                self.v_last = self.v.copy()

            vf = np.vdot(f, self.v)
            if vf > 0.0 and not is_uphill:
                self.v = (1.0 - self.a) * self.v + self.a * f / np.sqrt(
                    np.vdot(f, f)) * np.sqrt(np.vdot(self.v, self.v))
                if self.Nsteps > self.Nmin:
                    self.dt = min(self.dt * self.finc, self.dtmax)
                    self.a *= self.fa
                self.Nsteps += 1
            else:
                self.v[:] *= 0.0
                self.a = self.astart
                self.dt *= self.fdec
                self.Nsteps = 0

        self.v += self.dt * f

        # set displacement and angular moment to zero
        self.v = zero_rotation(self.atoms, self.v)

        dr = self.dt * self.v
        normdr = np.sqrt(np.vdot(dr, dr))
        if normdr > self.maxstep:
            dr = self.maxstep * dr / normdr
        r = atoms.get_positions()
        atoms.set_positions(r + dr)
        self.dump((self.v, self.dt))
