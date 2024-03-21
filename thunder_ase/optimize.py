"""
Due to the FIREBALL use Harris functional, the fmax as convergence criteria is too strict.
Herein, we rewrite some of the optimizer by using force rms as criteria.
"""
from os.path import isfile
from types import MethodType

import numpy as np
from ase.io.jsonio import read_json, write_json
from ase.optimize import (MDMin, FIRE, LBFGS, LBFGSLineSearch, BFGSLineSearch, BFGS,
                          GoodOldQuasiNewton, QuasiNewton, GPMin, Berny, ODE12r)
from ase.optimize.optimize import Optimizer, RestartError
from ase.parallel import world


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


# This is a new optimizer using RMS converge

for opter in (
        MDMin, FIRE, LBFGS, LBFGSLineSearch, BFGSLineSearch, BFGS, GoodOldQuasiNewton, QuasiNewton, GPMin, Berny,
        ODE12r):
    opter.converged = MethodType(rms_converged, opter)


class RattleBFGS(BFGS):

    """
    Since the force from Harris functional is not correct, we can mimic the quantum movement by
    using ase.Atoms.rattle, and then average the force.
    """
    def __init__(self, atoms, restart=None, logfile='-', trajectory=None,
                 maxstep=None, master=None, alpha=None, rattle=3):

        self.rattle = rattle

        BFGS.__init__(self, atoms, restart=restart, logfile=logfile, trajectory=trajectory,
                      maxstep=maxstep, master=master, alpha=alpha)

    def initialize(self):
        BFGS.initialize(self)
        if not isfile('force.rattle'):
            self.dump_f0(0.0)

    def todict(self):
        d = BFGS.todict(self)
        d.update(rattle=self.rattle)
        return d

    def get_forces(self):
        # read old force from file
        f0 = self.read_f0()
        # read new force from atoms
        f1 = self.atoms.get_forces() / self.rattle
        f = f0 + f1
        # write new force to file
        self.dump_f0(f)
        return f

    def dump_f0(self, data):
        if world.rank == 0 and self.rattle != 1:
            with open('force.rattle', 'w') as fd:
                write_json(fd, data)

    def read_f0(self):
        with open('force.rattle') as fd:
            try:
                f0 = read_json(fd, always_array=True)
            except Exception as ex:
                msg = ('Could not decode temp force file as JSON.  '
                       f'You may need to delete the restart file force.rattle')
                raise RestartError(msg) from ex
        return f0

    def step(self, f=None):
        atoms = self.atoms

        if f is None:
            f = self.get_forces()

        if self.nsteps % self.rattle == 0:
            # run BFGS step
            BFGS.step(self, f)
            # write 0.0 to new force.rattle
            self.dump_f0(0.0)
        else:
            # rattle atoms position with a very small displacement
            atoms.rattle(0.00001)

