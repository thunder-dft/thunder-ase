"""
Due to the FIREBALL use Harris functional, the fmax as convergence criteria is too strict.
Herein, we rewrite some of the optimizer by using force rms as criteria.
"""
from types import MethodType

import numpy as np
from ase.optimize import BFGS


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


BFGS.converged = MethodType(rms_converged, BFGS)
