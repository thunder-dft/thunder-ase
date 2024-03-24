import timeit
from ase import units
from ase.build import molecule
from ase.calculators.socketio import SocketIOCalculator
from ase.md.verlet import VelocityVerlet as NVE
from thunder_ase.fireball import Fireball


def run_with_socket():
    # Construct Structure
    atoms = molecule('C60')

    # set Fdata dir
    Fdata_path = 'Fdata'

    # Sockets
    unixsocket = 'thunder-ase'

    max_step = 10
    kwargs = {
        'ipi': 1,
        'nstepf': max_step + 1,  # max step
        'inet': 0,
        'host': unixsocket,
        'scf_tolerance_set': 1E-8,
    }
    fireball = Fireball(command='fireball-ase.9.x', Fdata_path=Fdata_path, **kwargs)
    dyn = NVE(atoms, trajectory='md.traj', logfile='md.log', timestep=0.5 * units.fs)

    with SocketIOCalculator(fireball, log=None, unixsocket=unixsocket) as calc:
        atoms.calc = calc
        dyn.run(max_step)


print(timeit.timeit(run_with_socket, number=1))
