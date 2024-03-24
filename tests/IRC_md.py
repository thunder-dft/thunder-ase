import numpy as np
from scan_ts import FixCombBond
import ase.io
from ase import units
from ase.md.nvtberendsen import NVTBerendsen as NVT
from ase.md.velocitydistribution import ZeroRotation
from ase.optimize import BFGS
from ase.calculators.socketio import SocketIOCalculator
from thunder_ase.fireball import Fireball
import os

filename = 'Bu_55_ren.xyz'
atoms0 = ase.io.read(filename)
atoms0.positions += 10.0

init_length = 3.0
final_length = 7.0
delta_l = 0.01
nl = int(abs(init_length - final_length) / delta_l)
xx = np.linspace(init_length, final_length, nl)

# Run MD 0.25 * 800 = 0.2 ps
max_step = 800
dt = 0.25  # timestep

# FIREBALL settings
# set Fdata dir
Fdata_path = 'Fdata'
# Sockets
unixsocket = 'thunder-ase16'
kwargs = {
    'ipi': 1,
    'nstepf': max_step + 1,  # max step
    'inet': 0,
    'host': unixsocket,
    'scf_tolerance_set': 1E-8,
    'iwriteout_cdcoeffs': 1,
    'iwriteout_charges': 1,
}
fireball = Fireball(command='fireball-ase.9.x', Fdata_path=Fdata_path, **kwargs)
fireball.prefix = ''
atoms = atoms0.copy()

# optimize structure
dyn = BFGS(atoms)

with SocketIOCalculator(fireball, log=None, unixsocket=unixsocket) as calc:
    atoms.calc = calc
    dyn.run(fmax=0.2)

ase.io.write('Bu_55_opt.xyz', atoms)
# atoms.cell = [20, 20, 20]

for idx, il in enumerate(xx):
    # Define reaction coordinate
    liste = []
    liste += [[9, 3, 1.0]]  # bond 1
    liste += [[10, 2, 1.0]]  # bond 2
    cc = FixCombBond(il, liste)
    cc.initialize(atoms)
    # adjust structure according to the initial value (in this example the initial value is not modified).
    atoms.set_positions(cc.xyz)
    # set all constraints, keeping old constraints
    atoms.set_constraint([cc] + atoms.constraints)

    workdir = '{:0>3d}'.format(idx)
    if not os.path.isdir(workdir):
        os.makedirs(workdir)
    output = workdir + '/md-{:.3f}'.format(il)
    dyn = NVT(atoms,
              temperature_K=230,
              taut=10 * dt * units.fs,
              logfile=output + '.log',
              timestep=dt * units.fs,
              )

    with SocketIOCalculator(fireball, log=None, unixsocket=unixsocket) as calc:
        atoms.calc = calc
        for istep, _ in enumerate(dyn.irun(max_step)):
            ZeroRotation(atoms, preserve_temperature=True)  # Sets the total angular momentum to zero
            if (istep + 1) % 30 == 0:
                ase.io.write(output + ".xyz", atoms, append=True)
            if (istep + 1) % 100 == 0:
                calc.launch_client.calc.read_results()
                mwfn_name = output + "-{:0>5d}.mwfn".format(istep)
                calc.launch_client.calc.write_mwfn(filename=mwfn_name)
                # write multiwfn input for LUMO and HOMO
                for orbital in ('h', 'l'):
                    multiwfn_input = [
                        '5',
                        '4',
                        orbital,
                        '2',
                        '2',
                        '0',
                        'q',
                    ]
                    input_name = output + "-{:0>5d}-{}.inp".format(istep, orbital)
                    with open(input_name, 'w') as f:
                        f.write('\n'.join(multiwfn_input))

                    print('Calculate HOMO/LUMO functions ...')
                    os.system('Multiwfn {} < {} > /dev/null'.format(mwfn_name, input_name))
                    os.rename('MOvalue.cub', output + "-{:0>5d}-{}.cub".format(istep, orbital))

    if os.path.exists('001.xyz'):
        os.remove('001.xyz')
    if os.path.exists('001.json'):
        os.remove('001.json')
    if os.path.exists('001.log'):
        os.remove('001.log')
    atoms.set_constraint()
    atoms.calc = None
