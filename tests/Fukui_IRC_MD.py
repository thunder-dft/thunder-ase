import numpy as np
from thunder_ase.fireball import Fireball
import ase.io
import os

init_length = 3.0
final_length = 7.0
delta_l = 0.01
nl = int(abs(init_length - final_length) / delta_l)
xx = np.linspace(init_length, final_length, nl)

Fdata_path = '../Fdata'

for idx, il in enumerate(xx):
    workdir = '{:0>3d}'.format(idx)
    os.chdir(workdir)
    xyzfile = 'md-{:.3f}.xyz'.format(il)

    for qstate in [0, -1, 1]:
        atoms = ase.io.read(xyzfile, index='-1')
        if os.path.isfile('001.CHARGES'):
            os.remove('001.CHARGES')

        kwargs = {
            'iconstraint_l': 1,
            'iconstraint_rcm': 0,  # don't shift molecule to COM
            'iwriteout_charges': 1,  # Writing out the charges.
            'iwriteout_cdcoeffs': 1,
            'efermi_T': 200.0,
            'max_scf_iterations_set': 100,
            'scf_tolerance_set': 0.00000001,
            'beta_set': 0.04,
            'qstate': qstate,
        }

        fireball = Fireball(command='fireball-ase.9.x', Fdata_path=Fdata_path, **kwargs)
        atoms.calc = fireball
        e0 = atoms.get_potential_energy()
        print("Calculate struture {} with {} electron, energy is {}.".format(idx, qstate, e0))
        fireball.write_mwfn(filename="md-{:.3f}_q{:d}.mwfn".format(il, qstate))

    # write multiwfn input
    multiwfn_input = [
        '22',
        '3',
        'md-{:.3f}_q0.mwfn'.format(il),
        'md-{:.3f}_q1.mwfn'.format(il),
        'md-{:.3f}_q-1.mwfn'.format(il),
        '2',
        '5',
        '6',
        '7',
        '8',
        '0',
        '0',
        'q',
    ]
    input_name = 'Fukui_multiwfn.inp'
    with open(input_name, 'w') as f:
        f.write('\n'.join(multiwfn_input))

    print('Calculate Fukui functions ...')
    os.system('Multiwfn {} < {} > /dev/null'.format('md-{:.3f}_q0.mwfn'.format(il), input_name))

    os.chdir('../')
