from ase.build import bulk

from thunder_ase.fireball import Fireball

atoms = bulk('Si', 'diamond', a=5.459)
Fdata_path = '/home/ren/Programs/Github/thunder-ase/data/Fdata-McWEDA-0.15-3SN.Sis4.8p5.35'
kwargs = {'kpt_size': [3, 3, 3],
          'iwriteout_charges': 1,  # Writing out the charges.
          'iwriteout_cdcoeffs': 1,  # Writing out orbital info.
          'taurelax': 5.0,
          'efermi_T': 200.0,
          'ifix_CHARGES': 0,
          'max_scf_iterations_set': 100,
          'scf_tolerance_set': 0.00000001,
          'beta_set': 0.04,
          }

calc = Fireball(command='/home/ren/bin/lightning.3.x',
                Fdata_path=Fdata_path,
                **kwargs)
atoms.set_calculator(calc)

e0 = atoms.get_potential_energy()
efermi = atoms.calc.get_fermi_level()

print("The energy is {:.3f} eV.".format(e0))
print("The Fermi Level is {:.3f} eV.".format(efermi))
calc.write_mwfn()
print("Done!")
