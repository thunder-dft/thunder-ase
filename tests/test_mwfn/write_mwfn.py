from ase.build import molecule
import ase
from ase.visualize import view

from thunder_ase.fireball import Fireball

#atoms = molecule('C2H2')
#del(atoms[2:])
atoms = ase.Atoms('Fe2', positions=[[0, 0, 0], [0, 0, 2.5]])
Fdata_path = '/home/ren/Programs/Github/thunder-ase/data/Fdata-McWEDA-0.15-3SN.Hs3.75.Cs4.00p4.45.Os3.35p3.80-3SNP.Fes5.30p5.30d4.80'

#atoms = ase.Atoms('Si')
#Fdata_path = '/home/ren/Programs/Github/thunder-ase/data/Fdata-McWEDA-0.15-3SN.Sis4.8p5.35'

kwargs = {'iwriteout_charges': 1,  # Writing out the charges.
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
