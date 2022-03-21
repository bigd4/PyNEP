from pynep.calculate import NEP
from ase.io import read
import numpy as np

a = read('POSCAR')
calc = NEP({'Te': 0, 'Pb': 1}, "nep.txt")
a.set_calculator(calc)
e = a.get_potential_energy()
print(e)
f_py = a.get_forces()
f_c = np.loadtxt("force_analytical.out")
np.testing.assert_almost_equal(f_py, f_c, decimal=5)
