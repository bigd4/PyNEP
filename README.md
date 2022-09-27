# PyNEP 
[![Documentation Status](https://readthedocs.org/projects/pynep/badge/?version=latest)](https://pynep.readthedocs.io/en/latest/)

**[PyNEP](https://pynep.readthedocs.io/en/latest/)** is a python interface of the machine learning potential **NEP** used in **[GPUMD](https://github.com/brucefan1983/GPUMD)**.
## Features

- ase calculator of NEP
- descriptor and latent descriptor calculation of atoms
- load and dump for GPUMD dataset
- phonopy calculation of NEP (need phonopy and spglib)
- structures select
    + Farthest Point Sampling
## Installation

### Requirements


|  Package  | version |
|  ----  | ----  |
| [Python](https://www.python.org/) | >=     3.8 |
| [NumPy](https://docs.scipy.org/doc/numpy/reference/) | <       1.22.0 |
|[SciPy](https://docs.scipy.org/doc/scipy/reference/)|>=     1.1|
|[ase](https://wiki.fysik.dtu.dk/ase/index.html)|>=     3.18.0|


### By pip 

```shell
$ pip install git+https://github.com/bigd4/PyNEP.git
```

### By setup.py 

```shell
$ git clone --recursive https://github.com/bigd4/PyNEP.git
$ cd PyNEP
$ python setup.py install
```

### From Source

```shell
$ git clone --recursive https://github.com/bigd4/PyNEP.git
$ cd PyNEP/nep_cpu
$ mkdir build
$ cd build
$ cmake .. && make
$ cp nep.so ../../PyNEP
```

Add `PyNEP` to your [`PYTHONPATH`](https://wiki.fysik.dtu.dk/ase/install.html#envvar-PYTHONPATH) environment variable in your `~/.bashrc` file.

```shell
$ export PYTHONPATH=<path-to-PyNEP-package>:$PYTHONPATH
```

## Usage

```python
from ase.build import bulk
atoms = bulk('C', 'diamond', cubic=True)

# calculate energy and forces
from pynep.calculate import NEP
calc = NEP('nep.txt')
atoms = bulk('C', 'diamond', cubic=True)
atoms.set_calculator(calc)
energy = atoms.get_potential_energy()
forces = atoms.get_forces()
stress = atoms.get_stress()  # stress in ase is different from virial in GPUMD

# calculate descriptors and latent descriptors
des = calc.get_property('descriptor', atoms)
lat = calc.get_property('latent', atoms)

# load and dump GPUMD data
from pynep.io import load_nep, dump_nep
dump_nep('C.in', [atoms])
atoms = load_nep('C.in')[0]

# calculate band strucuture, dos and thermal properties (need spglib and phonopy)
from pynep.phono import PhonoCalc
phono_calc = PhonoCalc(calc)
phono_calc.calculate(atoms)
```

