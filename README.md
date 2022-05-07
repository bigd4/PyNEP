# PyNEP - Python tools for NEP

## Features

- ase calculator of NEP
- descriptor and latent descriptor calculation of atoms
- load and dump for GPUMD dataset
- phonopy calculation of NEP (need phonopy)
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
$ cd pynep
$ python setup.py install
```

### From Source

```shell
$ git clone --recursive https://github.com/bigd4/PyNEP.git
$ cd pynep/nep_cpu
$ mkdir build
$ cd build
$ cmake .. && make
$ cp nep.so ../../pynep
```

Add `pynep` to your [`PYTHONPATH`](https://wiki.fysik.dtu.dk/ase/install.html#envvar-PYTHONPATH) environment variable in your `~/.bashrc` file.

```shell
$ export PYTHONPATH=<path-to-pynep-package>:$PYTHONPATH
```

## Usage

```python
from pynep.calculate import NEP
import numpy as np
from ase.build import bulk

a = bulk('TePb', 'rocksalt', 6.6)
calc = NEP({'Te': 0, 'Pb': 1}, "nep.txt")
a.set_calculator(calc)
e = a.get_potential_energy()
f = a.get_forces()
```

