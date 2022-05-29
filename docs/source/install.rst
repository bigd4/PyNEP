Installation
============

Requirements
------------

========================== =========
Package                                                 version
========================== =========
`Python`_                  >= 3.8
`NumPy`_                   < 1.22.0
`SciPy`_                   >= 1.1
`ase`_                     >= 3.18.0
========================== =========

.. _Python: https://www.python.org/
.. _Numpy: https://docs.scipy.org/doc/numpy/reference/
.. _SciPy: https://docs.scipy.org/doc/scipy/reference/
.. _ase: https://wiki.fysik.dtu.dk/ase/index.html


By pip
------

.. code:: shell

   $ pip install git+https://github.com/bigd4/PyNEP.git

By setup.py
-----------

.. code:: shell

   $ git clone --recursive https://github.com/bigd4/PyNEP.git
   $ cd pynep
   $ python setup.py install

From Source
-----------

.. code:: shell

   $ git clone --recursive https://github.com/bigd4/PyNEP.git
   $ pip install -r requirements.txt
   $ cd pynep/nep_cpu
   $ mkdir build
   $ cd build
   $ cmake .. && make
   $ cp nep.so ../../pynep

Add ``PyNEP`` to your ``PYTHONPATH`` environment variable in your ``~/.bashrc`` file.

.. code:: shell

   $ export PYTHONPATH=<path-to-pynep-package>:$PYTHONPATH
