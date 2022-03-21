import os, sys
from setuptools import setup, find_packages
from setuptools.command.build_ext import build_ext
from pynep import __version__

DIR = os.path.abspath(os.path.dirname(__file__))
sys.path.append(os.path.join(DIR, "pybind11"))
from pybind11.setup_helpers import Pybind11Extension
del sys.path[-1]

# Avoid a gcc warning below:
# cc1plus: warning: command line option ‘-Wstrict-prototypes’ is valid
# for C/ObjC but not for C++
class BuildExt(build_ext):
    def build_extensions(self):
        if '-Wstrict-prototypes' in self.compiler.compiler_so:
            self.compiler.compiler_so.remove('-Wstrict-prototypes')
        if '-Wsign-compare' in self.compiler.compiler_so:
            self.compiler.compiler_so.remove('-Wsign-compare')
            self.compiler.compiler_so.append('-Wno-sign-compare')
        super().build_extensions()


#generatenew
module_Nep = Pybind11Extension('pynep.nep',
        sources = ['nep_cpu/src/pynep.cpp', 'nep_cpu/src/nep.cpp'],
        extra_compile_args=['-std=c++11'],
        )

with open('README.md') as f:
    long_description = f.read()

setup(
    name="pynep",
    version=__version__,
    author="Fan ZheYong, Wang Junjie",
    author_email="141120108@smail.nju.edu",
    url="https://git.nju.edu.cn/gaaooh/magus",
    packages=find_packages(),
    python_requires=">=3.8",
    install_requires=[
        "numpy<1.22.0",     # numba not support numpy >= 1.22.0
        "ase>=3.18",
    ],
    license="MIT",
    description="PyNEP: Python tools for NEP",
    long_description=long_description,
    long_description_content_type="text/markdown",
    ext_modules=[module_Nep], 
    cmdclass={'build_ext': BuildExt},
    # entry_points={"console_scripts": [""]},
)
