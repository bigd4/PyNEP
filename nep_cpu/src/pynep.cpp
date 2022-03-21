/*
    Copyright 2017 Zheyong Fan, Ville Vierimaa, Mikko Ervasti, and Ari Harju
    This file is part of GPUMD.
    GPUMD is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    GPUMD is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
    You should have received a copy of the GNU General Public License
    along with GPUMD.  If not, see <http://www.gnu.org/licenses/>.
*/

/*----------------------------------------------------------------------------80
Usage:
    Compile:
        g++ -O3 main.cpp nep.cpp
    run:
        ./a.out
------------------------------------------------------------------------------*/

#include "nep.h"
#include "utility.h"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <time.h>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;

const int MN = 1000;

struct Atom {
  int N;
  std::vector<int> num_cells, type, NN_radial, NL_radial, NN_angular, NL_angular;
  std::vector<double> box, ebox, position, r12, potential, force, virial;
};

class NepCalculator
{
  public:
    NepCalculator(std::string);
    void setAtoms(int, std::vector<int>, std::vector<double>, std::vector<double>);
    void calculate();
    std::vector<double> getPotentialEnergy();
    std::vector<double> getForces();
    std::vector<double> getVirials();
  private:
    Atom atom;
    NEP3 calc;
    std::string model_file;
};

NepCalculator::NepCalculator(std::string _model_file)
{
  model_file = _model_file;
}

void NepCalculator::setAtoms(
  int _N,
  std::vector<int> _type,
  std::vector<double> _box,
  std::vector<double> _position)
{
  Atom _atom;
  _atom.N = _N;
  _atom.num_cells.resize(3);
  _atom.box = _box;
  _atom.box.resize(18);
  _atom.ebox.resize(18);
  get_inverse(_atom.box.data());

  _atom.type = _type;
  _atom.position = _position;

  _atom.NN_radial.resize(_atom.N);
  _atom.NL_radial.resize(_atom.N * MN);
  _atom.NN_angular.resize(_atom.N);
  _atom.NL_angular.resize(_atom.N * MN);
  _atom.r12.resize(_atom.N * MN * 6);
  _atom.position.resize(_atom.N * 3);
  _atom.potential.resize(_atom.N);
  _atom.force.resize(_atom.N * 3);
  _atom.virial.resize(_atom.N * 9);

  atom = _atom;
  calc = NEP3(model_file);
}

void NepCalculator::calculate()
{
  find_neighbor_list_small_box(
    calc.paramb.rc_radial, calc.paramb.rc_angular, atom.N, atom.box, atom.position,
    atom.num_cells, atom.ebox, atom.NN_radial, atom.NL_radial, atom.NN_angular, atom.NL_angular,
    atom.r12);

  calc.compute(
    atom.NN_radial, atom.NL_radial, atom.NN_angular, atom.NL_angular, atom.type, atom.r12,
    atom.potential, atom.force, atom.virial);
}

std::vector<double> NepCalculator::getPotentialEnergy()
{
  return atom.potential;
}

std::vector<double> NepCalculator::getForces()
{
  return atom.force;
}

std::vector<double> NepCalculator::getVirials()
{
  return atom.virial;
}

PYBIND11_MODULE(nep, m){
    m.doc() = "nep";
    py::class_<NepCalculator>(m, "NepCalculator")
		.def(py::init<std::string>())
		.def("setAtoms", &NepCalculator::setAtoms)
		.def("calculate", &NepCalculator::calculate)
		.def("getPotentialEnergy", &NepCalculator::getPotentialEnergy)
    .def("getForces", &NepCalculator::getForces)
    .def("getVirials", &NepCalculator::getVirials)
		;
}
