#include "nep.h"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <time.h>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;

struct Atom {
  int N;
  std::vector<int> type;
  std::vector<double> box, position, potential, force, virial;
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
  calc = NEP3(model_file);
}

void NepCalculator::setAtoms(
  int _N,
  std::vector<int> _type,
  std::vector<double> _box,
  std::vector<double> _position)
{
  Atom _atom;
  _atom.N = _N;
  _atom.box = _box;
  _atom.type = _type;
  _atom.position = _position;

  _atom.position.resize(_atom.N * 3);
  _atom.potential.resize(_atom.N);
  _atom.force.resize(_atom.N * 3);
  _atom.virial.resize(_atom.N * 9);

  atom = _atom;
}

void NepCalculator::calculate()
{
  calc.compute(atom.type, atom.box, atom.position, atom.potential, atom.force, atom.virial);
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
