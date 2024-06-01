#include "nep.h"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <time.h>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/iostream.h>
namespace py = pybind11;

struct Atom {
  int N;
  std::vector<int> type;
  std::vector<double> box, position, potential, force, virial, descriptor, latent;
};

class NepCalculator
{
  public:
    NepCalculator(std::string);
    void setAtoms(int, std::vector<int>, std::vector<double>, std::vector<double>);
    void calculate();
    py::dict info;
    std::vector<double> getPotentialEnergy();
    std::vector<double> getForces();
    std::vector<double> getVirials();
    std::vector<double> getDescriptors();
    std::vector<double> getLatent();
  private:
    Atom atom;
    NEP3 calc;
    std::string model_file;
    bool HAS_CALCULATED=false;
};

NepCalculator::NepCalculator(std::string _model_file)
{
  model_file = _model_file;
  py::scoped_ostream_redirect stream(
    std::cout,                               // std::ostream&
    py::module_::import("sys").attr("stdout") // Python output
    );
  calc = NEP3(model_file);
  info["version"] = calc.paramb.version;
  info["zbl"] = calc.zbl.enabled;
  info["radial_cutoff"] = calc.paramb.rc_radial;
  info["angular_cutoff"] = calc.paramb.rc_angular;
  info["n_max_radial"] = calc.paramb.n_max_radial;
  info["n_max_angular"] = calc.paramb.n_max_angular;
  info["basis_size_radial"] = calc.paramb.basis_size_radial;
  info["basis_size_angular"] = calc.paramb.basis_size_angular;
  info["l_max_3body"] = calc.paramb.L_max;
  info["num_node"] = calc.annmb.dim;
  info["num_para"] = calc.annmb.num_para;
  info["element_list"] = calc.element_list;
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

  _atom.potential.resize(_atom.N);
  _atom.force.resize(_atom.N * 3);
  _atom.virial.resize(_atom.N * 9);
  _atom.descriptor.resize(_atom.N * calc.annmb.dim);
  _atom.latent.resize(_atom.N * calc.annmb.num_neurons1);
  atom = _atom;
  HAS_CALCULATED = false;
}

void NepCalculator::calculate()
{
  if (!HAS_CALCULATED){
    calc.compute(atom.type, atom.box, atom.position, atom.potential, atom.force, atom.virial);
    HAS_CALCULATED = true;
  }
}

std::vector<double> NepCalculator::getPotentialEnergy()
{
  calculate();
  return atom.potential;
}

std::vector<double> NepCalculator::getForces()
{
  calculate();
  return atom.force;
}

std::vector<double> NepCalculator::getVirials()
{
  calculate();
  return atom.virial;
}

std::vector<double> NepCalculator::getDescriptors()
{
  calc.find_descriptor(atom.type, atom.box, atom.position, atom.descriptor);
  return atom.descriptor;
}

std::vector<double> NepCalculator::getLatent()
{
  calculate();
  calc.find_latent_space(atom.type, atom.box, atom.position, atom.latent);
  return atom.latent;
}

PYBIND11_MODULE(nep, m){
    m.doc() = "nep";
    py::class_<NepCalculator>(m, "NepCalculator")
		.def(py::init<std::string>())
    .def_readonly("info", &NepCalculator::info)
		.def("setAtoms", &NepCalculator::setAtoms)
		.def("calculate", &NepCalculator::calculate)
		.def("getPotentialEnergy", &NepCalculator::getPotentialEnergy)
    .def("getForces", &NepCalculator::getForces)
    .def("getVirials", &NepCalculator::getVirials)
    .def("getDescriptors", &NepCalculator::getDescriptors)
    .def("getLatent", &NepCalculator::getLatent)
		;
}
