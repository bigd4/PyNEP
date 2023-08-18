import numpy as np
from ase.calculators.calculator import Calculator, all_changes, PropertyNotImplementedError
from .nep import NepCalculator


class NEP(Calculator):

    """ASE calculator for NEP (see https://github.com/brucefan1983/GPUMD)
    supported properties: energy, forces, stress, descriptor, latent descriptor

    Examples:

    Use C_2022_NEP3.txt to calculate properties of diamond 
    (from https://github.com/brucefan1983/GPUMD/tree/master/potentials/nepCalculate)

    energy, forces and stress:

    >>> calc = NEP('C_2022_NEP3.txt')
    >>> atoms = bulk('C', 'diamond', cubic=True)
    >>> atoms.set_calculator(calc)
    >>> energy = atoms.get_potential_energy()
    >>> forces = atoms.get_forces()
    >>> stress = atoms.get_stress()

    descriptor and latent descriptor:

    >>> des = calc.get_property('descriptor', atoms)
    >>> lat = calc.get_property('latent', atoms)
    """

    implemented_properties = [
        "energy",
        "energies",
        "forces",
        "stress",
        "descriptor",
        "latent",
    ]

    def __init__(self, model_file="nep.txt", **kwargs) -> None:
        """Initialize calculator

        Args:
            model_file (str, optional): filename of nep model. Defaults to "nep.txt".
        """
        Calculator.__init__(self, **kwargs)
        self.calc = NepCalculator(model_file)
        self.type_dict = {e: i for i, e in enumerate(
            self.calc.info["element_list"])}

    def __repr__(self):
        info = self.calc.info
        ret = "NEP {} calculator with {} symbols: ".format(
            info['version'], len(self.type_dict))
        info.pop("version")
        for key in self.type_dict:
            ret += key + " "
        ret += "\n" + "-" * 30
        for key, value in info.items():
            ret += "\n  {}: {}".format(key.ljust(20, ' '), value)
        ret += "\n" + "-" * 30 + "\n"
        return ret

    def calculate(
        self,
        atoms=None,
        properties=None,
        system_changes=all_changes,
    ):
        if properties is None:
            properties = self.implemented_properties

        Calculator.calculate(self, atoms, properties, system_changes)
        symbols = self.atoms.get_chemical_symbols()
        _type = [self.type_dict[k] for k in symbols]
        _box = atoms.cell.transpose(1, 0).reshape(-1).tolist()
        _position = atoms.get_positions().transpose(1, 0).reshape(-1).tolist()
        self.calc.setAtoms(len(atoms), _type, _box, _position)
        self.calc.calculate()
        self.results["energy"] = np.sum(self.calc.getPotentialEnergy())
        self.results["energies"] = np.array(self.calc.getPotentialEnergy())
        self.results["forces"] = np.array(
            self.calc.getForces()).reshape(3, -1).transpose(1, 0)

        if "stress" in properties:
            virial = np.sum(
                np.array(self.calc.getVirials()).reshape(9, -1), axis=1)
            if sum(atoms.get_pbc()) > 0:
                stress = -0.5 * (virial.copy() +
                                 virial.copy().T) / atoms.get_volume()
                self.results['stress'] = stress.flat[[0, 4, 8, 5, 2, 1]]
            else:
                raise PropertyNotImplementedError

        if "descriptor" in properties:
            self.results['descriptor'] = np.array(
                self.calc.getDescriptors()).reshape(-1, len(atoms)).transpose(1, 0)

        if "latent" in properties:
            self.results['latent'] = np.array(
                self.calc.getLatent()).reshape(-1, len(atoms)).transpose(1, 0)


class JointCalculator(Calculator):
    
    implemented_properties = [
        "energy",
        "forces",
        "stress",
    ]

    def __init__(self, *args, **kwargs) -> None:
        Calculator.__init__(self, **kwargs)
        self.calc_list = args

    def calculate(
        self,
        atoms=None,
        properties=None,
        system_changes=all_changes,
    ):
        if properties is None:
            properties = self.implemented_properties

        Calculator.calculate(self, atoms, properties, system_changes)

        self.results["energy"] = 0.
        self.results["forces"] = np.zeros((len(atoms), 3))
        self.results['stress'] = np.zeros(6)
        for calc in self.calc_list:
            self.results["energy"] += calc.get_potential_energy(self.atoms)
            self.results["forces"] += calc.get_forces(self.atoms)
            if "stress" in properties:
                self.results['stress'] += calc.get_stress(self.atoms)
