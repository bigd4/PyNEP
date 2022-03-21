import numpy as np
from ase.calculators.calculator import Calculator, all_changes
from .nep import NepCalculator


class NEP(Calculator):

    implemented_properties = ["energy", "forces"]

    def __init__(self, type_dict, model_file="nep.txt", **kwargs) -> None:
        Calculator.__init__(self, **kwargs)
        self.calc = NepCalculator(model_file)
        self.type_dict = type_dict

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
        _box = atoms.cell.reshape(-1).tolist()
        _position = atoms.get_positions().transpose(1, 0).reshape(-1).tolist()
        self.calc.setAtoms(len(atoms), _type, _box, _position)
        self.calc.calculate()
        self.results['energy'] = np.sum(self.calc.getPotentialEnergy())
        self.results['forces'] = np.array(self.calc.getForces()).reshape(3, -1).transpose(1, 0)
