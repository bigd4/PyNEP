import numpy as np
from ase.calculators.calculator import Calculator, all_changes, PropertyNotImplementedError
from .nep import NepCalculator


class NEP(Calculator):

    implemented_properties = [
        "energy", 
        "forces", 
        "stress", 
        "descriptor",
        ]

    def __init__(self, type_dict, model_file="nep.txt", **kwargs) -> None:
        Calculator.__init__(self, **kwargs)
        self.calc = NepCalculator(model_file)
        self.type_dict = type_dict

    def __repr__(self):
        info = self.calc.info
        ret = "NEP {} calculator with {} symbols: ".format(info['version'], len(self.type_dict))
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
        self.results["forces"] = np.array(self.calc.getForces()).reshape(3, -1).transpose(1, 0)

        if "stress" in properties:
            virial = np.sum(np.array(self.calc.getVirials()).reshape(9, -1), axis=1)
            if sum(atoms.get_pbc()) > 0:
                stress = -0.5 * (virial.copy() + virial.copy().T) / atoms.get_volume()
                self.results['stress'] = stress.flat[[0, 4, 8, 5, 2, 1]]
            else:
                raise PropertyNotImplementedError

        if "descriptor" in properties:
            self.results['descriptor'] = np.array(self.calc.getDescriptors()).reshape(-1, len(atoms)).transpose(1, 0)
