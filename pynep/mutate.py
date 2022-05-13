import abc
import numpy as np
from ase.data import covalent_radii, atomic_numbers
from ase.geometry import get_distances


class Mutation(abc.ABC):
    def __init__(self, d_ratio=0.8, min_dis_mat=None, attempt_number=200):
        self.d_ratio = d_ratio
        self.min_dis_mat = min_dis_mat
        self.attempt_number = attempt_number
        self.descriptor = self.__class__.__name__

    def get_tolerance_distances(self, atoms):
        if self.min_dis_mat is None:
            unique_symbol = set(atoms.get_chemical_symbols())
            radius = {symbol: self.d_ratio * covalent_radii[atomic_numbers[symbol]] for symbol in unique_symbol}
            min_dis_mat = {(symbol1, symbol2): radius[symbol1] + radius[symbol2] for symbol1 in unique_symbol
                                                                                 for symbol2 in unique_symbol}
        else:
            min_dis_mat = self.min_dis_mat
        tolerance_distances = [[min_dis_mat[(symbol1, symbol2)] for symbol1 in atoms.get_chemical_symbols()] 
                                                                for symbol2 in atoms.get_chemical_symbols()]
        return np.array(tolerance_distances)

    def generate(self, atoms, number=10):
        new_frames = []
        for _ in range(number):
            new_atoms = self.mutate(atoms)
            if new_atoms is not None:
                new_frames.append(new_atoms)
        return new_frames

    @abc.abstractmethod
    def mutate(self, atoms):
        pass

    @staticmethod
    def check_distance(atoms, indices, tolerance_distances):
        distances = get_distances(atoms.positions, atoms.positions[indices], atoms.cell, atoms.pbc)[1]
        for i, j in enumerate(indices):
            distances[j, i] = 100
        if np.any(distances < tolerance_distances[:, indices]):
            return False
        return True


class Combine(Mutation):
    def __init__(self, *op_list):
        self.op_list = op_list
        self.descriptor = self.__class__.__name__
    
    def mutate(self, atoms):
        a = atoms.copy()
        for op in self.op_list:
            a = op.mutate(a)
            if a is None:
                return None
        return a


class Rattle(Mutation):
    def __init__(self, probability=0.8, rattle_range=2, d_ratio=0.8, min_dis_mat=None, attempt_number=200):
        super().__init__(d_ratio=d_ratio, min_dis_mat=min_dis_mat, attempt_number=attempt_number)
        self.probability = probability
        self.rattle_range = rattle_range

    @staticmethod
    def pos_add_sphere(rattle_strength):
        r = rattle_strength * np.random.rand()**(1/3)
        theta = np.random.uniform(low=0, high=2*np.pi)
        phi = np.random.uniform(low=0, high=np.pi)
        pos_add = r * np.array([np.cos(theta)*np.sin(phi),
                                np.sin(theta)*np.sin(phi),
                                np.cos(phi)])
        return pos_add

    def mutate(self, atoms):
        a = atoms.copy()
        tolerance_distances = self.get_tolerance_distances(a)
        Natoms = len(a)
        indices_to_rattle = np.arange(0, Natoms)[np.random.rand(Natoms) < self.probability]
        indices_to_rattle = np.random.permutation(indices_to_rattle)
        if len(indices_to_rattle) == 0:
            return None
        success_rattle = 0
        for i in indices_to_rattle:
            posi_0 = np.copy(a.positions[i])
            for j in range(self.attempt_number):
                rattle_range = self.rattle_range - 0.8 * self.rattle_range / self.attempt_number * j
                pos_add = self.pos_add_sphere(rattle_range)
                a.positions[i] += pos_add
                if self.check_distance(a, [i], tolerance_distances):
                    success_rattle += 1
                    break 
                else:
                    a.positions[i] = posi_0
        if success_rattle > 0.5 * len(indices_to_rattle):
            return a
        else:
            return None


class Strain(Mutation):
    def __init__(self, cell_cut=1., sigma=0.1, d_ratio=0.8, min_dis_mat=None, attempt_number=200):
        super().__init__(d_ratio=d_ratio, min_dis_mat=min_dis_mat, attempt_number=attempt_number)
        self.cell_cut = cell_cut
        self.sigma = sigma

    def mutate(self, atoms):
        a = atoms.copy()
        tolerance_distances = self.get_tolerance_distances(a)
        for i in range(self.attempt_number):
            cell_cut = self.cell_cut - 0.8 * self.cell_cut / self.attempt_number * i
            lat_gauss = np.clip(np.random.normal(0, self.sigma, 6), -self.sigma, self.sigma) * cell_cut
            strain = np.array([
                [1+lat_gauss[0], lat_gauss[1]/2, lat_gauss[2]/2],
                [lat_gauss[1]/2, 1+lat_gauss[3], lat_gauss[4]/2],
                [lat_gauss[2]/2, lat_gauss[4]/2, 1+lat_gauss[5]]
                ])
            new_cell = atoms.cell[:] @ strain
            a.set_cell(new_cell, scale_atoms=True)
            if self.check_distance(a, np.arange(len(a)), tolerance_distances):
                return a
        return None


class Swap(Mutation):
    def __init__(self, swap_probability=0.8, d_ratio=0.8, min_dis_mat=None, attempt_number=200):
        super().__init__(d_ratio=d_ratio, min_dis_mat=min_dis_mat, attempt_number=attempt_number)
        self.swap_probability = swap_probability

    def mutate(self, atoms):
        a = atoms.copy()
        num_swaps = np.random.randint(1, max(int(self.swap_probability * len(atoms)), 2))
        unique_symbols = list(set(atoms.get_chemical_symbols()))
        tolerance_distances = self.get_tolerance_distances(a)
        if len(unique_symbols) < 2:
            return None
        success_change = 0
        for _ in range(num_swaps):
            for n in range(self.attempt_number):
                s1, s2 = np.random.choice(unique_symbols, 2, replace=False)
                s1_list = [i for i in range(len(a)) if a[i].symbol == s1]
                s2_list = [i for i in range(len(a)) if a[i].symbol == s2]
                i = np.random.choice(s1_list)
                j = np.random.choice(s2_list)
                a[i].symbol, a[j].symbol = a[j].symbol, a[i].symbol
                if self.check_distance(a, [i, j], tolerance_distances):
                    success_change += 1
                    break 
                else:
                    a[i].symbol, a[j].symbol = a[j].symbol, a[i].symbol
        if success_change < 0.6 * num_swaps:
            return None
        return a


class ChangeAtomType(Mutation):
    def __init__(self, change_probability=0.8, d_ratio=0.8, min_dis_mat=None, attempt_number=200):
        super().__init__(d_ratio=d_ratio, min_dis_mat=min_dis_mat, attempt_number=attempt_number)
        self.change_probability = change_probability

    def mutate(self, atoms):
        a = atoms.copy()
        unique_symbols = list(set(atoms.get_chemical_symbols()))
        tolerance_distances = self.get_tolerance_distances(a)
        if len(unique_symbols) < 2:
            return None
        success_change = 0
        for _ in range(self.attempt_number):
            for i in range(len(a)):
                if np.random.rand() < self.change_probability:
                    s1 = a[i].symbol
                    s2 = np.random.choice([s for s in unique_symbols if s != s1])
                    a[i].symbol = s2
                if self.check_distance(a, [i], tolerance_distances):
                    success_change += 1
                    break 
                else:
                    a[i].symbol = s1
        return a
