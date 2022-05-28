import abc
import numpy as np
from ase.data import covalent_radii, atomic_numbers
from ase.geometry import get_distances


class Mutation(abc.ABC):
    """Base class of Muatation
    """
    def __init__(self, d_ratio=0.8, min_dis_mat=None, attempt_number=200):
        """Initialize Mutation Object

        Args:
            d_ratio (float, optional): Use d_ratio * covalent_radii to calculate min_dis_mat.
                If min_dis_mat is given, this will be ignored.
                Defaults to 0.8. 
            min_dis_mat (dict , optional): min distance between two atoms.
                Such as {('C', 'C'): 1.} 
                Defaults to None. 
            attempt_number (int, optional): Max times to try to mutate. 
                Defaults to 200.
        """
        self.d_ratio = d_ratio
        self.min_dis_mat = min_dis_mat
        self.attempt_number = attempt_number
        self.descriptor = self.__class__.__name__

    def get_tolerance_distances(self, atoms):
        """get min distance between atoms

        Args:
            atoms (Atoms): atoms to constuct the tolerance distances

        Returns:
            2d list: min distances matrix
        """
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
        """generate atoms

        Args:
            atoms (Atoms): atoms to be mutated
            number (int, optional): number of new atoms to be generated. Defaults to 10.

        Returns:
            A list of atoms obeject: Generated atoms
        """
        new_frames = []
        for _ in range(number):
            new_atoms = self.mutate(atoms)
            if new_atoms is not None:
                new_frames.append(new_atoms)
        return new_frames

    @abc.abstractmethod
    def mutate(self, atoms):
        """mutate the atoms, every subclass must have this method

        Args:
            atoms (Atoms): atoms to be mutated
        """
        pass

    @staticmethod
    def check_distance(atoms, indices, tolerance_distances):
        """check distances of atoms

        Args:
            atoms (Atoms): atoms to be check
            indices (1d list): indices of atoms to be cheak
            tolerance_distances (allow distance between atoms): 2d list

        Returns:
            bool: If the atoms meet the requirements
        """
        distances = get_distances(atoms.positions, atoms.positions[indices], atoms.cell, atoms.pbc)[1]
        for i, j in enumerate(indices):
            distances[j, i] = 100
        if np.any(distances < tolerance_distances[:, indices]):
            return False
        return True


class Combine(Mutation):
    """Class to combine a series of mutations
    """
    def __init__(self, *op_list):
        """Combile mutations

        Example:
        >>> m1 = Rattle()
        >>> m2 = Strain()
        >>> new = Combine(m1, m2)
        """
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
    """Randomly rattle atoms
    """
    def __init__(self, probability=0.8, rattle_range=2, d_ratio=0.8, min_dis_mat=None, attempt_number=200):
        """

        Args:
            probability (float, optional): rattle probability of each atom
            rattle_range (int, optional): atoms will only rattle in this range. Defaults to 2.
        """
        super().__init__(d_ratio=d_ratio, min_dis_mat=min_dis_mat, attempt_number=attempt_number)
        self.probability = probability
        self.rattle_range = rattle_range

    @staticmethod
    def pos_add_sphere(rattle_strength):
        """randomly new positions in a sphere

        Args:
            rattle_strength (float): radius of the sphere

        Returns:
            array: positions to add
        """
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
    """Randomly change the cell
    """
    def __init__(self, cell_cut=1., sigma=0.1, d_ratio=0.8, min_dis_mat=None, attempt_number=200):
        """

        Args:
            cell_cut (_type_, optional): max value of strain matrix. Defaults to 1..
            sigma (float, optional): sigma of gaussian. Defaults to 0.1.
        """
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
    """Swap atoms
    """
    def __init__(self, swap_probability=0.8, d_ratio=0.8, min_dis_mat=None, attempt_number=200):
        """

        Args:
            swap_probability (float, optional): Max swap number is swap_probability * N_atoms. 
                Defaults to 0.8.
        """
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
    """Change type of atom
    """
    def __init__(self, change_probability=0.8, d_ratio=0.8, min_dis_mat=None, attempt_number=200):
        """

        Args:
            change_probability (float, optional): change probability of each atom. 
                Defaults to 0.8.
        """
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
