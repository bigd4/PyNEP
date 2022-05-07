import numpy as np
import spglib
import ase
from ase import Atoms
from ase.units import kJ, mol
from ase.cell import Cell
from phonopy import Phonopy
from phonopy.structure.atoms import PhonopyAtoms
from phonopy.phonon.band_structure import get_band_qpoints_and_path_connections
from phonopy.interface.phonopy_yaml import PhonopyYaml


def plot_band_structure(band_dict):
    import matplotlib.pyplot as plt
    from mpl_toolkits.axes_grid1 import ImageGrid    
    labels_path = band_dict['labels_path']
    frequencies = band_dict['frequencies']
    distances = band_dict['distances']
    fig = plt.figure()
    axs = ImageGrid(fig, 111, nrows_ncols=(1, len(labels_path)), axes_pad=0.2, label_mode="L")

    max_freq = max([np.max(fq) for fq in frequencies])
    max_dist = distances[-1][-1]
    xscale = max_freq / max_dist * 1.5
    distances_scaled = [d * xscale for d in distances]

    n = 0
    axs[0].set_ylabel("Frequency", fontsize=14)
    for i, path in enumerate(labels_path):
        axs[i].spines['bottom'].set_linewidth(1.5)
        axs[i].spines['left'].set_linewidth(1.5)
        axs[i].spines['right'].set_linewidth(1.5)
        axs[i].spines['top'].set_linewidth(1.5)
        axs[i].tick_params(labelsize=14)
        xticks = [distances_scaled[n][0]]
        for label in path[:-1]:
            xticks.append(distances_scaled[n][-1])
            axs[i].plot([distances_scaled[n][-1], distances_scaled[n][-1]], 
                        [0, max_freq],
                        linewidth=2,
                        linestyle=":",
                        c='grey')
            axs[i].plot(distances_scaled[n], 
                        frequencies[n], 
                        linewidth=2,
                        c='g')
            n += 1
        axs[i].set_xlim(xticks[0], xticks[-1])
        axs[i].set_xticks(xticks)
        axs[i].set_xticklabels(path)
        axs[i].plot([xticks[0], xticks[-1]], 
                    [0, 0], 
                    linewidth=1,
                    c='black')
    plt.savefig('phono.png')
    return axs

def ase2phono(atoms):
    return PhonopyAtoms(symbols=atoms.get_chemical_symbols(),
                        cell=atoms.cell.array,
                        scaled_positions=atoms.get_scaled_positions())

def phono2ase(cell):
    return Atoms(symbols=cell.get_chemical_symbols(),
                 scaled_positions=cell.get_scaled_positions(),
                 cell=cell.get_cell(),
                 pbc=True)

class PhonoCalc:
    def __init__(self, calc, dim='Auto', mesh=(10, 10, 10), t_step=10, t_max=0., t_min=0.):
        self.calc = calc
        self.dim = dim
        self.mesh = mesh
        self.t_step = t_step
        self.t_max = t_max
        self.t_min = t_min

    def calculate(self, atoms, properties=['band', 'dos', 'pdos', 'thermal', 'ZPE']):
        results = {}
        try:
            print('Calculating force constants...')
            results['force_constants'] = self.get_force_constants(atoms)
            phpy_yaml = PhonopyYaml(settings=
                    {'force_sets': True,
                     'displacements': True,
                     'force_constants': True,
                     'born_effective_charge': True,
                     'dielectric_constant': True})
            phpy_yaml.set_phonon_info(self.phonon)
            atoms.info['phono_info'] = phpy_yaml._yaml
        except:
            print('Fail to collect force constants')
            return atoms
        if 'band' in properties:
            print('Calculating band structure...')
            results['band_dict'] = self.get_band_structure(atoms)
            plot_band_structure(results['band_dict'])
        if 'dos' in properties:
            print('Calculating total dos...')
            results['dos_dict'] = self.get_dos()
        if 'pdos' in properties:
            print('Calculating partial dos...')
            results['pdos_dict'] = self.get_pdos()
        if 'thermal' in properties:
            print('Calculating thermal properties...')
            results['tp_dict'] = self.get_thermal(self.t_min, self.t_step, self.t_max)
        if 'ZPE' in properties:
            print('Calculating ZPE...')
            results['ZPE'] = self.get_thermal(0., 1., 1.)['free_energy'][0] * kJ / mol
        atoms.info.update(results)
        return atoms

    def get_supercell_matrix(self, atoms):
        if isinstance(self.dim, str):
            if self.dim == 'Auto':
                lengths = atoms.cell.lengths()
                return np.round(10 / lengths).astype('int')
        else:
            return self.dim

    def get_force_constants(self, atoms):
        unitcell = ase2phono(atoms)
        supercell_matrix = self.get_supercell_matrix(atoms)
        primitive_lattice = spglib.find_primitive(atoms)[0]
        primitive_matrix = primitive_lattice @ np.linalg.inv(atoms.cell)
        self.phonon = Phonopy(
            unitcell=unitcell, 
            supercell_matrix=supercell_matrix,
            primitive_matrix=primitive_matrix)
        self.phonon.generate_displacements(distance=0.01)
        supercells = self.phonon.get_supercells_with_displacements()
        set_of_forces = []
        for cell in supercells:
            forces = self.calc.get_forces(phono2ase(cell))
            forces -= np.mean(forces, axis=0)
            set_of_forces.append(forces)
        set_of_forces = np.array(set_of_forces)
        self.phonon.produce_force_constants(forces=set_of_forces)
        return self.phonon.force_constants

    def get_dos(self):
        self.phonon.run_mesh(self.mesh)
        self.phonon.run_total_dos(use_tetrahedron_method=True)
        return self.phonon.get_total_dos_dict()

    def get_pdos(self):
        self.phonon.run_mesh(self.mesh, with_eigenvectors=True, is_mesh_symmetry=False)
        self.phonon.run_projected_dos()
        return self.phonon.get_projected_dos_dict()

    def get_thermal(self, t_min, t_step, t_max):
        self.phonon.run_mesh(self.mesh)
        self.phonon.run_thermal_properties(t_min=t_min,
                                           t_step=t_step, 
                                           t_max=t_max)
        return self.phonon.get_thermal_properties_dict()

    def get_band_structure(self, atoms):
        cell = Cell(spglib.find_primitive(atoms)[0])
        special_points = cell.get_bravais_lattice().get_special_points()
        labels_path = ase.dft.kpoints.parse_path_string(cell.bandpath().path)
        labels, path = [], []
        for label_path in labels_path:
            p = []
            for l in label_path:
                labels.append(l)
                p.append(special_points[l].tolist())
            path.append(p)
        qpoints, connections = get_band_qpoints_and_path_connections(path, npoints=51)
        self.phonon.run_band_structure(qpoints, path_connections=connections, labels=labels)
        bands_dict = self.phonon.get_band_structure_dict()
        bands_dict['labels_path'] = labels_path
        return bands_dict
 
