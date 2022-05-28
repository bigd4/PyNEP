import numpy as np
from ase.calculators.singlepoint import SinglePointCalculator
from ase import Atoms


def load_nep(filename):
    """Read Atoms objects from file

    Args:
        filename (str): Name of the file to read from

    Returns:
        A list of Atoms objects: frames
    """
    frames, n_atoms, has_virial = [], [], []
    with open(filename, 'r') as f:
        line = f.readline()
        n_frames = int(line)
        for _ in range(n_frames):
            line = f.readline()
            n_atoms.append(int(line.split()[0]))
            has_virial.append(int(line.split()[1]))
        for i in range(n_frames):
            d = {}
            line = f.readline()
            info = list(map(float, line.split()))
            d['energy'] = info[0]
            if has_virial[i] == 1:
                xx, yy, zz, xy, yz, zx = info[1:]
                virial = np.array([[xx, xy, zx], [xy, yy, yz], [zx, yz, zz]])
            line = f.readline()
            cell = np.array(list(map(float, line.split()))).reshape(3, 3)
            if has_virial[i] == 1:
                d['stress'] = -virial / np.linalg.det(cell)
            symbols, positions, d['forces'] = [], [], []
            for _ in range(n_atoms[i]):
                line = f.readline()
                symbols.append(line.split()[0])
                positions.append(list(map(float, line.split()[1: 4])))
                d['forces'].append(list(map(float, line.split()[4: 7])))
            atoms = Atoms(cell=cell, positions=positions, symbols=symbols, pbc=True)
            calc = SinglePointCalculator(atoms, **d)
            atoms.info.update(d)
            atoms.set_calculator(calc)
            frames.append(atoms)
    return frames 


def dump_nep(filename, frames):
    """dump data in nep format

    Args:
        filename (str): Name of the file to read from
        frames (list): a list of Atoms objects.
    """
    for atoms in frames:
        for p in ['energy', 'forces', 'stress']:
            if p not in atoms.info:
                try:
                    atoms.info[p] = atoms.calc.get_property(p)
                except:
                    pass
    with open(filename, 'w') as f:
        n_frames = len(frames)
        f.write(str(n_frames) + '\n')
        for atoms in frames:
            has_virial = int('stress' in atoms.info)
            f.write('{} {} \n'.format(len(atoms), has_virial))
        for atoms in frames:
            ret = str(atoms.info['energy'])
            if 'stress' in atoms.info:
                if len(atoms.info['stress']) == 6:
                    virial = -atoms.info['stress'][[0, 1, 2, 5, 3, 4]] * atoms.get_volume()
                else:
                    virial = -atoms.info['stress'][[0, 1, 2, 0, 1, 2], [0, 1, 2, 1, 2, 0]] * atoms.get_volume()
                for v in virial:
                    ret += ' ' + str(v)
            ret += '\n{:.8e} {:.8e} {:.8e} {:.8e} {:.8e} {:.8e} {:.8e} {:.8e} {:.8e}\n'.format(*atoms.get_cell().reshape(-1))
            s = atoms.get_chemical_symbols()
            p = atoms.get_positions()
            forces = atoms.info['forces']
            for i in range(len(atoms)):
                ret += '{:2} {:>15.8e} {:>15.8e} {:>15.8e} {:>15.8e} {:>15.8e} {:>15.8e}\n'.format(s[i], *p[i], *forces[i])
            f.write(ret)
