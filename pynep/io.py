import numpy as np
from ase.calculators.singlepoint import SinglePointCalculator
from collections import OrderedDict
from ase import Atoms
import re

# copy from dpdata-0.2.6:  dpdata/xyz/quip_gap_xyz.py
# class QuipGapxyzSystems(object): handle_single_xyz_frame
# dpdata: https://github.com/deepmodeling/dpdata
def Proc_block(lines):
    atom_num = int(lines[0].strip('\n').strip())
    if len(lines) != atom_num + 2:
        raise RuntimeError("format error, atom_num=={}, {}!=atom_num+2".format(atom_num, len(lines)))
    data_format_line = lines[1].strip('\n').strip()+str(' ')
    field_value_pattern= re.compile(r'(?P<key>\S+)=(?P<quote>[\'\"]?)(?P<value>.*?)(?P=quote)\s+')
    prop_pattern = re.compile(r'(?P<key>\w+?):(?P<datatype>[a-zA-Z]):(?P<value>\d+)')

    data_format_list= [kv_dict.groupdict() for kv_dict in field_value_pattern.finditer(data_format_line)]
    field_dict = {}
    for item in data_format_list:
        field_dict[item['key'].lower()]=item['value']

    Properties = field_dict['properties']
    prop_list = [kv_dict.groupdict() for kv_dict in prop_pattern.finditer(Properties)]

    data_lines = []
    for line in lines[2:]:
        data_lines.append(list(filter(bool, line.strip().split())))
    data_array = np.array(data_lines)
    used_colomn = 0

    type_array = None
    coords_array = None
    Z_array = None
    force_array = None
    virials = None
    for kv_dict in prop_list:
        if kv_dict['key'] == 'species':
            if kv_dict['datatype'] != 'S':
                raise RuntimeError("datatype for species must be 'S' instead of {}".format(kv_dict['datatype']))
            field_length = int(kv_dict['value'])
            type_array = data_array[:,used_colomn:used_colomn+field_length].flatten()
            used_colomn += field_length
            continue
        elif kv_dict['key'] == 'pos':
            if kv_dict['datatype'] != 'R':
                raise RuntimeError("datatype for pos must be 'R' instead of {}".format(kv_dict['datatype']))
            field_length = int(kv_dict['value'])
            coords_array = data_array[:,used_colomn:used_colomn+field_length]
            used_colomn += field_length
            continue
        elif kv_dict['key'] == 'Z':
            if kv_dict['datatype'] != 'I':
                raise RuntimeError("datatype for pos must be 'R' instead of {}".format(kv_dict['datatype']))
            field_length = int(kv_dict['value'])
            Z_array = data_array[:,used_colomn:used_colomn+field_length].flatten()
            used_colomn += field_length
            continue
        elif kv_dict['key'] == 'force':
            if kv_dict['datatype'] != 'R':
                raise RuntimeError("datatype for pos must be 'R' instead of {}".format(kv_dict['datatype']))
            field_length = int(kv_dict['value'])
            force_array = data_array[:,used_colomn:used_colomn+field_length]
            used_colomn += field_length
            continue
        else:
            raise RuntimeError("unknown field {}".format(kv_dict['key']))

    type_num_dict = OrderedDict()
    atom_type_list = []
    type_map = {}
    temp_atom_max_index = 0
    if type_array is None:
        raise RuntimeError("type_array can't be None type, check .xyz file")
    for ii in type_array:
        if ii not in type_map:
            type_map[ii] = temp_atom_max_index
            temp_atom_max_index += 1
            temp_atom_index = type_map[ii]
            atom_type_list.append(temp_atom_index)
            type_num_dict[ii] = 1
        else:
            temp_atom_index = type_map[ii]
            atom_type_list.append(temp_atom_index)
            type_num_dict[ii] += 1
    type_num_list = []
    for atom_type,atom_num in type_num_dict.items():
        type_num_list.append((atom_type,atom_num))
    type_num_array = np.array(type_num_list)
    if field_dict.get('virial', None):
        virials = np.array(list(filter(bool,field_dict['virial'].split(' ')))).reshape(3,3).astype('float32')
    else:
        virials = None

    cell = np.array(np.array(list(filter(bool,field_dict['lattice'].split(' ')))).reshape(3,3)).astype('float32')
    positions = np.array(coords_array).astype('float32')
    symbols = np.array(type_array).astype(str)

    d = {}
    d['energy'] = np.array(field_dict['energy']).astype('float32')
    if virials is not None:
        d['stress'] = -virials / np.linalg.det(cell)
    d['forces'] = np.array(force_array).astype('float32')

    atoms = Atoms(cell=cell, positions=positions, symbols=symbols, pbc=True)
    calc = SinglePointCalculator(atoms, **d)
    atoms.info.update(d)
    atoms.set_calculator(calc)

    return atoms


def load_nep(filename, ftype="nep"):
    """Read Atoms objects from file

    Args:
        filename (str): Name of the file to read from
        ftype (str): Type of the file to read from

    Returns:
        A list of Atoms objects: frames
    """
    if ftype == "nep":

        frames, n_atoms, has_virial = [], [], []
        with open(filename, 'r') as f:
            line = f.readline()
            n_frames = int(line)
            for _ in range(n_frames):
                line = f.readline()
                print(line)
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
                d['forces'] = np.array(d['forces'])
                atoms = Atoms(cell=cell, positions=positions, symbols=symbols, pbc=True)
                calc = SinglePointCalculator(atoms, **d)
                atoms.info.update(d)
                atoms.set_calculator(calc)
                frames.append(atoms)

    elif ftype == "exyz":

        fin = open(filename, 'r')
        n_atoms = []
        p3 = re.compile(r'^\s*(\d+)\s*')
        block_xyz = []
        while True:
            line = fin.readline()
            if not line:
                break
            if p3.match(line):
                atom_num = int(p3.match(line).group(1))
                n_atoms.append(atom_num)
                lines = []
                lines.append(line)
                for ii in range(atom_num+1):
                    lines.append(fin.readline())
                block_xyz.append(lines)
        fin.close()

        frames = []
        n_frames = len(n_atoms)
        for i in range(n_frames):
            atoms = Proc_block(block_xyz[i])
            frames.append(atoms)

    return frames


def dump_nep(filename, frames, ftype="nep", dvi_virial=False):
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

    if ftype == "nep":

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

    # copy from nep2xyz.py
    elif ftype == "exyz":

        has_virial = ""
        no_virial = ""
        for atoms in frames:

            this_infos = str(int(len(atoms))) + "\n"
            this_infos += "energy=" + str(atoms.info['energy']) + " "
            this_infos += "config_type=FromPyNEP "
            this_infos += "pbc=\"T T T\" "
            if 'stress' in atoms.info:
                #print(atoms.info['stress'])
                #print(len(atoms.info['stress']))
                if len(atoms.info['stress']) == 6:
                    virial = -atoms.info['stress'][[0, 5, 4, 5, 1, 3, 4, 3, 2]] * atoms.get_volume()
                else:
                    virial = -atoms.info['stress'].reshape(-1) * atoms.get_volume()
                this_infos += "virial=\"" + " ".join(list(map(str, virial))) + "\" "
            this_infos += "Lattice=\"" + " ".join(list(map(str, atoms.get_cell().reshape(-1)))) + "\" "
            this_infos += "Properties=species:S:1:pos:R:3:force:R:3\n"

            for j in range(int(len(atoms))):

                s = atoms.get_chemical_symbols()
                p = atoms.get_positions()
                forces = atoms.info['forces']
                this_infos += '{:2} {:>15.8e} {:>15.8e} {:>15.8e} {:>15.8e} {:>15.8e} {:>15.8e}\n'.format(
                               s[j], *p[j], *forces[j])
            if 'stress' in atoms.info:
                has_virial += this_infos
            else:
                no_virial += this_infos

        if dvi_virial == True:
            fho = open(f"has_virial_{filename}", 'w')
            fho.write(has_virial)
            fho.close()
            fno = open(f"no_virial_{filename}", 'w')
            fno.write(no_virial)
            fno.close()
        else:
            fo = open(filename, 'w')
            fo.write(has_virial)
            fo.write(no_virial)
            fo.close()
