"""
Mutation
===========================

This example shows how to mutate structures including Rattle, Strain, Swap, ChangeAtomType
"""
from pynep.mutate import Combine, Rattle, Strain, Swap, ChangeAtomType
from ase.io import read, write
from pynep.io import load_nep


frames = read('data.traj', '::10')
print("read {} structures".format(len(frames)))
rattle_mutate = Rattle(probability=0.8,               # rattle probability of each atom
                       rattle_range=1.5,              # atoms will only rattle in this range
                       d_ratio=0.7,                   # min distance controlled by d_ratio * covalent_radii
                       min_dis_mat=None,              # min distance dict such as {('C', 'C'): 1.} 
                                                      # if this is given, d_ratio will be ignored 
                       attempt_number=20)
strain_mutate = Strain(cell_cut=1.,                   # max value of strain matrix
                       sigma=0.1,                     # sigma of gaussian
                       d_ratio=0.7, 
                       min_dis_mat={('C', 'C'): 1.},  # d_ratio=0.7 will be ignored
                       attempt_number=20)

mutate = Combine(rattle_mutate, strain_mutate)        # you can combine a series of mutations
newframes = []
for atoms in frames:
    newframes.extend(mutate.generate(atoms, 10))
print("generate {} mutated structures".format(len(newframes)))
write('new.traj', newframes)


# swap mutation swap atoms with different type, so the formula won't change
frames = load_nep('PbTe.in')[::10]
print("read {} structures".format(len(frames)))

swap_mutate = Swap(swap_probability=0.1,   # max swap number is swap_probability * N_atoms
                   d_ratio=0.7, 
                   min_dis_mat=None, 
                   attempt_number=20)
newframes = []
for atoms in frames:
    newframes.extend(swap_mutate.generate(atoms, 5))
print("generate {} mutated structures".format(len(newframes)))
write('swap.traj', newframes)

# change mutation just randomly change each atom's symbol to another so the formula may change
change_mutate = ChangeAtomType(change_probability=0.5,  # change probability of each atom
                               d_ratio=0.7, 
                               min_dis_mat=None, 
                               attempt_number=20)
newframes = []
for atoms in frames:
    newframes.extend(change_mutate.generate(atoms, 2))
print("generate {} mutated structures".format(len(newframes)))
write('change.traj', newframes)
