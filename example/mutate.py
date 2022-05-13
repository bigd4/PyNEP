from pynep.mutate import Combine, Rattle, Strain
from ase.io import read, write


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
