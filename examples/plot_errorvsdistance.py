"""
Relationship between error and latent distance
===============================================

Distances in latent space does not provide an error estimate in the units of the 
property being predicted. So we calibrate the error estimate by ftting the predictive variance to 
a simple conditional Gaussian distribution of the error similar to 
`A quantitative uncertainty metric controls error in neural network-driven chemical discovery 
<https://doi.org/10.1039/C9SC02298H>`_

Here we took carbon for example:
"""

from cProfile import label
from pynep.io import load_nep
from pynep.calculate import NEP
import numpy as np
from scipy.spatial.distance import cdist
from scipy.optimize import minimize
import matplotlib.pyplot as plt


calc = NEP('C_2022_NEP3.txt')
train_data = load_nep('train.in')
test_data = load_nep('test.in')
train_lat = np.concatenate([calc.get_property('latent', atoms) for atoms in train_data])
#%%
# We use 10% testset to calculate the force error and the min distance to train data of each atoms

def split_dataset(dataset, ratio=0.1):
    indices = np.arange(len(dataset))
    np.random.shuffle(indices)
    num = int(len(dataset) * ratio)
    return dataset[:num], dataset[num:]
d1, d2 = split_dataset(test_data)

d, e = [], []
for atoms in d1:
    nep_forces = calc.get_forces(atoms)
    dft_forces = atoms.info['forces']
    error = np.sqrt(np.sum((nep_forces - dft_forces) ** 2, axis=1)).tolist()
    lat = calc.get_property('latent', atoms)
    mind = np.min(cdist(lat, train_lat), axis=1).tolist()
    d.extend(mind)
    e.extend(error)
d = np.array(d)           # min distances to trainset
e = np.array(e)           # errors between nep and dft

#%%
# Here we suppose :math:`\sigma = a + b \cdot d + c \cdot d^2`, and a, b, c > 0. 
# Maximize :math:`\mathcal{N}(\epsilon \vert 0, \sigma)` is same to minimize 
# :math:`2 \cdot ln(\sigma) + d^2/\sigma^2` 

def loss(args):
    a, b, c = args
    sigma = np.abs(a + b * d + c * d ** 2)
    l = np.sum(2 * np.log(sigma) + (e / sigma) ** 2)
    return l

cons = (
    {'type': 'ineq', 'fun': lambda x: x - 0.0001},
)

x0 = np.array([1, 1, 1])
res = minimize(loss, x0, constraints=cons, method='SLSQP')

#%%
# results
d_, e_ = [], []
for atoms in d2[::10]:
    nep_forces = calc.get_forces(atoms)
    dft_forces = atoms.info['forces']
    error = np.sqrt(np.sum((nep_forces - dft_forces) ** 2, axis=1)).tolist()

    lat = calc.get_property('latent', atoms)
    mind = np.min(cdist(lat, train_lat), axis=1).tolist()

    d_.extend(mind)
    e_.extend(error)


d_ = np.array(d_)
e_ = np.array(e_)
sigma_ = res.x[0] + res.x[1] * d_ + res.x[2] * d_ ** 2
s1 = len(e_[e_ < sigma_]) / len(e_)
s2 = len(e_[e_ < 2 * sigma_]) / len(e_)
plt.scatter(d_, e_, s=4)

dd = np.linspace(0, max(d_) + 0.1, 100)
sigma = res.x[0] + res.x[1] * dd + res.x[2] * dd ** 2
plt.fill_between(dd, sigma, facecolor='blue', alpha=0.3,
                 label="{:.2%}".format(s1))
plt.fill_between(dd, 2 * sigma, sigma, facecolor='yellow', alpha=0.3,
                 label="{:.2%}".format(s2 - s1))
plt.legend()
