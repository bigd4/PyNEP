# TODO
# KDTree
import numpy as np
from scipy.spatial.distance import cdist


class FarthestPointSample:
    """Farthest point sampler

    Example:

    1. Select points from random 2d points:

    >>> data = np.random.randn(100000, 2)
    >>> selector = FarthestPointSample(min_distance=0.05)
    >>> indices = selector.select(data)

    2. Select atoms with structure descriptors

    Suppose we already have frames to be selected and a NEP calculator. 
    In fact, you can get descriptors by any other method, such as SOAP.

    >>> des = np.array([np.mean(calc.get_property('descriptor', atoms), axis=0) for atoms in frames])
    # Use average of each atomic descriptors to get structure descriptors, shape: (Nframes, Ndescriptor)
    >>> sampler = FarthestPointSample()
    >>> indices = sampler.select(des, [])
    >>> selected_structures = [frames[i] for  i in indices]

    3. Select atoms with atomic latent descriptors
    
    >>> lat = np.concatenate([calc.get_property('latent', atoms) for atoms in frames])
    # shape: (Natoms, Nlatent)
    >>> comesfrom = np.concatenate([i] * len(atoms) for i, atoms in enumerate(frames)])
    # shape: (Natoms, )  the ith data in lat belong to the atoms: frames[comesfrom[i]]
    >>> sampler = FarthestPointSample()
    >>> indices = [comesfrom[i] for i in sampler.select(lat, [])]
    >>> indices = set(indices)  # remove repeated indices because two points may come from the same structure
    >>> selected_structures = [frames[i] for  i in indices]

    """
    def __init__(self, min_distance=0.1, metric='euclidean', metric_para={}):
        """Initial the sampler

        Args:
            min_distance (float, optional): minimum distance between selected data. Defaults to 0.1.
            metric (str, optional): metric of distance between data. 
                Defaults to 'euclidean'. Any metric can be used by 'scipy.spatial.distance.cdist' 
                such as 'cosine', 'minkowski' can also be used.
            metric_para (dict, optional): Extra arguments to metric. 
                Defaults to {}.

            More information about metric can be found: https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.distance.cdist.html
        """
        self.min_distance = min_distance
        self.metric = metric
        self.metric_para = {}

    def select(self, new_data, now_data=[], min_distance=None, min_select=1, max_select=None):
        """Select those data fartheset from given data

        Args:
            new_data (2d list or array): A series of points to be selected
            now_data (2d list or array): Points already in the dataset. 
                Defaults to []. (No existed data)
            min_distance (float, optional): 
                If distance between two points exceeded the minimum distance, stop the selection. 
                Defaults to None (use the self.min_distance)
            min_select (int, optional): Minimal numbers of points to be selected. This may cause
                some distance between points less than given min_distance.
                Defaults to 1.
            max_select (int, optional): Maximum numbers of points to be selected. 
                Defaults to None. (No limitation)

        Returns:
            A list of int: index of selected points
        """
        min_distance = min_distance or self.min_distance
        max_select = max_select or len(new_data)
        to_add = []
        if len(new_data) == 0:
            return to_add
        if len(now_data) == 0:
            to_add.append(0)
            now_data.append(new_data[0])
        distances = np.min(cdist(new_data, now_data, metric=self.metric, **self.metric_para), axis=1)

        while np.max(distances) > min_distance or len(to_add) < min_select:
            i = np.argmax(distances)
            to_add.append(i)
            if len(to_add) >= max_select:
                break
            distances = np.minimum(distances, cdist([new_data[i]], new_data, metric=self.metric)[0])
        return to_add
