# TODO
# KDTree
# atomic descriptor
import numpy as np
from scipy.spatial.distance import cdist


class FarthestPointSample:
    def __init__(self, min_distance=0.1, metric='euclidean'):
        self.min_distance = min_distance
        self.metric = metric

    def select(self, new_data, now_data=[], min_distance=None, min_select=1, max_select=None):
        min_distance = min_distance or self.min_distance
        max_select = max_select or len(new_data)
        to_add = []
        if len(new_data) == 0:
            return to_add
        if len(now_data) == 0:
            to_add.append(0)
            now_data.append(new_data[0])
        distances = np.min(cdist(new_data, now_data, metric=self.metric), axis=1)

        while np.max(distances) > min_distance or len(to_add) < min_select:
            i = np.argmax(distances)
            to_add.append(i)
            if len(to_add) >= max_select:
                break
            distances = np.minimum(distances, cdist([new_data[i]], new_data, metric=self.metric)[0])
        return to_add
