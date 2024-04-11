import numpy as np


def inside_starting_region(p_i: np.array, c_sp: np.array, r_sp: float):
    '''
    Function to test whether a point is inside the starting region
    :param p_i: point to be tested
    :param c_sp: center of the region
    :param r_sp: radius of the region
    :return: True if the point if inside the region, False otherwise
    '''
    if np.linalg.norm(p_i - c_sp) < r_sp:
        return True
    else:
        return False


def samples_outside_region(points: np.array, center: np.array, radius: float, n_points: int):
    '''
    :param points:
    :param center:
    :param radius:
    :param n_points:
    :return:
    '''
    n_outside = 0
    for point in points:
        if np.linalg.norm(point - center) > radius:
            n_outside += 1
    if n_outside >= n_points:
        return True
    else:
        return False


def get_groups_from_indexes(indexes):
    '''
    Divide indexes into groups of sequential indexes
    :param indexes: numpy array of indexes
    :return: groups of sequential indexes
    '''
    diffs = np.diff(indexes) != 1
    idx = np.nonzero(diffs)[0] + 1
    groups = np.split(indexes, idx)
    return groups

