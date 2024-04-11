'''

'''

import numpy as np
from scipy.spatial import Delaunay
from utils import get_groups_from_indexes, samples_outside_region

import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec


def _center_fin_pos(fin_pos_xy):
    min_ = np.min(fin_pos_xy[:, 1])
    max_ = np.max(fin_pos_xy[:, 1])
    fin_pos_xy[:, 1] = fin_pos_xy[:, 1] - 0.5 * (min_ + max_)

    fin_pos = np.zeros(fin_pos_xy.shape)

    for i, xy in enumerate(fin_pos_xy):
        if np.abs(xy[0]) > np.abs(xy[1]):
            fin_pos[i, 0] = xy[0]
        else:
            fin_pos[i, 1] = xy[1]

    return fin_pos


def _get_baseline_trial_number_fin_pos(fin_pos: np.array) -> list[np.array]:
    # fin_pos is a numpy array of dimension (n, 2)
    # if the second column is zeros, there are only two targets
    cond_ = fin_pos[:, 1] == 0
    if np.sum(cond_.astype(int)) == fin_pos.shape[0]:
        # Two targets
        # x | -x
        number_baseline_trials = np.full((2, fin_pos.shape[0]), None)
    else:
        # Four Targets
        # x | -x | y | -y
        number_baseline_trials = np.full((4, fin_pos.shape[0]), None)
    # i = 0
    for num, pos in enumerate(fin_pos):
        if np.abs(pos[0]) > np.abs(pos[1]):
            # Target along x
            if pos[0] > 0:
                number_baseline_trials[0, num] = num
            else:
                number_baseline_trials[1, num] = num
        else:
            # Target along y
            if pos[1] > 0:
                number_baseline_trials[2, num] = num
            else:
                number_baseline_trials[3, num] = num
    # Get rid of Nones
    list_number_bst = []
    for i in range(number_baseline_trials.shape[0]):
        row = number_baseline_trials[i, :]
        row = row[row is not None]
        list_number_bst.append(row)

    return list_number_bst


def _get_points_convex_hull(number_baseline_trials: list[np.array],
                            baseline_trials_dim_1: np.array,
                            baseline_trials_dim_2: np.array,
                            n_std: float) -> list[np.array]:

    points_convex_hull = [[] for _ in range(len(number_baseline_trials))]
    for direction, n_bl_trials in enumerate(number_baseline_trials):
        # For each dimension get the mean and std
        # Get the points
        m_dim1 = np.mean(baseline_trials_dim_1[n_bl_trials], axis=0)
        m_dim2 = np.mean(baseline_trials_dim_2[n_bl_trials], axis=0)

        s_dim1 = n_std * np.std(baseline_trials_dim_1[n_bl_trials], axis=0)
        s_dim2 = n_std * np.std(baseline_trials_dim_2[n_bl_trials], axis=0)

        if direction == 0: # RIGHT
            points_convex_hull[direction].append(np.column_stack([m_dim1 - s_dim1, m_dim2 + s_dim2]))
            points_convex_hull[direction].append(np.column_stack([m_dim1 + s_dim1, m_dim2 - s_dim2]))
        elif direction == 1: # LEFT
            points_convex_hull[direction].append(np.column_stack([m_dim1 - s_dim1, m_dim2 - s_dim2]))
            points_convex_hull[direction].append(np.column_stack([m_dim1 + s_dim1, m_dim2 + s_dim2]))

        # TODO: DOUBLE CHECK THIS BY PLOTTING
        elif direction == 2: # UP
            points_convex_hull[direction].append(np.column_stack([m_dim1 - s_dim1, m_dim2 + s_dim2]))
            points_convex_hull[direction].append(np.column_stack([m_dim1 + s_dim1, m_dim2 + s_dim2]))
        else: # DOWN
            points_convex_hull[direction].append(np.column_stack([m_dim1 - s_dim1, m_dim2 - s_dim2]))
            points_convex_hull[direction].append(np.column_stack([m_dim1 + s_dim1, m_dim2 - s_dim2]))


def _get_convex_hull():
    pass


def set_zones_changes_of_mind(data_bl_trials_x, data_bl_trials_y, fin_pos_x, n_std=1.5):
    '''

    :param data_bl_trials_x:
    :param data_bl_trials_y:
    :param fin_pos_x:
    :param n_std:
    :return:
    '''
    # n_std = 1.5

    left_bl_trials = []
    right_bl_trials = []

    for i in range(len(fin_pos_x)):
        if fin_pos_x[i] < 0:
            left_bl_trials.append(i)
        elif fin_pos_x[i] > 0:
            right_bl_trials.append(i)

    x_left_baselines = []
    y_left_baselines = []
    for trial_number in left_bl_trials:
        # Get each one of the trials in the list
        x_left_baselines.append(data_bl_trials_x[trial_number])
        y_left_baselines.append(data_bl_trials_y[trial_number])

    x_right_baselines = []
    y_right_baselines = []
    for trial_number in right_bl_trials:
        # Get each one of the trials in the list
        x_right_baselines.append(data_bl_trials_x[trial_number])
        y_right_baselines.append(data_bl_trials_y[trial_number])

    x_left_mean = np.mean(np.array(x_left_baselines), axis=0)
    y_left_mean = np.mean(np.array(y_left_baselines), axis=0)
    x_left_std = n_std * np.std(np.array(x_left_baselines), axis=0)
    y_left_std = n_std * np.std(np.array(y_left_baselines), axis=0)

    left_points_convex_hull = []
    for x_m, y_m, x_s, y_s in zip(x_left_mean, y_left_mean, x_left_std, y_left_std):
        left_points_convex_hull.append(np.array([x_m + x_s, y_m + y_s]))
        left_points_convex_hull.append(np.array([x_m - x_s, y_m - y_s]))

    left_hull = Delaunay(left_points_convex_hull)

    x_right_mean = np.mean(np.array(x_right_baselines), axis=0)
    y_right_mean = np.mean(np.array(y_right_baselines), axis=0)
    x_right_std = n_std * np.std(np.array(x_right_baselines), axis=0)
    y_right_std = n_std * np.std(np.array(y_right_baselines), axis=0)

    right_points_convex_hull = []
    for x_m, y_m, x_s, y_s in zip(x_right_mean, y_right_mean, x_right_std, y_right_std):
        right_points_convex_hull.append(np.array([x_m - x_s, y_m + y_s]))
        right_points_convex_hull.append(np.array([x_m + x_s, y_m - y_s]))

    right_hull = Delaunay(right_points_convex_hull)

    return left_hull, right_hull, np.array(left_points_convex_hull), np.array(right_points_convex_hull)


def get_change_of_mind(left_hull, right_hull, x, y, n_points=5, center=np.array([0, 0]), radius=0):
    # n_points = 5
    # If right, check left zone
    if x[-1] > 0:
        mask_left = np.full(x.shape, True)
        mask_right = np.full(x.shape, True)
        points_t = np.column_stack((x, y))
        idx_left = left_hull.find_simplex(points_t) < 0
        idx_right = right_hull.find_simplex(points_t) < 0

        mask_left[idx_left] = False
        mask_right[idx_right] = False

        indexes_left = np.arange(0, x.size, 1)
        groups_left = get_groups_from_indexes(indexes_left[mask_left])

        indexes_right = np.arange(0, x.size, 1)
        idx_right = indexes_right[mask_right]

        for idx_left in groups_left:

            # if idx_left.size >= n_points:
            # TODO: Double check this function
            # print(points_t[idx_left])
            if samples_outside_region(points_t[idx_left], center, radius, n_points):

                i = 0
                for element in x[idx_right]:
                    if element in x[idx_left]:
                        i += 1

                if idx_left.size - i >= n_points:
                    return 1
        return 0

    # If left, check right zone
    elif x[-1] < 0:
        mask_right = np.full(x.shape, True)
        mask_left = np.full(x.shape, True)
        points_t = np.column_stack((x, y))
        idx_right = right_hull.find_simplex(points_t) < 0
        idx_left = left_hull.find_simplex(points_t) < 0

        mask_right[idx_right] = False
        mask_left[idx_left] = False

        indexes_right = np.arange(0, x.size, 1)
        groups_right = get_groups_from_indexes(indexes_right[mask_right])

        indexes_left = np.arange(0, x.size, 1)
        idx_left = indexes_left[mask_left]

        for idx_right in groups_right:

            # if idx_right.size >= n_points:
            # TODO: Double check this function
            if samples_outside_region(points_t[idx_right], center, radius, n_points):

                i = 0
                for element in x[idx_left]:
                    if element in x[idx_right]:
                        i += 1

                if idx_right.size - i >= n_points:
                    return 1
        return 0

    else:
        # Case in which the x in zeros
        return 0


