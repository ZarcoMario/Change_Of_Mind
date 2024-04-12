import numpy as np
from scipy.spatial import Delaunay
from utils import get_groups_from_indexes, samples_outside_region


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
        # Points inside
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


