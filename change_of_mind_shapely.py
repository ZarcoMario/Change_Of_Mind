import numpy as np
from shapely.geometry import Polygon
from shapely import Point
from geopandas import GeoSeries

from utils import get_groups_from_indexes, samples_outside_region

import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec


def _center_fin_pos(fin_pos_xy: np.array):
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
    '''
    This method assumes a '+'-like arrangement of the targets in the case of 4 of them
    :param fin_pos:
    :return:
    '''
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
        row = row[row != None].astype(int)
        list_number_bst.append(row)

    return list_number_bst


def _get_points_and_polygon(number_baseline_trials: list[np.array],
                            baseline_trials_dim_1: np.array,
                            baseline_trials_dim_2: np.array,
                            n_std: float) -> tuple[list, list]:
    '''

    :param number_baseline_trials:
    :param baseline_trials_dim_1:
    :param baseline_trials_dim_2:
    :param n_std:
    :return:
    '''
    points_list = []
    polygons_list = []
    for direction, n_bl_trials in enumerate(number_baseline_trials):

        m_dim1 = np.mean(baseline_trials_dim_1[n_bl_trials], axis=0)
        m_dim2 = np.mean(baseline_trials_dim_2[n_bl_trials], axis=0)

        s_dim1 = n_std * np.std(baseline_trials_dim_1[n_bl_trials], axis=0)
        s_dim2 = n_std * np.std(baseline_trials_dim_2[n_bl_trials], axis=0)

        # The limits used for RIGHT AND LEFT come from a MATLAB file Chris provided (afaik)
        # I did not have time to verify the limist for UP AND DOWN. So, change them appropriately
        if direction == 0: # RIGHT
            # NOTICE THAT THE SECOND ARRAY HAS BEEN FLIPPED SO THAT THE POINTS FOLLOW AN ODER TO CREATE THE POLYGON
            points_ = np.append(np.column_stack((m_dim1 - s_dim1, m_dim2 + s_dim2)),
                                np.column_stack((m_dim1 + s_dim1, m_dim2 - s_dim2))[::-1], axis=0)
            points_list.append(points_)
            polygons_list.append(Polygon(points_))

        elif direction == 1: # LEFT
            points_ = np.append(np.column_stack((m_dim1 - s_dim1, m_dim2 - s_dim2)),
                                np.column_stack((m_dim1 + s_dim1, m_dim2 + s_dim2))[::-1], axis=0)
            points_list.append(points_)
            polygons_list.append(Polygon(points_))

        elif direction == 2: # UP
            points_ = np.append(np.column_stack((m_dim1 - s_dim1, m_dim2 + s_dim2)),
                                np.column_stack((m_dim1 + s_dim1, m_dim2 + s_dim2))[::-1], axis=0)
            points_list.append(points_)
            polygons_list.append(Polygon(points_))

        elif direction == 3: # DOWN
            points_ = np.append(np.column_stack((m_dim1 - s_dim1, m_dim2 - s_dim2)),
                                np.column_stack((m_dim1 + s_dim1, m_dim2 - s_dim2))[::-1], axis=0)
            points_list.append(points_)
            polygons_list.append(Polygon(points_))

        # fig = plt.figure()
        # gs = GridSpec(1, 1)
        #
        # ax = fig.add_subplot(gs[0, 0])
        # ax.grid(True)
        # ax.plot(m_dim1, m_dim2, '.')
        # # ax.set_xlim(-0.4, 0.4)
        # # ax.set_ylim(0.6, 1.4)
        # points = points_list[direction]
        # ax.plot(points[:, 0], points[:, 1], '.')
        # for blt_dim1, blt_dim2 in zip(baseline_trials_dim_1[n_bl_trials], baseline_trials_dim_2[n_bl_trials]):
        #     ax.plot(blt_dim1, blt_dim2, '-', alpha=0.5, color='grey')
        #
        # x, y = polygons_list[direction].exterior.xy
        # ax.plot(x, y)
        # plt.show()

    return points_list, polygons_list


def set_zones_changes_of_mind(baseline_trials_dim_1: np.array,
                              baseline_trials_dim_2: np.array,
                              fin_pos: np.array,
                              n_std: float = 1.5):
    fin_pos = _center_fin_pos(fin_pos)

    # Call _get_baseline_trial_number_fin_pos
    number_baseline_trials = _get_baseline_trial_number_fin_pos(fin_pos)

    # Call _get_points_and_polygon
    points_list, polygons_list = _get_points_and_polygon(number_baseline_trials,
                                                         baseline_trials_dim_1,
                                                         baseline_trials_dim_2,
                                                         n_std)
    return polygons_list, points_list


def get_changes_of_mind_two_targets(polygons_list: list[Polygon],
                                    data_dim1: np.array,
                                    data_dim2: np.array,
                                    n_points: int = 10,
                                    center: np.array = np.array([0, 0]),
                                    radius: float = 0) -> int:
    # s = GeoSeries(map(Point, zip(baseline_trials_dim_1[5], baseline_trials_dim_2[5])))
    # res = s.within(polygon)
    # ax.plot(baseline_trials_dim_1[5], baseline_trials_dim_2[5], 'r.')
    # print(res)

    # If right, check left zone
    if data_dim1[-1] > 0:

        s = GeoSeries(map(Point, zip(data_dim1, data_dim2)))
        mask_right = s.within(polygons_list[0])
        mask_left = s.within(polygons_list[1])

        indexes_left = np.arange(0, data_dim1.size, 1)
        groups_left = get_groups_from_indexes(indexes_left[mask_left])

        indexes_right = np.arange(0, data_dim1.size, 1)
        idx_right = indexes_right[mask_right]

        for idx_left in groups_left:

            if samples_outside_region(np.column_stack((data_dim1, data_dim2))[idx_left],
                                      center, radius, n_points):

                i = 0
                for element in data_dim1[idx_right]:
                    if element in data_dim1[idx_left]:
                        i += 1

                if idx_left.size - i >= n_points:
                    return 1

        return 0

    # If left, check right zone
    elif data_dim1[-1] < 0:

        s = GeoSeries(map(Point, zip(data_dim1, data_dim2)))
        mask_right = s.within(polygons_list[0])
        mask_left = s.within(polygons_list[1])

        indexes_right = np.arange(0, data_dim1.size, 1)
        groups_right = get_groups_from_indexes(indexes_right[mask_right])

        indexes_left = np.arange(0, data_dim1.size, 1)
        idx_left = indexes_left[mask_left]

        for idx_right in groups_right:

            if samples_outside_region(np.column_stack((data_dim1, data_dim2))[idx_right],
                                      center, radius, n_points):

                i = 0
                for element in data_dim1[idx_left]:
                    if element in data_dim1[idx_right]:
                        i += 1

                if idx_right.size - i >= n_points:
                    return 1

        return 0

    else:

        return 0


def _change_of_mind(target_direction, other_directions, polygons_list,
                    data_dim1, data_dim2,
                    n_points,
                    center, radius):
    
    s = GeoSeries(map(Point, zip(data_dim1, data_dim2)))
    indexes = np.arange(0, data_dim1.size, 1)
    idx_target_direction = indexes[s.within(polygons_list[target_direction])]

    groups_other_directions = []
    for direction in other_directions:
        groups_other_directions.append(
            get_groups_from_indexes(indexes[s.within(polygons_list[direction])])
        )

    res = np.zeros(other_directions.size)
    j = 0
    for group in groups_other_directions:
        for idx_other_direction in group:

            if samples_outside_region(np.column_stack((data_dim1, data_dim2))[idx_other_direction],
                                      center, radius, n_points):
                i = 0
                for element in data_dim1[idx_target_direction]:
                    if element in data_dim1[idx_other_direction]:
                        i += 1

                if idx_other_direction.size - i >= n_points:
                    # return 1
                    res[j] = 1
                    break

        j += 1

    return res


def get_changes_of_mind_four_targets(polygons_list: list[Polygon],
                                    data_dim1: np.array,
                                    data_dim2: np.array,
                                    data_direction: str,
                                    n_points: int = 5,
                                    center: np.array = np.array([0, 0]),
                                    radius: float = 0) -> int:

    if data_direction == "right":
        res = _change_of_mind(0, np.array([1]), polygons_list, data_dim1, data_dim2, n_points, center, radius)
        print(res)
    elif data_direction == "left":
        res = _change_of_mind(1, np.array([0]), polygons_list, data_dim1, data_dim2, n_points, center, radius)
        print(res)
    elif data_direction == "up":
        pass
    elif data_direction == "down":
        pass
    else:
        print("Incorrect Data Direction")
        return -1





