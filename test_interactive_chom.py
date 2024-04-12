import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from align_trajectories import align_trial_start

from change_of_mind_shapely import _get_baseline_trial_number_fin_pos, _get_points_and_polygon
from change_of_mind_shapely import set_zones_changes_of_mind, get_changes_of_mind_two_targets
from change_of_mind_shapely import get_changes_of_mind_four_targets

dataset = "VR-S1"
path_ = dataset + r"\P01\trial_results.csv"
path_norm_trj = dataset + r"\P01\Normalized_Trajectories"


results = pd.read_csv(path_, usecols=["fin_pos_x", "fin_pos_y"])
#print(results)
fin_pos_xy = results[["fin_pos_x", "fin_pos_y"]].to_numpy()[0:16]
# print(fin_pos_xy)
# min_ = np.min(fin_pos_xy[:, 1])
# max_ = np.max(fin_pos_xy[:, 1])
# fin_pos_xy[:, 1] = fin_pos_xy[:, 1] - 0.5 * (min_ + max_)

n_baseline_trials = 16

data_bl_trials_dim1, data_bl_trials_dim2 = [], []

for trial_number in range(1, n_baseline_trials + 1):

    path_trial = path_norm_trj + r"\normalized_trajectory_T" + str(trial_number).zfill(3) + ".csv"
    # print(path_trial)
    normalized_data = pd.read_csv(path_trial, usecols=['t', 'x', 'y', 'z'])

    t = normalized_data['t'].to_numpy()
    x = normalized_data['x'].to_numpy()
    y = normalized_data['y'].to_numpy()
    z = normalized_data['z'].to_numpy()

    x, z = align_trial_start(x, z, [0, 0.38])

    data_bl_trials_dim1.append(x)
    data_bl_trials_dim2.append(z)

polygons_list, points_list = set_zones_changes_of_mind(np.array(data_bl_trials_dim1), np.array(data_bl_trials_dim2),
                                                       fin_pos_xy, 1.5)

for trial_number in range(249, 249 + 1):

    path_trial = path_norm_trj + r"\normalized_trajectory_T" + str(trial_number).zfill(3) + ".csv"
    normalized_data = pd.read_csv(path_trial, usecols=['t', 'x', 'y', 'z'])

    t = normalized_data['t'].to_numpy()
    x = normalized_data['x'].to_numpy()
    y = normalized_data['y'].to_numpy()
    z = normalized_data['z'].to_numpy()

    x, z = align_trial_start(x, z, [0, 0.38])
    #
    # chom = get_changes_of_mind_two_targets(polygons_list, x, z, 10, [0, 0.38], 0.005)
    # if chom:
    #     print(trial_number)
    #     plt.plot(x, z, '.')
    #     plt.show()

    if x[-1] < 0:
        get_changes_of_mind_four_targets(polygons_list, x, z, "left", 10, [0, 0.38], 0.005)
        plt.plot(x, z, '.')
        plt.show()

