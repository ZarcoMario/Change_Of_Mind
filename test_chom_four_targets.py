import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from align_trajectories import align_trial_start

from change_of_mind_shapely import set_zones_changes_of_mind, get_changes_of_mind_four_targets

dataset = "VR-F2"
path_ = dataset + r"\P01\trial_results.csv"
path_norm_trj = dataset + r"\P01\Normalized_Trajectories"


results = pd.read_csv(path_, usecols=["fin_pos_x", "fin_pos_y"])
fin_pos_xy = results[["fin_pos_x", "fin_pos_y"]].to_numpy()[0:16]
n_baseline_trials = 16

data_bl_trials_dim1, data_bl_trials_dim2 = [], []

for trial_number in range(1, n_baseline_trials + 1):

    path_trial = path_norm_trj + r"\normalized_trajectory_T" + str(trial_number).zfill(3) + ".csv"

    normalized_data = pd.read_csv(path_trial, usecols=['t', 'x', 'y', 'z'])

    t = normalized_data['t'].to_numpy()
    x = normalized_data['x'].to_numpy()
    y = normalized_data['y'].to_numpy()
    z = normalized_data['z'].to_numpy()

    x, y = align_trial_start(x, y, [0, 0.38])

    data_bl_trials_dim1.append(x)
    data_bl_trials_dim2.append(y)

polygons_list, points_list = set_zones_changes_of_mind(np.array(data_bl_trials_dim1), np.array(data_bl_trials_dim2),
                                                       fin_pos_xy, 1.5)

direction = pd.read_csv(path_, usecols=["central_stimulus_direction"])["central_stimulus_direction"].to_numpy()

for trial_number in range(n_baseline_trials + 1, 224  + 1):

    print(trial_number, direction[trial_number - 1])

    path_trial = path_norm_trj + r"\normalized_trajectory_T" + str(trial_number).zfill(3) + ".csv"
    normalized_data = pd.read_csv(path_trial, usecols=['t', 'x', 'y', 'z'])

    t = normalized_data['t'].to_numpy()
    x = normalized_data['x'].to_numpy()
    y = normalized_data['y'].to_numpy()
    z = normalized_data['z'].to_numpy()

    x, y = align_trial_start(x, y, [0, 0.38])

    res = get_changes_of_mind_four_targets(polygons_list, x, y, direction[trial_number - 1], 10, [0, 0.38], 0.005)
    print(res)

    if np.any(res == 1):

        fig = plt.figure()
        gs = GridSpec(1, 1)
        ax = fig.add_subplot(gs[0, 0])
        ax.grid(True)
        ax.plot(x, y, '.')

        for polygon_ in polygons_list:
            x, y = polygon_.exterior.xy
            ax.plot(x, y, 'k-')

        ax.set_xlim(-0.4, 0.4)
        ax.set_ylim(0.0, 0.8)
        plt.show()

