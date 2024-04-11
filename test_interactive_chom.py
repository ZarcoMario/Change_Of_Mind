import pandas as pd
import numpy as np
from change_of_mind import _get_baseline_trial_number_fin_pos


path_ = r"C:\Users\mzar066\Downloads\F2_P01\S001\trial_results.csv"
results = pd.read_csv(path_, usecols=["fin_pos_x", "fin_pos_y"])
#print(results)
fin_pos_xy = results[["fin_pos_x", "fin_pos_y"]].to_numpy()[0:16]
# print(fin_pos_xy)
min_ = np.min(fin_pos_xy[:, 1])
max_ = np.max(fin_pos_xy[:, 1])
fin_pos_xy[:, 1] = fin_pos_xy[:, 1] - 0.5 * (min_ + max_)
# print(fin_pos_xy)

fin_pos = np.zeros(fin_pos_xy.shape)

trial_num = 1

for xy in fin_pos_xy:
    if np.abs(xy[0]) > np.abs(xy[1]):
        fin_pos[trial_num - 1, 0] = xy[0]
    else:
        fin_pos[trial_num - 1, 1] = xy[1]
    trial_num += 1
print(fin_pos)

_get_baseline_trial_number_fin_pos(fin_pos)