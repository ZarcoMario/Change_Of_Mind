import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import numpy as np
import pandas as pd
from scipy.spatial import ConvexHull
from align_trajectories import align_trial_start
from change_of_mind import set_zones_changes_of_mind


class InteractiveTrajectories:

    def __init__(self, path_folder, file_name, l_bound, u_bound, plane, path_results, path_save_folder, n_targets):
        self.path_folder = path_folder # Path to normalized trajectories
        self.file_name = file_name # Generic name of the file
        self.l_bound = l_bound
        self.u_bound = u_bound
        self.plane = plane
        self.path_results = path_results # Path to folder that contains the final position of the tracker
        self.path_save_folder = path_save_folder
        self.trials = range(self.l_bound, self.u_bound + 1, 1)

        # Hard-coding some parameters here
        # Initial coordinate to align the trajectories
        self.n_std = 1.5
        # Initial coordinate to align the trajectories
        self.start_p = np.array([0, 0.38])

        self.baseline_trials_dim_1 = []
        self.baseline_trials_dim_2 = []
        self.fin_pos = np.zeros((u_bound - l_bound + 1, 2))
        try:
            assert n_targets in [2, 4], "The number of targets must be 2 or 4"
        except AssertionError as msg:
            print(msg)
        self.n_targets = n_targets

        self.fig = plt.figure(figsize=(16, 8))
        gs = GridSpec(1, 1)
        self.ax = self.fig.add_subplot(gs[0, 0])
        self.ax.grid(True)
        self.lines = []
        self.lines_color = []
        self.lines_convexhull = []
        self.picked = np.full(u_bound - l_bound + 1, True)
        self.plot_trajectories()
        self.connect()
        plt.show()

    def plot_trajectories(self):
        # Get each trajectory
        for trial_num in self.trials:
            data = pd.read_csv(self.path_folder + self.file_name + str(trial_num).zfill(3) + ".csv",
                               usecols=[self.plane[0], self.plane[1]])
            data_dim_1, data_dim_2 = data[self.plane[0]].to_numpy(), data[self.plane[1]].to_numpy()

            # Align Trajectories
            data_dim_1, data_dim_2 = align_trial_start(data_dim_1, data_dim_2, self.start_p)

            self.baseline_trials_dim_1.append(data_dim_1)
            self.baseline_trials_dim_2.append(data_dim_2)

            line, = self.ax.plot(data_dim_1, data_dim_2, '.-', alpha=0.5)
            line.set_picker(True)
            line.set_pickradius(5)
            self.lines.append(line)
            self.lines_color.append(line.get_color())

            if self.n_targets == 2:
                # Two-target version (assuming the change along x on Unity)
                self.fin_pos[trial_num - 1, 0] = data_dim_1[-1]
            else:
                # Four-target condition
                if np.abs(data_dim_1[-1]) > np.abs(data_dim_2[-1]):
                    # Target along x
                    self.fin_pos[trial_num - 1, 0] = data_dim_1[-1]
                elif np.abs(data_dim_1[-1]) < np.abs(data_dim_2[-1]):
                    # Target along y
                    self.fin_pos[trial_num - 1, 1] = data_dim_2[-1]
                else:
                    # I am not checking validity here. Is there any other odd case?
                    pass

        # Convert to numpy so that they can be masked
        self.baseline_trials_dim_1 = np.array(self.baseline_trials_dim_1)
        self.baseline_trials_dim_2 = np.array(self.baseline_trials_dim_2)

        if self.n_targets == 2:

            # Update the zones for changes of mind
            left_hull, right_hull, left_points_convex_hull, right_points_convex_hull = set_zones_changes_of_mind(
                self.baseline_trials_dim_1, self.baseline_trials_dim_2, self.fin_pos[:, 0], n_std=self.n_std)

            # Redraw the zones
            # Plot Left Convex Hull
            for simplex in ConvexHull(left_points_convex_hull).simplices:
                line, = self.ax.plot(left_points_convex_hull[simplex, 0], left_points_convex_hull[simplex, 1], 'k-')
                self.lines_convexhull.append(line)

            # Plot Right Convex Hull
            for simplex in ConvexHull(right_points_convex_hull).simplices:
                line, = self.ax.plot(right_points_convex_hull[simplex, 0], right_points_convex_hull[simplex, 1], 'k-')
                self.lines_convexhull.append(line)

        return

    def connect(self):
        self.fig.canvas.mpl_connect('pick_event', self.on_click)
        self.fig.canvas.mpl_connect('close_event', self.on_close)
        return

    def on_click(self, event):
        for i, line in enumerate(self.lines):
            if event.artist == line:
                if self.picked[i]:
                    line.set_color("grey")
                    self.picked[i] = False
                else:
                    line.set_color(self.lines_color[i])
                    self.picked[i] = True

        if self.n_targets == 2:

            # Update the zones for changes of mind
            left_hull, right_hull, left_points_convex_hull, right_points_convex_hull = set_zones_changes_of_mind(
                self.baseline_trials_dim_1[self.picked], self.baseline_trials_dim_2[self.picked],
                self.fin_pos[self.picked, 0], n_std=self.n_std)

            # Redraw the zones
            for line in self.lines_convexhull:
                line.remove()
            self.lines_convexhull = []

            # Plot Left Convex Hull
            for simplex in ConvexHull(left_points_convex_hull).simplices:
                line, = self.ax.plot(left_points_convex_hull[simplex, 0], left_points_convex_hull[simplex, 1], 'k-')
                self.lines_convexhull.append(line)

            # Plot Right Convex Hull
            for simplex in ConvexHull(right_points_convex_hull).simplices:
                line, = self.ax.plot(right_points_convex_hull[simplex, 0], right_points_convex_hull[simplex, 1], 'k-')
                self.lines_convexhull.append(line)

        self.fig.canvas.draw()
        return

    def on_close(self, event):
        # Save discarded trials
        # dict_res = {
        #     'trial_num': self.trials,
        #     'incorrect': self.picked.astype(int)
        # }
        # pd.DataFrame.from_dict(dict_res).to_csv(self.path_folder + r"\movement_termination_incorrect.csv", index=False)
        pass


if __name__ == "__main__":

    path_folder_ = r"D:\mzar066\Jobs 2024\Research_Assistant\VR-S1\Results_VR_S1\Normalized_Trajetcories\P04"
    file_name_ = r"\normalized_trajectory_T"
    lower_bound = 1
    upper_bound = 16
    plane_ = "xz"
    interactive = InteractiveTrajectories(path_folder_, file_name_, lower_bound, upper_bound, plane_, "", "", 2)
