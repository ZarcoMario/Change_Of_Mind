'''
align_trajectories_2d is a python version of the methods in MouseTrap R package
https://github.com/PascalKieslich/mousetrap
The original method can be found in
R\align.R
src\trajAlign.cpp -> trajAlign, trajAlign3d
src\RcppExports.cpp
WARNING: Default coordinates in align.R are isotropic.
    Thus, this method needs to be testes if it is to be used.
'''
import numpy as np


def align_trajectories_2d(xs, ys, start_p, end_p, start=True, end=False):  # coordinates=np.array([0, 0, 1, 1])):
    n_ = len(xs)
    xi_mean, yi_mean = 0, 0
    xf_mean, yf_mean = 0, 0
    for x_, y_ in zip(xs, ys):
        xi_mean += x_[0]
        yi_mean += x_[0]
        xf_mean += x_[-1]
        yf_mean += x_[-1]

    xi_mean = xi_mean / n_
    yi_mean = yi_mean / n_

    xf_mean = xf_mean / n_
    yf_mean = yf_mean / n_

    x_aligned = []
    y_aligned = []
    for x_, y_ in zip(xs, ys):
        if start:
            if end:
                x_new = ((x_ - x_[0]) / (x_[-1] - x_[0])) * (end_p[0] - start_p[0]) + start_p[0]
                y_new = ((y_ - y_[0]) / (y_[-1] - y_[0])) * (end_p[1] - start_p[1]) + start_p[1]
            else:
                x_new = ((x_ - x_[0]) / (xf_mean - x_[0])) * (end_p[0] - start_p[0]) + start_p[0]
                y_new = ((y_ - y_[0]) / (yf_mean - y_[0])) * (end_p[1] - start_p[1]) + start_p[1]
        else:
            if end:
                x_new = ((x_ - xi_mean) / (x_[-1] - xi_mean)) * (end_p[0] - start_p[0]) + start_p[0]
                y_new = ((y_ - yi_mean) / (y_[-1] - yi_mean)) * (end_p[1] - start_p[1]) + start_p[1]
            else:
                x_new = ((x_ - xi_mean) / (xf_mean - xi_mean)) * (end_p[0] - start_p[0]) + start_p[0]
                y_new = ((y_ - yi_mean) / (yf_mean - yi_mean)) * (end_p[1] - start_p[1]) + start_p[1]
        x_aligned.append(x_new)
        y_aligned.append(y_new)
    return x_aligned, y_aligned


def align_trial_start(x, y, start):
    return x - x[0] + start[0], y - y[0] + start[1]
