import numpy as np
import matplotlib as mpl

def w_cmap(W, levels=1024, shift=0, invert=False):
    max_color = np.array([0.020, 0.19, 0.38, 1.0])
    mid_color = np.array([1, 1, 1, 1.0])
    if invert:
        min_color = np.array([1, 0.70, 0.87, 1])
        neg_color = np.array([0.4, 0.0, 0.12, 1])
    else:
        min_color = np.array([0.4, 0.0, 0.12, 1])
        neg_color = np.array([1, 0.70, 0.87, 1])
    # get min and max values from Wigner function
    bounds = [W.min(), W.max()]
    # create empty array for RGBA colors
    adjust_RGBA = np.hstack((np.zeros((levels, 3)), np.ones((levels, 1))))
    zero_pos = np.round(levels * np.abs(shift - bounds[0])
                        / (bounds[1] - bounds[0]))
    num_pos = levels - zero_pos
    num_neg = zero_pos - 1
    # set zero values to mid_color
    adjust_RGBA[int(zero_pos)] = mid_color
    # interpolate colors
    for k in range(0, levels):
        if k < zero_pos:
            interp = k / (num_neg + 1.0)
            adjust_RGBA[k][0:3] = (1.0 - interp) * \
                min_color[0:3] + interp * neg_color[0:3]
        elif k > zero_pos:
            interp = (k - zero_pos) / (num_pos + 1.0)
            adjust_RGBA[k][0:3] = (1.0 - interp) * \
                mid_color[0:3] + interp * max_color[0:3]
    # create colormap
    wig_cmap = mpl.colors.LinearSegmentedColormap.from_list('wigner_cmap',
                                                            adjust_RGBA,
                                                            N=levels)
    return wig_cmap
