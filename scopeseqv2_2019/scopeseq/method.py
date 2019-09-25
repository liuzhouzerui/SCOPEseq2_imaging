import numpy as np
from scipy.spatial import distance


def register(position, reference, d_th, doublet_removal=False):
    """
    find the nearest position in another dataset
    :param position: data points [X,Y]
    :param reference: reference data points [X,Y]
    :param d_th: the distance threshold of near
    :param doublet_removal: if match to multiple reference set, then decided as not registered.
    :return:
        d_position: the nearest ref data point num
        d_inrange: the distance to the nearest ref data point
    """
    d_tmp = distance.cdist(position, reference, 'euclidean')
    d_position = np.argmin(d_tmp, axis=1)
    d_value = np.min(d_tmp, axis=1)
    if doublet_removal:
        double = (d_tmp < d_th).sum(axis=1) > 1
        d_value[double] = d_th
    d_inrange = (d_value < d_th)
    return d_position, d_inrange


def find_rotation_matrix(x, y):
    """
    :param x: 2d vector
    :param y: 2d vector
    :return: the rotation matrix from y to x
    """
    x = x / np.sqrt(np.sum(x**2))
    y = y / np.sqrt(np.sum(y**2))
    cos_xy = x[0]*y[0]+x[1]*y[1]
    sin_xy = x[0]*y[1]-x[1]*y[0]
    rotation_matrix = np.array([[cos_xy, - sin_xy], [sin_xy, cos_xy]])
    return rotation_matrix


def assign_obc(x, barcode_ref, no_signal_th=None, mode='all'):
    """
    obc num is 0 based
    :param x: probe intensity vector
    :param barcode_ref: reference obc pool
    :param no_signal_th: threshold for no signal
    :param mode: 'all' for iteration, 'max' for only first round
    :return: obc vector. or -1 for unmatched obc
    """
    if min(x) == -1:
        return -1
    if (no_signal_th is not None) and (max(x) < no_signal_th):
        return -1
    x_sorted = np.sort(x)
    x_dif = np.diff(x_sorted)/x_sorted[0:(len(x)-1)]
    x_th = x_sorted[np.argmax(x_dif)]
    barcode = (x > x_th) * 1
    i = 1
    if mode == 'all':
        while i < (len(x)-1):
            if ''.join(barcode.astype("str")) in barcode_ref.values:
                return np.where(barcode_ref.values == ''.join(barcode.astype("str")))[0][0]
            else:
                x_dif[np.argmax(x_dif)] = 0
                x_th = x_sorted[np.argmax(x_dif)]
                barcode = (x > x_th) * 1
                i = i + 1
        if ''.join(barcode.astype("str")) in barcode_ref.values:
            return np.where(barcode_ref.values == ''.join(barcode.astype("str")))[0][0]
        else:
            return -1
    if mode == 'max':
        if ''.join(barcode.astype("str")) in barcode_ref.values:
            return np.where(barcode_ref.values == ''.join(barcode.astype("str")))[0][0]
        else:
            return -1


def find_unique_match_position(target, ref):
    match = np.where(ref == target)[0]
    if match.size == 1:
        return match[0]
    else:
        return -1
