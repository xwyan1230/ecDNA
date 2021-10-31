import numpy as np
from skimage.morphology import remove_small_objects, medial_axis, extrema, binary_dilation, dilation
from skimage.measure import label, regionprops
from skimage.filters import sobel
from skimage import segmentation
import math
import random

"""
# ---------------------------------------------------------------------------------------------------
# FUNCTIONS for 0-AND-1 NP.ARRAY (BINARY IMAGE)
# ---------------------------------------------------------------------------------------------------

remove_small
    FUNCTION: remove objects smaller than the specified size
    SYNTAX:   remove_small(obj: np.array, min_size=10)

remove_large
    FUNCTION: remove objects larger than the specified size
    SYNTAX:   remove_large(obj: np.array, max_size=1000)

get_intensity
    FUNCTION: measure mean intensity time series for all given objects
    SYNTAX:   get_intensity(obj: np.array, pixels_tseries: list)

"""


def remove_small(obj: np.array, min_size=10):
    """
    Remove objects smaller than the specified size.

    Expects obj to be an integer image array with objects labeled with 1, and removes objects
    smaller than min_size.

    :param obj: np.array, 0-and-1
    :param min_size: int, optional (default: 10)
                The smallest allowable object size.
    :return: out: nd.array, 0-and-1, same shape and type as input obj

    """

    obj_bool = np.array(obj, bool)
    obj_mask = remove_small_objects(obj_bool, min_size)
    out = np.zeros_like(obj)
    out[obj_mask] = 1

    return out


def remove_large(obj: np.array, max_size=1000):
    """
    Remove objects larger than the specified size.

    Expects obj to be an integer image array with objects labeled with 1, and removes objects
    larger than max_size.

    :param obj: np.array, 0-and-1
    :param max_size: int, optional (default: 1000)
                The largest allowable object size.
    :return: out: np.array, 0-and-1, same shape and type as input obj
    """

    obj_bool = np.array(obj, bool)
    obj_mask = remove_small_objects(obj_bool, max_size)
    out = obj.copy()
    out[obj_mask] = 0

    return out


def get_intensity(label_obj: np.array, pixels_tseries: list):
    """
    Measure mean intensity time series for all given objects

    Usage examples:
    1) measure bleach spots/ctrl spots intensity series

    :param label_obj: np.array, 0-and-1 object mask
    :param pixels_tseries: list, pixels time series
                e.g. [pixels_t0, pixels_t1, pixels_t2, ...]
    :return: obj_int_tseries: list
                list of intensity time series
    """

    max_t = len(pixels_tseries)
    obj_int_tseries = [[] for _ in range(len(np.unique(label_obj))-1)]

    for t in range(0, max_t):
        # measure mean intensity for objects
        obj_props = regionprops(label_obj, pixels_tseries[t])
        for i in range(len(obj_props)):
            obj_int_tseries[i].append(obj_props[i].mean_intensity)

    return obj_int_tseries
