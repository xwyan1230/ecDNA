import numpy as np
from skimage.filters import threshold_otsu
from skimage.morphology import binary_dilation, binary_erosion, disk
import shared.objects as obj


"""
# ---------------------------------------------------------------------------------------------------
# FUNCTIONS for OBJECT SEGMENTATION
# ---------------------------------------------------------------------------------------------------

cell_seg_trans
    FUNCTION: perform cell segmentation from a transmitted light image
    SYNTAX:   cell_seg_trans(img_trans: np.array, max_size=1000, min_size=100, selem_num=8)

"""


def cell_seg_trans(img_trans: np.array, max_size=1000, min_size=100, selem_num=8):
    """
    Perform cell segmentation from a transmitted light image

    Algorithm description:
    Segment cells by identifying highlighted parts surrounded by dark boundaries. This algorithm will
    potentially miss lots of cells with possible false-positives of spaces between cells which is
    around similar size of a given cell.
    Perform otsu thresholding to identify high intensity region, serial erosion and dilation to better
    close cell regions, size selection to exclude large background region and small debris.
    This method has only been tested for images acquired with EVOS-M5000.

    :param img_trans: np.array
                        transmitted light image
    :param max_size: int, optional, default = 1000
                        maximum size of a cell
    :param min_size: int, optional, default = 100
                        minimum size of a cell
    :param selem_num: int, optional, default = 8
                        number of disk size, tested from 5-9 and 7-9 seem all fine
    :return: out: np.array, binary image
                same size as original transmitted light image
                1: cell region
                0: background
    """

    thresh_val = threshold_otsu(img_trans)
    out = img_trans > thresh_val

    selem = disk(selem_num)
    out = binary_erosion(out, selem)
    out = binary_dilation(out, selem)

    out = obj.remove_large(out, max_size=max_size)
    out = obj.remove_small(out, min_size=min_size)

    return out


