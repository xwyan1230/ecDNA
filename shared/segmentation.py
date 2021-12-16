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
    
cell_seg_fluorescent   
    FUNCTION: perform cell segmentation from a fluorescent image
    SYNTAX:   cell_seg_fluorescent(img: np.array, otsu_factor=1.5, maxima_threshold=1, max_size=1800, 
                                   min_size=300, circ_thresh=0.6)

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


def cell_seg_fluorescent(img: np.array, otsu_factor=1.5, maxima_threshold=1, max_size=1800,
                         min_size=300, circ_thresh=0.6):
    """
    Perform cell segmentation from a fluorescent image

    Algorithm description:
    Otsu segmentation followed by watershed and filtering

    :param img: np.array
                    fluorescent image
    :param otsu_factor: float
                    factor multiplied to otsu threshold used to perform otsu thresholding
    :param maxima_threshold: int
                    threshold for identify maxima during watershed
    :param max_size: float
                    maximum allowable object size
    :param min_size: float
                    minimum allowable object size
    :param circ_thresh: float
                    minimum allowable circularity
    :return:
    """
    thresh_val = threshold_otsu(img)
    out = img > thresh_val * otsu_factor
    out = obj.label_watershed(out, maxima_threshold)
    out = obj.label_remove_large(out, max_size)
    out = obj.label_remove_small(out, min_size)
    out = obj.label_remove_low_circ(out, circ_thresh)

    return out


