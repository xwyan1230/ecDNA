import numpy as np
from skimage.filters import threshold_otsu
from skimage.morphology import binary_dilation, binary_erosion, disk
import shared.objects as obj


"""
# ---------------------------------------------------------------------------------------------------
# FUNCTIONS for IMAGE PROCESSING
# ---------------------------------------------------------------------------------------------------

img_to_int
    FUNCTION: convert RGB color image into intensity image
    SYNTAX:   img_to_int(img: np.array)

"""


def img_to_int(img: np.array):
    """
    Convert RGB color image into intensity image

    :param img: color image acquired from EVOS-M5000, could be any channel including TRANS
    :return: img: intensity image
    """
    out = img
    if len(np.shape(img)) == 3:
        out = img[:, :, 0] + img[:, :, 1] + img[:, :, 2]
    return out
