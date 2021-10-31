import numpy as np
from skimage.filters import threshold_otsu
from skimage.morphology import binary_dilation, binary_erosion, disk


"""
# ---------------------------------------------------------------------------------------------------
# FUNCTIONS for OBJECT SEGMENTATION
# ---------------------------------------------------------------------------------------------------

remove_small
    FUNCTION: remove objects smaller than the specified size
    SYNTAX:   remove_small(obj: np.array, min_size=10)

remove_large
    FUNCTION: remove objects larger than the specified size
    SYNTAX:   remove_large(obj: np.array, max_size=1000)
"""


def cell_seg_trans(img_trans):

    thresh_val = threshold_otsu(img_trans)
    temp_img = img_trans > thresh_val

    s8 = disk(8)
    temp_img = binary_erosion(temp_img, s8)
    temp_img = binary_dilation(temp_img, s8)

    cell =


