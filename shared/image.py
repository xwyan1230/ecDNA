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

img_local_position
    FUNCTION: determine position of local image
    SYNTAX:   img_local_position(img: np.array, centroid: tuple, local_size: int)

img_local_seg
    FUNCTION: generate local segmented image
    SYNTAX:   img_local_seg(img: np.array, position: list, label_number: int)

edge_from_seg
    FUNCTION: detect edge from 0,1 image
    SYNTAX:   edge_from_seg(img: np.array)

image_paste
    FUNCTION: paste image based on reference distance
    SYNTAX:   image_paste(paste_to_size: np.array, paste_from_img: np.array, distance: list)

image_deduction
    FUNCTION: image deduction (returns img1-img2, negative value replace as 0)
    SYNTAX:   image_deduction(img1: np.array, img2: np.array)
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


def img_local_position(img: np.array, centroid: tuple, local_size: int):
    """
    Determine position of local image

    Check if local image is within the original image, otherwise, substitute the boundary with original image

    :param img: np.array, original image
    :param centroid: tuple, center of the new local image
    :param local_size: int, half size of the new local image
    :return: return the four boundary of the new image
    """
    shape = img.shape
    left = int(centroid[0] - local_size) if centroid[0] > local_size else 0
    right = int(centroid[0] + local_size) if centroid[0] + local_size < shape[0] else shape[0]
    top = int(centroid[1] - local_size) if centroid[1] > local_size else 0
    bottom = int(centroid[1] + local_size) if centroid[1] + local_size < shape[1] else shape[1]

    return left, right, top, bottom


def img_local_seg(img: np.array, position: list, label_number: int):
    """
    Generate local segmented image

    :param img: np.array, original segmented image
    :param position: list, position of local segmented image
    :param label_number: int, number of label in original segmented image
    :return:
    """
    temp = img[position[0]:position[1], position[2]:position[3]]
    out = np.zeros_like(temp, dtype=float)
    out[temp == label_number] = 1

    return out


def edge_from_seg(img: np.array):
    """
    Detect edge from 0,1 image

    Purpose: the boundary of region 1 will be labeled as 1

    :param img: np.array, original 0, 1 image
    :return:
    """
    shape = img.shape
    out = np.zeros_like(img, dtype=float)
    for m in range(shape[0]-2):
        for n in range(shape[1]-2):
            if (img[m+1, n+1] == 1) & ((img[m, n] == 0)|(img[m, n+1] == 0)|(img[m, n+2] == 0)|(img[m+1, n] == 0)|
                                       (img[m+1, n+2] == 0)|(img[m+2, n] == 0)|(img[m+2, n+1] == 0)|(img[m+2, n+2] == 0)):
                out[m+1, n+1] = 1

    return out


def image_paste(paste_to_size: np.array, paste_from_img: np.array, distance: list):
    """
    Paste image based on reference distance

    :param paste_to_size: np.array, same size image as output
    :param paste_from_img: np.array, image to be pasted
    :param distance: [x,y] reference distance
    :return:
    """
    paste_to_img = np.zeros_like(paste_to_size)
    for i in range(paste_from_img.shape[0]):
        for j in range(paste_from_img.shape[1]):
            if paste_from_img[i][j] != 0:
                paste_to_img[i+distance[0]][j+distance[1]] = paste_from_img[i][j]

    return paste_to_img


def image_deduction(img1: np.array, img2: np.array):
    """
    Image deduction (returns img1-img2, negative value replace as 0)

    :param img1: np.array
    :param img2: np.array
    :return:
    """
    out = np.zeros_like(img1)
    for i in range(len(img1)):
        for j in range(len(img1[i])):
            if img1[i][j] > img2[i][j]:
                out[i][j] = img1[i][j] - img2[i][j]

    return out
