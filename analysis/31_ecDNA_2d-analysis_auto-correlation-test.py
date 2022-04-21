import matplotlib.pyplot as plt
from skimage.measure import label, regionprops
import shared.segmentation as seg
from skimage.morphology import binary_dilation, binary_erosion, disk, dilation, medial_axis
import shared.objects as obj
from scipy import ndimage
import pandas as pd
import math
import numpy as np
import shared.image as img
import random
import shared.math as mat
import napari
import shared.dataframe as dat

# input parameters
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20220407_sp8_DMandHSR/DM_singleZ/"
prefix = '20220407_DMandHSR_DM_singleZ'
total_fov = 6
sample = 'DM'
local_size = 150
rmax = 100
radial_interval = 1
radial_max = 120
relative_radial_interval = 0.01
average_image_size = 400

FISH_sum = np.zeros(shape=(average_image_size, average_image_size))
nuclear_sum = np.zeros(shape=(average_image_size, average_image_size))
nuclear_seg_sum = np.zeros(shape=(average_image_size, average_image_size))
center = [average_image_size/2, average_image_size/2]

data = pd.DataFrame(columns=['fov', 'nuclear_label', 'int_thresh', 'g', 'g_value'])

for fov in range(total_fov):
    print("Start analyzing FOV %s/%s" % (fov+1, total_fov))

    # load images
    img_nuclear = plt.imread('%s/%s_s%s_ch00.tif' % (master_folder, prefix, fov), format=None)
    img_FISH = plt.imread('%s/%s_s%s_ch02.tif' % (master_folder, prefix, fov), format=None)
    img_MYC = plt.imread('%s/%s_s%s_ch01.tif' % (master_folder, prefix, fov), format=None)

    # nuclear segmentation
    print("Start segmentation...")
    img_nuclear_seg = seg.nuclear_seg1(img_nuclear, local_factor=99, min_size=9000, max_size=50000)
    img_nuclear_seg = seg.filter_solidity(img_nuclear_seg)
    img_nuclear_seg_convex = seg.obj_to_convex(img_nuclear_seg)
    img_nuclear_seg_convex = dilation(img_nuclear_seg_convex, disk(3))

    # ecDNA segmentation
    _, img_FISH_seg = seg.find_organelle(img_FISH, 'na', extreme_val=18000, bg_val=seg.get_bg_int([img_FISH])[0],
                                         min_size=5, max_size=50000)
    img_FISH_seg = binary_dilation(img_FISH_seg, disk(10))
    for i in range(10):
        img_FISH_seg = binary_erosion(img_FISH_seg)
    img_FISH_seg = ndimage.binary_fill_holes(img_FISH_seg)
    img_FISH_seg = obj.remove_small(img_FISH_seg, min_size=20)

    # filter
    img_nuclear_seg_convex = seg.filter_mean_int(img_nuclear_seg_convex, img_FISH_seg, img_FISH, 20000)
    img_FISH_seg[img_nuclear_seg_convex == 0] = 0

    # props
    print("Start analyzing features...")

    FISH_props = regionprops(img_nuclear_seg_convex, img_FISH)
    nuclear_props = regionprops(img_nuclear_seg_convex, img_nuclear)
    MYC_props = regionprops(img_nuclear_seg_convex, img_MYC)

    for i in range(len(FISH_props)):
        """data = pd.DataFrame(columns=['fov', 'nuclear_label', 'int_thresh', 'dots', 'g', 'g_value'])
        print("Start analyzing nuclear %s/%s" % (i + 1, len(FISH_props)))
        # auto-correlation
        nuclear_centroid = FISH_props[i].centroid
        position = img.img_local_position(img_nuclear_seg_convex, nuclear_centroid, local_size)
        nuclear_seg = img.img_local_seg(img_nuclear_seg_convex, position, i + 1)
        img_FISH_temp = img_FISH.copy()
        FISH = img_FISH_temp[position[0]:position[1], position[2]:position[3]]
        vector = []
        vector_cum_weight = []
        weight = 0
        for int_thresh in np.arange(0, 65535, 3000):
            for dots in np.arange(0, 10000, 500):
                print("test int_thresh: %s, dots: %s" % (int_thresh, dots))
                for m in range(len(nuclear_seg)):
                    for n in range(len(nuclear_seg[0])):
                        if nuclear_seg[m][n] == 1:
                            vector.append([m, n])
                            if FISH[m][n] > int_thresh:
                                weight = weight + FISH[m][n] - int_thresh
                            vector_cum_weight.append(weight)
                random_dot = random.choices(vector, cum_weights=vector_cum_weight, k=dots)
                img_dot = np.zeros_like(nuclear_seg)
                for j in random_dot:
                    img_dot[j[0]][j[1]] = img_dot[j[0]][j[1]] + 1

                _, r, g, dg = mat.auto_correlation(img_dot, nuclear_seg, rmax)
                g_value = (g[1] + g[2] + g[3] + g[4] + g[5])*0.2

                data.loc[len(data.index)] = [fov, i, int_thresh, dots, g, g_value]"""
        dots = 10000
        print("Start analyzing nuclear %s/%s" % (i + 1, len(FISH_props)))
        # auto-correlation
        nuclear_centroid = FISH_props[i].centroid
        position = img.img_local_position(img_nuclear_seg_convex, nuclear_centroid, local_size)
        nuclear_seg = img.img_local_seg(img_nuclear_seg_convex, position, i + 1)
        img_FISH_temp = img_FISH.copy()
        FISH = img_FISH_temp[position[0]:position[1], position[2]:position[3]]
        vector = []
        vector_cum_weight = []
        weight = 0
        for int_thresh in np.arange(0, 65535, 3000):
            print("test int_thresh: %s, dots: %s" % (int_thresh, dots))
            for m in range(len(nuclear_seg)):
                for n in range(len(nuclear_seg[0])):
                    if nuclear_seg[m][n] == 1:
                        vector.append([m, n])
                        if FISH[m][n] > int_thresh:
                            weight = weight + FISH[m][n] - int_thresh
                        vector_cum_weight.append(weight)
            random_dot = random.choices(vector, cum_weights=vector_cum_weight, k=dots)
            img_dot = np.zeros_like(nuclear_seg)
            for j in random_dot:
                img_dot[j[0]][j[1]] = img_dot[j[0]][j[1]] + 1

            _, r, g, dg = mat.auto_correlation(img_dot, nuclear_seg, rmax)
            g_value = (g[1] + g[2] + g[3] + g[4] + g[5]) * 0.2

            data.loc[len(data.index)] = [fov, i, int_thresh, g, g_value]

data.to_csv('%s%s_auto_correlation.txt' % (master_folder, sample), index=False, sep='\t')

print("DONE!")




