from skimage.measure import label, regionprops
from skimage.morphology import binary_dilation
import shared.segmentation as seg
import shared.image as img
import shared.math as mat
import pandas as pd
import math
import matplotlib.pyplot as plt
import napari
import numpy as np

# input parameters
exp_name = '20220204_CtrlAndJQ1'
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20220204_CtrlAndJQ1_zstack/"
sample = 'JQ1'
local_size = 100
rmax = 60
total_FOV = 20
total_z = 16


data = pd.DataFrame(columns=['FOV', 'z', 'nuclear', 'nuclear_centroid', 'nuclear_area', 'nuclear_circ', 'FISH_mean_int',
                             'bg_nuclear', 'bg_FISH', 'g_FISH', 'dg_FISH', 'g_nuclear', 'dg_nuclear'])

for fov in range(total_FOV):
    fov_str = '0%s' % fov if fov < 10 else fov
    for z in range(total_z):
        z_str = '0%s' % z if z < 10 else z
        img_nuclear = plt.imread('%s%s/%s_new_NW_20pos_s%s_z%s_RAW_ch00.tif' % (master_folder, sample, sample, fov_str,
                                                                                z_str), format=None)
        img_FISH = plt.imread('%s%s/%s_new_NW_20pos_s%s_z%s_RAW_ch01.tif' % (master_folder, sample, sample, fov_str,
                                                                             z_str), format=None)

        bg_FISH = seg.get_bg_int([img_FISH])[0]
        bg_nuclear = seg.get_bg_int([img_nuclear])[0]
        print(bg_FISH)
        print(bg_nuclear)

        bg_img_nuclear = np.ones_like(img_nuclear) * int(bg_nuclear)
        img_nuclear_correct = img.image_deduction(img_nuclear, bg_img_nuclear)
        bg_img_FISH = np.ones_like(img_FISH) * int(bg_FISH)
        img_FISH_correct = img.image_deduction(img_FISH, bg_img_FISH)

        img_nuclear_seg = seg.nuclear_seg(img_nuclear)
        nuclear_props = regionprops(img_nuclear_seg)

        for i in range(len(nuclear_props)):
            print("Start image analysis fov %s/%s, z %s/%s, nuclear %s/%s..." % (fov+1, total_FOV, z + 1, total_z,
                                                                                 i + 1, len(nuclear_props)))
            if nuclear_props[i].area * 1.0 / nuclear_props[i].convex_area > 0.9:
                original_centroid = nuclear_props[i].centroid

                position = img.img_local_position(img_nuclear_seg, original_centroid, local_size)
                nuclear_seg = img.img_local_seg(img_nuclear_seg, position, i + 1)
                nuclear = img_nuclear_correct[position[0]:position[1], position[2]:position[3]]
                centroid = regionprops(label(nuclear_seg))[0].centroid
                FISH = img_FISH_correct[position[0]:position[1], position[2]:position[3]]

                nuclear_convex_local = regionprops(label(nuclear_seg))[0].convex_image
                centroid_convex = regionprops(label(nuclear_convex_local))[0].centroid
                nuclear_convex = img.image_paste(FISH, nuclear_convex_local, [int(centroid[0] - centroid_convex[0]),
                                                                              int(centroid[1] - centroid_convex[1])])

                nuclear_area = regionprops(label(nuclear_convex))[0].area
                FISH_mean_int = regionprops(label(nuclear_convex), FISH)[0].mean_intensity
                nuclear_circ = (4 * math.pi * nuclear_area) / (regionprops(label(nuclear_convex))[0].perimeter ** 2)

                _, r, g_FISH, dg_FISH = mat.auto_correlation(FISH, nuclear_seg, rmax)
                _, r, g_nuclear, dg_nuclear = mat.auto_correlation(nuclear, nuclear_seg, rmax)

                data.loc[len(data.index)] = \
                    [fov, z, i + 1, original_centroid, nuclear_area, nuclear_circ, FISH_mean_int, bg_nuclear, bg_FISH,
                     g_FISH, dg_FISH, g_nuclear, dg_nuclear]
            else:
                print("didn't pass convex filter")

data.to_csv('%sauto_correlation_%s.txt' % (master_folder, sample), index=False, sep='\t')

print("END!")