import matplotlib.pyplot as plt
from skimage.measure import label, regionprops
from skimage.morphology import medial_axis, dilation, erosion
import pandas as pd
import math
import numpy as np
import shared.image as img
import random
import shared.math as mat
import napari
import skimage.io as skio
import shared.dataframe as dat
import shared.objects as obj

# input parameters
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20220325_Natasha_THZ1/"
sample = 'DMSO'
pixel_size = 102  # nm (sp8 confocal 3144x3144:58.7, mischel 2048x2048:102)
cell_avg_size = 10  # um (Colo)
n_nuclear_convex_dilation = 3
nuclear_centroid_searching_range = 25  # pixel
local_size = 100

# SET UP PARAMETERS
total_fov = 1
total_z = 22

# LOAD IMAGE
im_stack_nuclear = skio.imread("%s%s/%s_RAW_ch00_fov3.tif" % (master_folder, sample, sample), plugin="tifffile")
im_stack_DNAFISH = skio.imread("%s%s/%s_RAW_ch01_fov3.tif" % (master_folder, sample, sample), plugin="tifffile")
im_stack_seg = skio.imread("%s%s/%s_seg_fov3.tif" % (master_folder, sample, sample), plugin="tifffile")

data_z = pd.read_csv('%s%s/%s_fov3_z.txt' % (master_folder, sample, sample), na_values=['.'], sep='\t')
data_z['centroid_nuclear'] = [dat.str_to_float(data_z['centroid_nuclear'][i]) for i in range(len(data_z))]

for i in range(len(data_z)):
    print("Analyzing nucleus %s/%s" % (i+1, len(data_z)))
    z = data_z['z'][i]
    label_nuclear = data_z['label_nuclear'][i]
    original_centroid_nuclear = data_z['centroid_nuclear'][i]

    img_nuclear_seg_convex = im_stack_seg[z]
    img_nuclear = im_stack_nuclear[z]
    img_DNAFISH = im_stack_DNAFISH[z]

    position = img.img_local_position(img_nuclear_seg_convex, original_centroid_nuclear, local_size)
    local_nuclear_seg_convex = img.img_local_seg(img_nuclear_seg_convex, position, label_nuclear)
    local_nuclear = img_nuclear.copy()
    local_nuclear = local_nuclear[position[0]:position[1], position[2]:position[3]]
    local_DNAFISH = img_DNAFISH.copy()
    local_DNAFISH = local_DNAFISH[position[0]:position[1], position[2]:position[3]]

    # radial distribution
    local_nuclear_props = regionprops(label(local_nuclear_seg_convex))
    local_nuclear_centroid = local_nuclear_props[0].centroid
    _, local_edge_distance_map = medial_axis(local_nuclear_seg_convex, return_distance=True)
    local_centroid_distance_map = img.distance_map_from_point(local_nuclear_seg_convex,
                                                              local_nuclear_centroid)
    local_centroid_distance_map[local_nuclear_seg_convex == 0] = 0
    local_edge_distance_map[local_nuclear_seg_convex == 0] = -1
    local_relative_r_map = local_centroid_distance_map / (
            local_centroid_distance_map + local_edge_distance_map)

    radial_distribution_relative_r_DNAFISH = \
        img.radial_distribution_from_distance_map(local_nuclear_seg_convex, local_relative_r_map,
                                                  local_DNAFISH, 0.01, 1)
    radial_distribution_relative_r_nuclear = \
        img.radial_distribution_from_distance_map(local_nuclear_seg_convex, local_relative_r_map,
                                                  local_nuclear, 0.01, 1)

    radial_distribution_relative_r_DNAFISH_smooth = dat.list_smooth(radial_distribution_relative_r_DNAFISH, 3)
    radial_distribution_relative_r_nuclear_smooth = dat.list_smooth(radial_distribution_relative_r_nuclear, 3)

    radial_subtract = np.array(radial_distribution_relative_r_DNAFISH_smooth)-np.array(radial_distribution_relative_r_nuclear_smooth)
    radial_subtract_center = np.mean(radial_subtract[0:40])
    radial_subtract_edge = np.mean(radial_subtract[40:80])

    # angle distribution
    local_angle_map = img.angle_map_from_point(local_nuclear_seg_convex, local_nuclear_centroid)
    angle_distribution_DNAFISH = img.radial_distribution_from_distance_map(local_nuclear_seg_convex, local_angle_map,
                                                                           local_DNAFISH, 1, 360)
    angle_distribution_nuclear = img.radial_distribution_from_distance_map(local_nuclear_seg_convex, local_angle_map,
                                                                           local_nuclear, 1, 360)

    angle_distribution_DNAFISH_smooth = dat.list_circle_smooth(angle_distribution_DNAFISH, 7)
    angle_distribution_nuclear_smooth = dat.list_circle_smooth(angle_distribution_nuclear, 7)
    angle_distribution_DNAFISH_smooth_centered, angle_distribution_nuclear_smooth_centered = \
        dat.list_peak_center_with_control(angle_distribution_DNAFISH_smooth, angle_distribution_nuclear_smooth)

    viewer = napari.view_image(local_DNAFISH)
    napari.run()
    ax = plt.plot(angle_distribution_DNAFISH_smooth_centered, c='g')
    plt.plot(angle_distribution_nuclear_smooth_centered, c='b')
    plt.show()


    # viewer = napari.view_image(local_nuclear, name='nuclear')
    # viewer.add_image(local_DNAFISH, name='FISH')
    # napari.run()