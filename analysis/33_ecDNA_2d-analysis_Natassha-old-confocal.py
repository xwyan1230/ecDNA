from skimage.measure import regionprops
import shared.segmentation as seg
import pandas as pd
import math
import numpy as np
import shared.image as img
import skimage.io as skio
import tifffile as tif

# input parameters
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/211102_3hr_JQ1washout/"
prefix = '211102_COLODM'
total_fov = 2
sample = '3hrJQ1_3hr1uMtriptolide_WO'
pixel_size = 40  # nm (Zeiss confocal scope)
cell_avg_size = 10  # um (Colo)
nuclear_size_range = [0.6, 1.5]  # used to filter nucleus
total_fov = 2

# SET UP PARAMETERS
# segmentation
local_factor_nuclear = 299  # needs to be odd number, rmax if (rmax % 2 == 1) else rmax+1
min_size_nuclear = (nuclear_size_range[0] * cell_avg_size * 1000/(pixel_size * 2)) ** 2 * math.pi
max_size_nuclear = (nuclear_size_range[1] * cell_avg_size * 1000/(pixel_size * 2)) ** 2 * math.pi

data = pd.DataFrame(columns=['FOV',
                             'z',
                             'label_nuclear',
                             'centroid_nuclear',
                             'area_nuclear',
                             'mean_intensity_MYC_DNAFISH_in_nucleus',
                             'total_intensity_MYC_DNAFISH_in_nucleus'])

# LOAD IMAGE
for fov in range(total_fov):
    im_stack = skio.imread("%s%s_%s_%s.tif" % (master_folder, prefix, sample, fov + 1), plugin="tifffile")
    nuclear_seg_convex = np.zeros(shape=(im_stack.shape[0], im_stack.shape[2], im_stack.shape[3]), dtype=np.uint16)
    for z in range(im_stack.shape[0]):
        print("Start analyzing FOV %s/%s, z %s/%s" % (fov+1, total_fov, z+1, im_stack.shape[0]))
        img_nuclear = im_stack[z][1]
        img_FISH = im_stack[z][0]
        img_nuclear_seg = seg.nuclear_seg1(img_nuclear, local_factor=local_factor_nuclear, min_size=min_size_nuclear,
                                           max_size=max_size_nuclear)
        img_nuclear_seg_convex = seg.obj_to_convex(img_nuclear_seg)
        nuclear_seg_convex[z] = img_nuclear_seg_convex
        FISH_props = regionprops(img_nuclear_seg_convex, img_FISH)
        nuclear_props = regionprops(img_nuclear_seg_convex, img_nuclear)

        for i in range(len(nuclear_props)):
            # nuclear morphology
            nuclear_label = FISH_props[i].label
            nuclear_centroid = FISH_props[i].centroid
            nuclear_area = FISH_props[i].area

            # ecDNA related
            FISH_mean_intensity_nuclear = FISH_props[i].mean_intensity
            FISH_total_intensity_nuclear = nuclear_area * FISH_mean_intensity_nuclear

            data.loc[len(data.index)] = [fov, z, nuclear_label, nuclear_centroid, nuclear_area,
                                         FISH_mean_intensity_nuclear, FISH_total_intensity_nuclear]

    # img.save_3d_image(nuclear_seg_convex, master_folder, 'nuclear_seg_convex_%s_%s' % (sample, fov+1))
    tif.imsave('%snuclear_seg_convex_%s_%s.tif' % (master_folder, sample, fov+1), nuclear_seg_convex)

data.to_csv('%s%s_intensity.txt' % (master_folder, sample), index=False, sep='\t')

print("DONE!")




