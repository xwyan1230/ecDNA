import shared.segmentation as seg
from skimage.morphology import disk, dilation
import math
import numpy as np
import skimage.io as skio
import tifffile as tif
import napari

# input parameters
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20220325_Natasha_THZ1/"
sample = 'DMSO'
pixel_size = 102  # nm (sp8 confocal 3144x3144:58.7, mischel 2048x2048:102)
cell_avg_size = 10  # um (Colo)
nuclear_size_range = [0.6, 1.5]  # used to filter nucleus
n_nuclear_convex_dilation = 3

# SET UP PARAMETERS
total_fov = 1
total_z = 22
# segmentation
local_factor_nuclear = 99  # ~99, needs to be odd number, rmax if (rmax % 2 == 1) else rmax+1
min_size_nuclear = (nuclear_size_range[0] * cell_avg_size * 1000/(pixel_size * 2)) ** 2 * math.pi
max_size_nuclear = (nuclear_size_range[1] * cell_avg_size * 1000/(pixel_size * 2)) ** 2 * math.pi

# LOAD IMAGE
im_stack_nuclear = skio.imread("%s%s/%s_RAW_ch00_fov3.tif" % (master_folder, sample, sample), plugin="tifffile")
im_stack_DNAFISH = skio.imread("%s%s/%s_RAW_ch01_fov3.tif" % (master_folder, sample, sample), plugin="tifffile")

# IMAGING ANALYSIS

# Perform nuclear segmentation
for fov in range(total_fov):
    print("Start nuclear segmentation FOV %s/%s" % (fov+1, total_fov))
    # img_nuclear = im_stack_nuclear[fov]
    # img_DNAFISH = im_stack_DNAFISH[fov]
    img_nuclear = im_stack_nuclear
    img_DNAFISH = im_stack_DNAFISH
    nuclear_seg_convex = np.zeros(shape=(img_nuclear.shape[0], img_nuclear.shape[1], img_nuclear.shape[2]),
                                  dtype=np.uint16)
    for z in range(total_z):
        # nuclear segmentation
        img_nuclear_seg = seg.nuclear_seg(img_nuclear[z], local_factor=local_factor_nuclear, min_size=min_size_nuclear,
                                          max_size=max_size_nuclear)
        # DNAFISH_seg, _ = seg.find_organelle(img_DNAFISH[z], 'na', extreme_val=500,
        #                                     bg_val=seg.get_bg_int([img_DNAFISH[z]])[0], min_size=0, max_size=10000)
        # viewer = napari.view_image(img_DNAFISH[z], blending='additive', colormap='green')
        # viewer.add_image(DNAFISH_seg, blending='additive', colormap='red')
        # napari.run()
        img_nuclear_seg_convex = seg.obj_to_convex(img_nuclear_seg)
        img_nuclear_seg_convex = dilation(img_nuclear_seg_convex, disk(n_nuclear_convex_dilation))
        nuclear_seg_convex[z] = img_nuclear_seg_convex

    # tif.imsave('%s%s/%s_seg_fov3.tif' % (master_folder, sample, sample), nuclear_seg_convex)

viewer = napari.view_image(im_stack_nuclear, name='nuclear')
viewer.add_image(nuclear_seg_convex, name='seg')
napari.run()