from skimage.measure import regionprops
import shared.segmentation as seg
from skimage.morphology import disk, dilation
import math
import numpy as np
import skimage.io as skio
import tifffile as tif

# input parameters
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20220420_sp8_DMandBRD4_plate/DM_singleZ_imageJ/"
prefix = 'DM_singleZ_25pos_RAW'
sample = 'DM'
"""master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20220420_sp8_DMandBRD4_plate/BRD4KO_singleZ_imageJ/"
prefix = 'DM-BRD4KO_singleZ_8pos_RAW'
sample = 'BRD4KO'"""
pixel_size = 58.7  # nm (sp8 confocal 3144x3144)
cell_avg_size = 10  # um (Colo)
nuclear_size_range = [0.6, 1.5]  # used to filter nucleus
solidity_threshold_nuclear = 0.9
n_nuclear_convex_dilation = 3
avg_background_MYC_DNAFISH = 8000

# LOAD IMAGE
im_stack_nuclear = skio.imread("%s%s_ch00.tif" % (master_folder, prefix), plugin="tifffile")
im_stack_DNAFISH = skio.imread("%s%s_ch02.tif" % (master_folder, prefix), plugin="tifffile")

# SET UP PARAMETERS
total_fov = im_stack_nuclear.shape[0]
# segmentation
local_factor_nuclear = 99  # ~99, needs to be odd number, rmax if (rmax % 2 == 1) else rmax+1
min_size_nuclear = (nuclear_size_range[0] * cell_avg_size * 1000/(pixel_size * 2)) ** 2 * math.pi
max_size_nuclear = (nuclear_size_range[1] * cell_avg_size * 1000/(pixel_size * 2)) ** 2 * math.pi

# Generate background corrected DNAFISH image
im_stack_DNAFISH_corrected = np.zeros(shape=(im_stack_nuclear.shape[0], im_stack_nuclear.shape[1],
                                             im_stack_nuclear.shape[2]), dtype=np.uint16)
for fov in range(total_fov):
    print("Start generating corrected DNAFISH image FOV %s/%s" % (fov + 1, total_fov))
    img_DNAFISH = im_stack_DNAFISH[fov]
    img_DNAFISH_corrected = img_DNAFISH.copy()
    img_DNAFISH_corrected = np.array(img_DNAFISH_corrected).astype(float)-avg_background_MYC_DNAFISH
    img_DNAFISH_corrected[img_DNAFISH_corrected < 0] = 0
    img_DNAFISH_corrected = img_DNAFISH_corrected.astype(int)
    im_stack_DNAFISH_corrected[fov] = img_DNAFISH_corrected

tif.imsave('%sDNAFISH_corrected_%s_%s.tif' % (master_folder, sample, avg_background_MYC_DNAFISH), im_stack_DNAFISH_corrected)

print("DONE!")