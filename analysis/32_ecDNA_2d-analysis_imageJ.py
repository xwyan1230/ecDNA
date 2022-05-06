from skimage.measure import regionprops
import shared.segmentation as seg
from skimage.morphology import disk, dilation
import math
import numpy as np
import skimage.io as skio
import tifffile as tif

# input parameters
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20220407_sp8_DMandHSR/DM_singleZ/"
prefix = '20220407_DMandHSR_DM_singleZ'
sample = 'DM'
pixel_size = 58.7  # nm (sp8 confocal 3144x3144)
cell_avg_size = 10  # um (Colo)
nuclear_size_range = [0.6, 1.5]  # used to filter nucleus
solidity_threshold_nuclear = 0.9
n_nuclear_convex_dilation = 3

# LOAD IMAGE
im_stack_nuclear = skio.imread("%s%s_nuclear.tif" % (master_folder, sample), plugin="tifffile")
im_stack_DNAFISH = skio.imread("%s%s_DNAFISH.tif" % (master_folder, sample), plugin="tifffile")

# SET UP PARAMETERS
total_fov = im_stack_nuclear.shape[0]
# segmentation
local_factor_nuclear = 99  # ~99, needs to be odd number, rmax if (rmax % 2 == 1) else rmax+1
min_size_nuclear = (nuclear_size_range[0] * cell_avg_size * 1000/(pixel_size * 2)) ** 2 * math.pi
max_size_nuclear = (nuclear_size_range[1] * cell_avg_size * 1000/(pixel_size * 2)) ** 2 * math.pi

# IMAGING ANALYSIS

# Perform nuclear segmentation
nuclear_seg_convex = np.zeros(shape=(im_stack_nuclear.shape[0], im_stack_nuclear.shape[1], im_stack_nuclear.shape[2]),
                              dtype=np.uint16)
mean_intensity_background_MYC_DNAFISH = []
for fov in range(total_fov):
    print("Start nuclear segmentation FOV %s/%s" % (fov+1, total_fov))
    img_nuclear = im_stack_nuclear[fov]
    img_DNAFISH = im_stack_DNAFISH[fov]

    # nuclear segmentation
    img_nuclear_seg = seg.nuclear_seg1(img_nuclear, local_factor=local_factor_nuclear, min_size=min_size_nuclear,
                                       max_size=max_size_nuclear)
    img_nuclear_seg = seg.filter_solidity(img_nuclear_seg, threshold=solidity_threshold_nuclear)
    img_nuclear_seg_convex = seg.obj_to_convex(img_nuclear_seg)
    img_nuclear_seg_convex = dilation(img_nuclear_seg_convex, disk(n_nuclear_convex_dilation))
    nuclear_seg_convex[fov] = img_nuclear_seg_convex

    # measure DNAFISH background
    """DNAFISH_props_temp = regionprops(img_nuclear_seg_convex, img_DNAFISH)
    mean_intensity_background_MYC_DNAFISH = \
        mean_intensity_background_MYC_DNAFISH \
        + [DNAFISH_props_temp[j].mean_intensity for j in range(len(DNAFISH_props_temp))]
avg_background_MYC_DNAFISH = sum(mean_intensity_background_MYC_DNAFISH)*1.0/len(mean_intensity_background_MYC_DNAFISH)
print(avg_background_MYC_DNAFISH)"""
tif.imsave('%s%s_nuclear_seg_convex.tif' % (master_folder, sample), nuclear_seg_convex)

"""# Generate background corrected DNAFISH image
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

tif.imsave('%sDNAFISH_corrected_%s.tif' % (master_folder, sample), im_stack_DNAFISH_corrected)"""

print("DONE!")




