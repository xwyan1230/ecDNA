import matplotlib.pyplot as plt
from skimage.measure import label, regionprops, regionprops_table
import shared.segmentation as seg
import shared.image as img
from skimage.morphology import binary_dilation, binary_erosion, disk, dilation
import napari
import shared.objects as obj
from scipy import ndimage
import pandas as pd

# input parameters
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20220407_sp8_DMandHSR/DM_singleZ/"
prefix = '20220407_DMandHSR_DM_singleZ'
local_size = 100
# load images
img_nuclear = plt.imread('%s/%s_s0_ch00.tif' % (master_folder, prefix), format=None)
img_FISH = plt.imread('%s/%s_s0_ch02.tif' % (master_folder, prefix), format=None)
# nuclear segmentation
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
img_FISH_seg[img_nuclear_seg_convex == 0] = 0

# filter
img_nuclear_seg_convex = seg.filter_mean_int(img_nuclear_seg_convex, img_FISH_seg, img_FISH, 20000)
# props

viewer = napari.Viewer()
viewer.add_image(img_nuclear_seg_convex, colormap='viridis', blending='additive')
viewer.add_image(img_FISH_seg, blending='additive')
napari.run()



