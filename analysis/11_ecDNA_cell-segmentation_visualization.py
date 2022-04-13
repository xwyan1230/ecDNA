import matplotlib.pyplot as plt
from skimage.measure import label, regionprops, regionprops_table
import shared.segmentation as seg
import shared.image as img
from skimage.morphology import binary_dilation, binary_erosion, disk
import napari
import pandas as pd

# input parameters
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20220407_sp8_DMandHSR/HSR_singleZ/"
prefix = '20220407_DMandHSR_HSR_singleZ'
local_size = 100

img_nuclear = plt.imread('%s/%s_s5_ch00.tif' % (master_folder, prefix), format=None)
img_FISH = plt.imread('%s/%s_s5_ch02.tif' % (master_folder, prefix), format=None)

img_nuclear_seg = seg.nuclear_seg1(img_nuclear, local_factor=199, min_size=9000, max_size=40000)
# 9000, 40000
"""FISH_seg, img_FISH_seg = seg.find_organelle(img_FISH, 'na', extreme_val=18000, bg_val=seg.get_bg_int([img_FISH])[0], min_size=20, max_size=50000)
img_FISH_seg = binary_dilation(img_FISH_seg, disk(8))
for i in range(8):
    img_FISH_seg = binary_erosion(img_FISH_seg)"""
#img_FISH_seg = binary_erosion(img_FISH_seg, disk(8))
"""for i in range(8):
    img_FISH_seg = binary_dilation(img_FISH_seg)
for i in range(8):
    img_FISH_seg = binary_erosion(img_FISH_seg)"""
props = regionprops(img_nuclear_seg)
nuclear_props = pd.DataFrame(regionprops_table(img_nuclear_seg, img_nuclear,
                                               properties=['label', 'area', 'convex_area', 'major_axis_length',
                                                           'minor_axis_length', 'centroid', 'eccentricity',
                                                           'max_intensity', 'min_intensity', 'mean_intensity',
                                                           'solidity']))

nuclear_props.to_csv('%stest.txt' % master_folder, index=False, sep='\t')
viewer = napari.Viewer()
viewer.add_image(img_nuclear)
viewer.add_image(img_FISH)
viewer.add_image(img_nuclear_seg)
viewer.add_points([219, 481], size=10, face_color='red')
napari.run()

"""for i in range(len(nuclear_props)):
    print("Start image analysis nuclear %s/%s..." % (i+1, len(nuclear_props)))
    original_centroid = nuclear_props[i].centroid
    position = img.img_local_position(img_nuclear_seg, original_centroid, local_size)
    nuclear_seg = img.img_local_seg(img_nuclear_seg, position, i+1)
    centroid = regionprops(label(nuclear_seg))[0].centroid
    nuclear = img_nuclear[position[0]:position[1], position[2]:position[3]]
    FISH = img_FISH[position[0]:position[1], position[2]:position[3]]
    if nuclear_props[i].area * 1.0 / nuclear_props[i].convex_area > 0.95:
        FISH_seg, _ = seg.find_organelle(FISH, 'na', extreme_val=500, bg_val=seg.get_bg_int([FISH])[0], min_size=0, max_size=10000)
        FISH_seg = binary_dilation(FISH_seg)
        FISH_seg[nuclear_seg == 0] = 0
        viewer = napari.Viewer()
        viewer.add_image(nuclear)
        viewer.add_image(FISH)
        viewer.add_image(nuclear_seg)
        viewer.add_image(FISH_seg)
        napari.run()
    else:
        print("Nuclear does not pass convex filter.")"""

print("END!")