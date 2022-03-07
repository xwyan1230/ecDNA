import matplotlib.pyplot as plt
from skimage.measure import label, regionprops
import shared.segmentation as seg
import shared.image as img
from skimage.morphology import binary_dilation
import napari

# input parameters
exp_name = '20220204_CtrlAndJQ1'
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20220204_CtrlAndJQ1/"
prefix = 'ctrl_new_NW_20pos_Processed001'
local_size = 100

img_nuclear = plt.imread('%s%s/%s_s00_RAW_ch00.tif' % (master_folder, prefix, prefix), format=None)
img_FISH = plt.imread('%s%s/%s_s00_RAW_ch01.tif' % (master_folder, prefix, prefix), format=None)

img_nuclear_seg = seg.nuclear_seg(img_nuclear)
nuclear_props = regionprops(img_nuclear_seg)

viewer = napari.Viewer()
viewer.add_image(img_nuclear)
viewer.add_image(img_FISH)
viewer.add_image(img_nuclear_seg)
napari.run()

for i in range(len(nuclear_props)):
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
        print("Nuclear does not pass convex filter.")

print("END!")