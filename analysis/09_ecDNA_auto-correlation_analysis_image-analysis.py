import pandas as pd
import matplotlib.pyplot as plt
from skimage.measure import regionprops
import shared.segmentation as seg
import shared.image as img
import shared.math as mat

# input parameters
exp_name = '20220204_CtrlAndJQ1'
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20220204_CtrlAndJQ1/"
prefix = 'JQ1_new_NW_20pos_Processed001'
local_size = 100
rmax = 60
total_FOV = 20
sample = 'JQ1'

data = pd.DataFrame(columns=['FOV', 'nuclear', 'r', 'g_FISH', 'dg_FISH', 'g_nuclear', 'dg_nuclear'])

for fov in range(total_FOV):
    if fov < 10:
        img_ctrl_nuclear = plt.imread('%s%s/%s_s0%s_RAW_ch00.tif' % (master_folder, prefix, prefix, fov), format=None)
        img_ctrl_FISH = plt.imread('%s%s/%s_s0%s_RAW_ch01.tif' % (master_folder, prefix, prefix, fov), format=None)
    else:
        img_ctrl_nuclear = plt.imread('%s%s/%s_s%s_RAW_ch00.tif' % (master_folder, prefix, prefix, fov), format=None)
        img_ctrl_FISH = plt.imread('%s%s/%s_s%s_RAW_ch01.tif' % (master_folder, prefix, prefix, fov), format=None)

    img_ctrl_nuclear_seg = seg.nuclear_seg(img_ctrl_nuclear)
    ctrl_nuclear_props = regionprops(img_ctrl_nuclear_seg)

    for i in range(len(ctrl_nuclear_props)):
        print("Start image analysis FOV %s/%s, nuclear %s/%s..." % (fov+1, total_FOV, i+1, len(ctrl_nuclear_props)))
        original_centroid = ctrl_nuclear_props[i].centroid
        position = img.img_local_position(img_ctrl_nuclear_seg, original_centroid, local_size)
        nuclear_seg = img.img_local_seg(img_ctrl_nuclear_seg, position, i+1)
        nuclear = img_ctrl_nuclear[position[0]:position[1], position[2]:position[3]]
        FISH = img_ctrl_FISH[position[0]:position[1], position[2]:position[3]]
        if ctrl_nuclear_props[i].area * 1.0 / ctrl_nuclear_props[i].convex_area > 0.95:
            _, r, g_FISH, dg_FISH = mat.auto_correlation(FISH, nuclear_seg, rmax)
            _, r, g_nuclear, dg_nuclear = mat.auto_correlation(nuclear, nuclear_seg, rmax)
            data.loc[len(data.index)] = [fov, i + 1, r, g_FISH, dg_FISH, g_nuclear, dg_nuclear]
        else:
            print("Nuclear does not pass convex filter.")

data.to_csv('%sauto_correlation_%s.txt' % (master_folder, sample), index=False, sep='\t')

print("END!")



