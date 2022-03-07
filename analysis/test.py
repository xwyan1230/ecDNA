"""import pandas as pd
import matplotlib.pyplot as plt
from skimage.measure import label, regionprops
import math
import shared.segmentation as seg
import shared.dataframe as dat
import shared.image as img
from scipy.io import savemat
import shared.math as mat

# input parameters
exp_name = '20220204_CtrlAndJQ1'
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20220204_CtrlAndJQ1/"
ctrl_prefix = 'ctrl_new_NW_20pos_Processed001'
local_size = 100
rmax = 60

img_ctrl_nuclear = plt.imread('%s%s/%s_s00_RAW_ch00.tif' % (master_folder, ctrl_prefix, ctrl_prefix), format=None)
img_ctrl_FISH = plt.imread('%s%s/%s_s00_RAW_ch01.tif' % (master_folder, ctrl_prefix, ctrl_prefix), format=None)
img_shape = img_ctrl_nuclear.shape

img_ctrl_nuclear_seg = seg.nuclear_seg(img_ctrl_nuclear)
ctrl_nuclear_props = regionprops(img_ctrl_nuclear_seg)

#for i in range(len(ctrl_nuclear_props)):
for i in range(1):
    print("Start image analysis nuclear %s/%s..." % (i+1, len(ctrl_nuclear_props)))
    original_centroid = ctrl_nuclear_props[i].centroid
    position = img.img_local_position(img_ctrl_nuclear_seg, original_centroid, local_size)
    nuclear_seg = img.img_local_seg(img_ctrl_nuclear_seg, position, i+1)
    centroid = regionprops(label(nuclear_seg))[0].centroid
    nuclear = img_ctrl_nuclear[position[0]:position[1], position[2]:position[3]]
    FISH = img_ctrl_FISH[position[0]:position[1], position[2]:position[3]]
    savemat("%sFISH_%s.mat" % (master_folder, i+1), {'FISH': FISH})
    savemat("%smask_%s.mat" % (master_folder, i+1), {'seg': nuclear_seg})
    x = mat.auto_correlation(FISH, nuclear_seg, 60, 1)
    print(x)
    print(len(x))"""

import numpy as np

from skimage.measure import label, regionprops
import shared.segmentation as seg
import shared.image as img
import shared.math as mat
import pandas as pd
import math
import matplotlib.pyplot as plt
import napari

# input parameters
exp_name = '20220204_CtrlAndJQ1'
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20220204_CtrlAndJQ1_zstack/"
sample = 'ctrl'
local_size = 100
rmax = 60
total_FOV = 20
total_z = 16

data = pd.DataFrame(columns=['FOV', 'z', 'nuclear', 'nuclear_centroid', 'nuclear_area', 'nuclear_circ', 'FISH_mean_int',
                             'g_FISH', 'dg_FISH', 'g_nuclear', 'dg_nuclear'])

for fov in range(total_FOV):
    fov_str = '0%s' % fov if fov < 10 else fov
    for z in range(total_z):
        z_str = '0%s' % z if z < 10 else z
        img_nuclear = plt.imread('%s%s/%s_new_NW_20pos_s%s_z%s_RAW_ch00.tif' % (master_folder, sample, sample, fov_str,
                                                                                z_str), format=None)
        img_FISH = plt.imread('%s%s/%s_new_NW_20pos_s%s_z%s_RAW_ch01.tif' % (master_folder, sample, sample, fov_str,
                                                                             z_str), format=None)

        img_nuclear_seg = seg.nuclear_seg(img_nuclear)

        viewer = napari.Viewer()
        viewer.add_image(img_nuclear)
        viewer.add_image(img_nuclear_seg)
        napari.run()