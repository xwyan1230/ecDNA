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

from skimage.measure import label, regionprops
import shared.segmentation as seg
import shared.image as img
import napari
import skimage.io as skio
import shared.math as mat
import pandas as pd
import math
import matplotlib.pyplot as plt
import numpy as np

# input parameters
exp_name = '20220301_NatashaFile'
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20220301_ecDNA_ctrlAndJQ1_NatashaFile/"
local_size = 250
rmax = 100
sample = 'DMSO_2'

imstack = skio.imread("%s210921_COLODM_washout_mycFISH_DMSO_2.tif" % master_folder, plugin="tifffile")

data = pd.DataFrame(columns=['FOV', 'z', 'nuclear', 'nuclear_centroid', 'nuclear_area', 'nuclear_circ', 'FISH_mean_int',
                             'g_FISH', 'dg_FISH', 'g_nuclear', 'dg_nuclear'])

#for z in range(imstack.shape[0]):
for z in range(1):
    z = 21
    img_nuclear = imstack[z][1]
    img_FISH = imstack[z][0]
    img_nuclear_seg = seg.nuclear_seg1(img_nuclear, local_factor=999, max_size=100000)
    nuclear_props = regionprops(img_nuclear_seg)

    for i in range(len(nuclear_props)):
        print("Start image analysis z %s/%s, nuclear %s/%s..." % (z+1, imstack.shape[0]+1, i+1, len(nuclear_props)))
        original_centroid = nuclear_props[i].centroid

        position = img.img_local_position(img_nuclear_seg, original_centroid, local_size)
        nuclear_seg = img.img_local_seg(img_nuclear_seg, position, i+1)
        nuclear = img_nuclear[position[0]:position[1], position[2]:position[3]]
        centroid = regionprops(label(nuclear_seg))[0].centroid
        FISH = img_FISH[position[0]:position[1], position[2]:position[3]]

        nuclear_convex_local = regionprops(label(nuclear_seg))[0].convex_image
        centroid_convex = regionprops(label(nuclear_convex_local))[0].centroid
        nuclear_convex = img.image_paste(FISH, nuclear_convex_local, [int(centroid[0]-centroid_convex[0]),
                                                                      int(centroid[1]-centroid_convex[1])])
        viewer = napari.Viewer()
        viewer.add_image(nuclear)
        viewer.add_image(nuclear_convex)
        viewer.add_image(FISH)
        napari.run()

        nuclear_area = regionprops(label(nuclear_convex))[0].area
        FISH_mean_int = regionprops(label(nuclear_convex), FISH)[0].mean_intensity
        nuclear_circ = (4 * math.pi * nuclear_area) / (regionprops(label(nuclear_convex))[0].perimeter ** 2)

        _, r, g_FISH, dg_FISH = mat.auto_correlation(FISH, nuclear_seg, rmax)
        _, r, g_nuclear, dg_nuclear = mat.auto_correlation(nuclear, nuclear_seg, rmax)
        data.loc[len(data.index)] = \
            [3, z, i + 1, original_centroid, nuclear_area, nuclear_circ, FISH_mean_int, g_FISH, dg_FISH, g_nuclear,
             dg_nuclear]

        r = np.arange(0, 101, 1)

        plt.subplots(figsize=(6, 4))
        plt.plot(r, g_FISH, color='#FFD700')
        plt.xlabel('r')
        plt.ylabel('g')
        plt.legend(loc=2, bbox_to_anchor=(0.02, 0.99))
        plt.savefig('%s/auto_correlation_%s.pdf' % (master_folder, i+1))
        plt.close()

    data.to_csv('%sauto_correlation_%s.txt' % (master_folder, sample), index=False, sep='\t')


print("END!")

