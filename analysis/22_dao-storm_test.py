import os
import skimage.io as skio
import storm_analysis.jupyter_examples.overlay_image as overlay_image
import storm_analysis.daostorm_3d.mufit_analysis as mfit
from skimage.measure import label, regionprops
import shared.segmentation as seg
import shared.image as img
import shared.math as mat
import matplotlib.pyplot as plt
import napari
import pandas as pd
import math
import imageio
import numpy as np
from sklearn.cluster import DBSCAN
from sklearn import metrics
from sklearn.datasets import make_blobs
from sklearn.preprocessing import StandardScaler

master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20220328_dao-storm_auto-correlationAndDBSCAN/"
local_size = 250
rmax = 100

data = pd.DataFrame(columns=['i', 'g_FISH'])

img_nuclear = plt.imread('%stest_DAPI.tif' % master_folder, format=None)
img_FISH = plt.imread('%stest_FISH.tif' % master_folder, format=None)
img_nuclear_seg = seg.nuclear_seg1(img_nuclear, local_factor=999, max_size=100000)

"""viewer = napari.Viewer()
viewer.add_image(img_nuclear)
viewer.add_image(img_nuclear_seg)
napari.run()"""

nuclear_props = regionprops(img_nuclear_seg)

for i in range(len(nuclear_props)):
    print(i)
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
    # plt.imsave('%stemp_%s.tiff' % (master_folder, i), FISH)
    # imageio.imwrite('%stemp%s.tiff' % (master_folder, i), FISH)

    if os.path.exists('%stemp%s.hdf5' % (master_folder, i)):
        os.remove('%stemp%s.hdf5' % (master_folder, i))
    mfit.analyze('%stemp%s.tiff' % (master_folder, i), '%stemp%s.hdf5' % (master_folder, i),
                 "%stesting.xml" % master_folder)
    # img.overlayImage('%stemp%s.tiff' % (master_folder, i), '%stemp%s.hdf5' % (master_folder, i), 0,
    #                  master_folder, 'temp%s' % i)

    FISH_daostorm = img.scatter_dot_from_hdf5(np.zeros_like(FISH), nuclear_seg, '%stemp%s.hdf5' % (master_folder, i))
    # img.display_dot('%stemp%s.tiff' % (master_folder, i), FISH_daostorm)

    _, r, g_FISH, dg_FISH = mat.auto_correlation(FISH_daostorm, nuclear_seg, rmax)

    data.loc[len(data.index)] = [i, g_FISH]

    plt.subplots(figsize=(6, 4))
    plt.plot(r, g_FISH, color='#40E0D0')
    plt.axhline(y=1, color='#FF4500', linestyle='--')
    plt.xlabel('r')
    plt.ylabel('g')
    plt.savefig('%stemp%s_auto_correlation.pdf' % (master_folder, i))
    plt.ylim([-0.5, 20.5])
    plt.savefig('%stemp%s_auto_correlation_part.pdf' % (master_folder, i))
    plt.close()

    """viewer = napari.Viewer()
    viewer.add_image(FISH)
    viewer.add_image(FISH_daostorm)
    napari.run()"""

plt.subplots(figsize=(6, 4))
for i in range(len(data)):
    plt.plot(r, data['g_FISH'][i], label='%s' % i)
plt.axhline(y=1, color='#FF4500', linestyle='--')
plt.xlabel('r')
plt.ylabel('g')
plt.legend()
plt.savefig('%sauto_correlation.pdf' % master_folder)
plt.ylim([-0.5, 20.5])
plt.savefig('%sauto_correlation_part.pdf' % master_folder)
plt.close()



