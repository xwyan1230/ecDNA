from skimage.measure import label, regionprops
import shared.segmentation as seg
import shared.image as img
import napari
import skimage.io as skio
import shared.math as mat
import pandas as pd
import math
import imageio
import os
import storm_analysis.daostorm_3d.mufit_analysis as mfit
import numpy as np
import matplotlib.pyplot as plt

# input parameters
exp_name = '20220301_NatashaFile'
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20220301_ecDNA_ctrlAndJQ1_NatashaFile/"
local_size = 250
rmax = 100
sample = 'DMSO_3'

imstack = skio.imread("%s210921_COLODM_washout_mycFISH_%s.tif" % (master_folder, sample), plugin="tifffile")

data = pd.DataFrame(columns=['FOV', 'z', 'nuclear', 'nuclear_centroid', 'nuclear_area', 'nuclear_circ', 'FISH_mean_int',
                             'g_FISH', 'dg_FISH', 'g_nuclear', 'dg_nuclear'])

for z in range(imstack.shape[0]):
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

        nuclear_area = regionprops(label(nuclear_convex))[0].area
        FISH_mean_int = regionprops(label(nuclear_convex), FISH)[0].mean_intensity
        nuclear_circ = (4 * math.pi * nuclear_area) / (regionprops(label(nuclear_convex))[0].perimeter ** 2)

        plt.imsave('%sFISH-color_z%s_i%s.tiff' % (master_folder, z, i), FISH)
        imageio.imwrite('%sFISH_z%s_i%s.tiff' % (master_folder, z, i), FISH)
        plt.imsave('%snuclear-color_z%s_i%s.tiff' % (master_folder, z, i), nuclear)
        imageio.imwrite('%snuclear_z%s_i%s.tiff' % (master_folder, z, i), nuclear)

        if os.path.exists('%sFISH_z%s_i%s.hdf5' % (master_folder, z, i)):
            os.remove('%sFISH_z%s_i%s.hdf5' % (master_folder, z, i))
        mfit.analyze('%sFISH_z%s_i%s.tiff' % (master_folder, z, i), '%sFISH_z%s_i%s.hdf5' % (master_folder, z, i),
                     "%stesting.xml" % master_folder)
        img.overlayImage('%sFISH_z%s_i%s.tiff' % (master_folder, z, i), '%sFISH_z%s_i%s.hdf5' % (master_folder, z, i),
                         0, master_folder, 'FISH_z%s_i%s' % (z, i))
        if os.path.exists('%snuclear_z%s_i%s.hdf5' % (master_folder, z, i)):
            os.remove('%snuclear_z%s_i%s.hdf5' % (master_folder, z, i))
        mfit.analyze('%snuclear_z%s_i%s.tiff' % (master_folder, z, i), '%snuclear_z%s_i%s.hdf5' % (master_folder, z, i),
                     "%stesting.xml" % master_folder)
        img.overlayImage('%snuclear_z%s_i%s.tiff' % (master_folder, z, i), '%snuclear_z%s_i%s.hdf5' % (master_folder, z, i),
                         0, master_folder, 'nuclear_z%s_i%s' % (z, i))

        FISH_daostorm = img.scatter_dot_from_hdf5(np.zeros_like(FISH), nuclear_seg,
                                                  '%sFISH_z%s_i%s.hdf5' % (master_folder, z, i))
        nuclear_daostorm = img.scatter_dot_from_hdf5(np.zeros_like(nuclear), nuclear_seg,
                                                     '%snuclear_z%s_i%s.hdf5' % (master_folder, z, i))

        if np.sum(FISH_daostorm) != 0:
            _, r, g_FISH, dg_FISH = mat.auto_correlation(FISH_daostorm, nuclear_seg, rmax)
        else:
            g_FISH = [np.nan]
            dg_FISH = [np.nan]
        if np.sum(nuclear_daostorm) != 0:
            _, r, g_nuclear, dg_nuclear = mat.auto_correlation(nuclear_daostorm, nuclear_seg, rmax)
        else:
            g_nuclear = [np.nan]
            dg_nuclear = [np.nan]

        data.loc[len(data.index)] = \
            [3, z, i, original_centroid, nuclear_area, nuclear_circ, FISH_mean_int, g_FISH, dg_FISH, g_nuclear, dg_nuclear]

data.to_csv('%sauto_correlation_%s.txt' % (master_folder, sample), index=False, sep='\t')

print("END!")



"""viewer = napari.Viewer()
        viewer.add_image(nuclear)
        viewer.add_image(nuclear_convex, opacity=0.5, colormap='red')
        napari.run()"""



