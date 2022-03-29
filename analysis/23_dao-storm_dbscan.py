import os
import storm_analysis.daostorm_3d.mufit_analysis as mfit
from skimage.measure import label, regionprops
import shared.segmentation as seg
import shared.image as img
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from sklearn.cluster import DBSCAN
from sklearn import metrics


master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20220328_dao-storm_auto-correlationAndDBSCAN/"
local_size = 250
rmax = 100

data = pd.DataFrame(columns=['i', 'g_FISH'])

img_nuclear = plt.imread('%stest_DAPI.tif' % master_folder, format=None)
img_FISH = plt.imread('%stest_FISH.tif' % master_folder, format=None)
img_nuclear_seg = seg.nuclear_seg1(img_nuclear, local_factor=999, max_size=100000)

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

    X, labels_true = img.xy_lst_within_region(nuclear_seg, '%stemp%s.hdf5' % (master_folder, i))

    # Compute DBSCAN
    db = DBSCAN(eps=20, min_samples=15).fit(X)
    core_samples_mask = np.zeros_like(db.labels_, dtype=bool)
    core_samples_mask[db.core_sample_indices_] = True
    labels = db.labels_

    # Number of clusters in labels, ignoring noise if present.
    n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)
    n_noise_ = list(labels).count(-1)

    print('Estimated number of clusters: %d' % n_clusters_)
    print('Estimated number of noise points: %d' % n_noise_)
    print("Homogeneity: %0.3f" % metrics.homogeneity_score(labels_true, labels))
    print("Completeness: %0.3f" % metrics.completeness_score(labels_true, labels))
    print("V-measure: %0.3f" % metrics.v_measure_score(labels_true, labels))
    print("Adjusted Rand Index: %0.3f"
          % metrics.adjusted_rand_score(labels_true, labels))
    print("Adjusted Mutual Information: %0.3f"
          % metrics.adjusted_mutual_info_score(labels_true, labels))
    print("Silhouette Coefficient: %0.3f"
          % metrics.silhouette_score(X, labels))

    # Plot result
    # Black removed and is used for noise instead.
    unique_labels = set(labels)
    colors = [plt.cm.Spectral(each)
              for each in np.linspace(0, 1, len(unique_labels))]
    for k, col in zip(unique_labels, colors):
        if k == -1:
            # Black used for noise.
            col = [0, 0, 0, 1]

        class_member_mask = (labels == k)

        xy = X[class_member_mask & core_samples_mask]
        print(xy)
        plt.plot(xy[:, 1], xy[:, 0], 'o', markerfacecolor=tuple(col),
                 markeredgecolor='k', markersize=10)

        xy = X[class_member_mask & ~core_samples_mask]
        plt.plot(xy[:, 1], xy[:, 0], 'o', markerfacecolor=tuple(col),
                 markeredgecolor='k', markersize=3)

    plt.title('Estimated number of clusters: %d' % n_clusters_)
    plt.savefig('%stemp%s_dbscan.pdf' % (master_folder, i))
    plt.show()





