import os
import numpy as np
import napari
import matplotlib.pyplot as plt
from skimage.filters import try_all_threshold, threshold_otsu
from skimage.measure import label, regionprops
from skimage.morphology import binary_dilation, binary_erosion, disk, remove_small_objects
from scipy import ndimage

# input parameters
exp_name = '20211022_EVOS-M5000_ColoDM-Cas9_nucleofectionTest'
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20211022_ColoDM-Cas9_nucleofectionTest/"
data_source = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20211022_ColoDM-Cas9_nucleofectionTest/EVOS_M5000/"
ch1 = 'GFP'  # channel for intensity measurement
ch2 = 'TRANS'  # channel for segmentation, generally TRANS or BFP

files = [x for x in os.listdir("%s%s/" % (data_source, ch1))]

for i in files:
    img_ch1 = plt.imread('%s%s/%s' % (data_source, ch1, i))
    img_ch2 = plt.imread('%s%s/%s%s.tif' % (data_source, ch2, i[:-(len(ch1)+4)], ch2))


#viewer = napari.Viewer()
#viewer.add_image(img_GFP)
#viewer.add_image(img_trans)
#napari.run()