from skimage.measure import regionprops
import shared.segmentation as seg
from skimage.morphology import disk, dilation
import math
import numpy as np
import skimage.io as skio
import tifffile as tif
import napari

# input parameters
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20220325_Natasha_THZ1/"
sample = 'DMSO'
num_z = 22
ch = '02'

# LOAD IMAGE
im_stack = skio.imread("%s%s/%s_RAW_ch%s.tif" % (master_folder, sample, sample, ch), plugin="tifffile")

for i in range(int(im_stack.shape[0]/22)):
    im_stack_fov = im_stack[(i*22):(i*22+22)]
    tif.imsave('%s%s/%s_RAW_ch%s_fov%s.tif' % (master_folder, sample, sample, ch, i), im_stack_fov)

# viewer = napari.Viewer()
# viewer.add_image(im_stack1)
# napari.run()