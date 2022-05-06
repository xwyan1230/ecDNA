import matplotlib.pyplot as plt
import numpy as np
import tifffile as tif

# input parameters
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20220407_sp8_DMandHSR/HSR_singleZ/"
prefix = '20220407_DMandHSR_HSR_singleZ'
sample = 'HSR'
total_fov = 6
n_pixel = 3144

im_stack_nuclear = np.zeros(shape=(total_fov, n_pixel, n_pixel), dtype=np.uint16)
im_stack_DNAFISH = np.zeros(shape=(total_fov, n_pixel, n_pixel), dtype=np.uint16)
im_stack_IF = np.zeros(shape=(total_fov, n_pixel, n_pixel), dtype=np.uint16)
for fov in range(total_fov):
    # load images
    im_stack_nuclear[fov] = plt.imread('%s/%s_s%s_ch00.tif' % (master_folder, prefix, fov), format=None)
    im_stack_DNAFISH[fov] = plt.imread('%s/%s_s%s_ch02.tif' % (master_folder, prefix, fov), format=None)
    im_stack_IF[fov] = plt.imread('%s/%s_s%s_ch01.tif' % (master_folder, prefix, fov), format=None)
tif.imsave('%s%s_nuclear.tif' % (master_folder, sample), im_stack_nuclear)
tif.imsave('%s%s_DNAFISH.tif' % (master_folder, sample), im_stack_DNAFISH)
tif.imsave('%s%s_IF.tif' % (master_folder, sample), im_stack_IF)

print("DONE!")
