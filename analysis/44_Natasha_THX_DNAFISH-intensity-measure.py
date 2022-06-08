import skimage.io as skio
import pandas as pd
from skimage.measure import regionprops
import napari
import seaborn as sns
import matplotlib.pyplot as plt
import shared.image as img

# input parameters
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20220325_Natasha_THZ1/"
sample = 'DMSO'
pixel_size = 102  # nm (sp8 confocal 3144x3144:58.7, mischel 2048x2048:102)
cell_avg_size = 10  # um (Colo)
n_nuclear_convex_dilation = 3

# SET UP PARAMETERS
total_fov = 1
total_z = 22

# LOAD IMAGE
im_stack_nuclear = skio.imread("%s%s/%s_RAW_ch00_fov3.tif" % (master_folder, sample, sample), plugin="tifffile")
im_stack_DNAFISH = skio.imread("%s%s/%s_RAW_ch01_fov3.tif" % (master_folder, sample, sample), plugin="tifffile")
im_stack_seg = skio.imread("%s%s/%s_seg_fov3.tif" % (master_folder, sample, sample), plugin="tifffile")

data = pd.DataFrame(columns=['FOV',
                             'z',
                             'label_nuclear',
                             'centroid_nuclear',
                             'area_nuclear',
                             'mean_intensity_MYC_DNAFISH_in_nucleus',
                             'total_intensity_MYC_DNAFISH_in_nucleus'])

fov = 0
for z in range(total_z):
    print("Start analyzing FOV %s/%s, z %s/%s" % (fov+1, total_fov, z+1, total_z))
    img_nuclear = im_stack_nuclear[z]
    img_DNAFISH = im_stack_DNAFISH[z]
    img_nuclear_seg = im_stack_seg[z]
    DNAFISH_props = regionprops(img_nuclear_seg, img_DNAFISH)
    nuclear_props = regionprops(img_nuclear_seg, img_nuclear)

    for i in range(len(nuclear_props)):
        # nuclear morphology
        nuclear_label = DNAFISH_props[i].label
        nuclear_centroid = DNAFISH_props[i].centroid
        nuclear_area = DNAFISH_props[i].area

        # ecDNA related
        FISH_mean_intensity_nuclear = DNAFISH_props[i].mean_intensity
        FISH_total_intensity_nuclear = nuclear_area * FISH_mean_intensity_nuclear

        data.loc[len(data.index)] = [fov, z, nuclear_label, nuclear_centroid, nuclear_area,
                                     FISH_mean_intensity_nuclear, FISH_total_intensity_nuclear]

print("Analyzing FOV %s/%s" % (fov + 1, total_fov))
nuclear_centroid_searching_range = 25  # pixel
data_z = data
data_z['centroid_x'] = [data_z['centroid_nuclear'][i][0] for i in range(len(data_z))]
data_z['centroid_y'] = [data_z['centroid_nuclear'][i][1] for i in range(len(data_z))]
data_z_fov = data_z[data_z['FOV'] == fov]
z_analyze = int(total_z / 2)
data_z_fov_analyze = data_z_fov[data_z_fov['z'] == z_analyze]
n_nuclear = len(data_z_fov[data_z_fov['z'] == z_analyze])
for i in range(n_nuclear):
    print("Analyzing nucleus %s/%s" % (i+1, n_nuclear))
    data_temp = data_z_fov[(data_z_fov['centroid_x'] > data_z_fov_analyze['centroid_x'].tolist()[i] -
                            nuclear_centroid_searching_range)
                           & (data_z_fov['centroid_x'] < data_z_fov_analyze['centroid_x'].tolist()[i] +
                              nuclear_centroid_searching_range)
                           & (data_z_fov['centroid_y'] > data_z_fov_analyze['centroid_y'].tolist()[i] -
                              nuclear_centroid_searching_range)
                           & (data_z_fov['centroid_y'] < data_z_fov_analyze['centroid_y'].tolist()[i] +
                              nuclear_centroid_searching_range)]
    # print(data_temp)
    for j in range(len(data_temp)):
        z = data_temp['z'].tolist()[j]
        label_nuclear = data_temp['label_nuclear'].tolist()[j]
        original_centroid_nuclear = data_temp['centroid_nuclear'].tolist()[j]

        img_nuclear_seg_convex = im_stack_seg[z]
        img_nuclear = im_stack_nuclear[z]
        img_DNAFISH = im_stack_DNAFISH[z]

        local_size = 100
        position = img.img_local_position(img_nuclear_seg_convex, original_centroid_nuclear, local_size)
        local_nuclear_seg_convex = img.img_local_seg(img_nuclear_seg_convex, position, label_nuclear)
        local_nuclear = img_nuclear.copy()
        local_nuclear = local_nuclear[position[0]:position[1], position[2]:position[3]]
        local_DNAFISH = img_DNAFISH.copy()
        local_DNAFISH = local_DNAFISH[position[0]:position[1], position[2]:position[3]]
        local_DNAFISH[local_nuclear_seg_convex == 0] = 0
        # temp = DNAFISH_props[i].intensity_image
        viewer = napari.view_image(local_DNAFISH)
        napari.run()
        ax = plt.hist(local_DNAFISH.ravel(), bins=256)
        plt.ylim([0, 1000])
        plt.show()

# data.to_csv('%s%s/%s_fov3_intensity.txt' % (master_folder, sample, sample), index=False, sep='\t')

print("DONE!")


