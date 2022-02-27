import pandas as pd
import matplotlib.pyplot as plt
from skimage.measure import label, regionprops
import shared.segmentation as seg
import shared.image as img
from skimage.morphology import binary_dilation

# input parameters
exp_name = '20220204_CtrlAndJQ1'
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20220204_CtrlAndJQ1/"
prefix = 'ctrl_new_NW_20pos_Processed001'
local_size = 100
total_FOV = 20
sample = 'ctrl'

data = pd.DataFrame(columns=['FOV', 'nuclear', 'nuclear_area', 'nuclear_mean_int', 'ecDNA_number', 'ecDNA_area',
                             'ecDNA_mean_int'])

for fov in range(total_FOV):
    if fov < 10:
        img_nuclear = plt.imread('%s%s/%s_s0%s_RAW_ch00.tif' % (master_folder, prefix, prefix, fov), format=None)
        img_FISH = plt.imread('%s%s/%s_s0%s_RAW_ch01.tif' % (master_folder, prefix, prefix, fov), format=None)
    else:
        img_nuclear = plt.imread('%s%s/%s_s%s_RAW_ch00.tif' % (master_folder, prefix, prefix, fov), format=None)
        img_FISH = plt.imread('%s%s/%s_s%s_RAW_ch01.tif' % (master_folder, prefix, prefix, fov), format=None)

    img_nuclear_seg = seg.nuclear_seg(img_nuclear)
    nuclear_props = regionprops(img_nuclear_seg)

    for i in range(len(nuclear_props)):
        print("Start image analysis FOV %s/%s, nuclear %s/%s..." % (fov+1, total_FOV, i+1, len(nuclear_props)))
        original_centroid = nuclear_props[i].centroid
        position = img.img_local_position(img_nuclear_seg, original_centroid, local_size)
        nuclear_seg = img.img_local_seg(img_nuclear_seg, position, i+1)
        centroid = regionprops(label(nuclear_seg))[0].centroid
        nuclear = img_nuclear[position[0]:position[1], position[2]:position[3]]
        FISH = img_FISH[position[0]:position[1], position[2]:position[3]]
        if nuclear_props[i].area * 1.0 / nuclear_props[i].convex_area > 0.95:
            FISH_seg, _ = seg.find_organelle(FISH, 'na', extreme_val=500, bg_val=seg.get_bg_int([FISH])[0], min_size=0, max_size=10000)
            FISH_seg = binary_dilation(FISH_seg)
            FISH_seg[nuclear_seg == 0] = 0
            ecDNA_props = regionprops(label(FISH_seg), FISH)
            total_props = regionprops(label(nuclear_seg), FISH)
            total_area = total_props[0].area
            total_mean_int = total_props[0].mean_intensity
            ecDNA_n = len(ecDNA_props)
            ecDNA_area = [ecDNA_props[i].area for i in range(len(ecDNA_props))]
            ecDNA_mean_int = [ecDNA_props[i].mean_intensity for i in range(len(ecDNA_props))]
            data.loc[len(data.index)] = [fov, i + 1, total_area, total_mean_int, ecDNA_n, ecDNA_area, ecDNA_mean_int]
        else:
            print("Nuclear does not pass convex filter.")

data.to_csv('%sfeature_analysis_%s.txt' % (master_folder, sample), index=False, sep='\t')

print("END!")