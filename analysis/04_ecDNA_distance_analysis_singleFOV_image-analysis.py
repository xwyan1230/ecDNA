import pandas as pd
import matplotlib.pyplot as plt
from skimage.measure import label, regionprops
import math
import shared.segmentation as seg
import shared.dataframe as dat
import shared.image as img

# input parameters
exp_name = '20220204_CtrlAndJQ1'
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20220204_CtrlAndJQ1/"
ctrl_prefix = 'ctrl_new_NW_20pos_Processed001'
treatment_prefix = 'JQ1_new_NW_20pos_Processed001'
local_size = 100

img_ctrl_nuclear = plt.imread('%s%s/%s_s00_RAW_ch00.tif' % (master_folder, ctrl_prefix, ctrl_prefix), format=None)
img_ctrl_FISH = plt.imread('%s%s/%s_s00_RAW_ch01.tif' % (master_folder, ctrl_prefix, ctrl_prefix), format=None)
img_treatment_nuclear = plt.imread('%s%s/%s_s00_RAW_ch00.tif' % (master_folder, treatment_prefix, treatment_prefix),
                                   format=None)
img_treatment_FISH = plt.imread('%s%s/%s_s00_RAW_ch01.tif' % (master_folder, treatment_prefix, treatment_prefix),
                                format=None)
img_shape = img_ctrl_nuclear.shape

img_ctrl_nuclear_seg = seg.nuclear_seg(img_ctrl_nuclear)
ctrl_nuclear_props = regionprops(img_ctrl_nuclear_seg)

data = pd.DataFrame(columns=['nuclear', 'distance', 'k', 'border', 'intensity_FISH', 'intensity_nucleus',
                             'mean_intensity_FISH', 'mean_intensity_nucleus'])

for i in range(len(ctrl_nuclear_props)):
    print("Start image analysis nuclear %s/%s..." % (i+1, len(ctrl_nuclear_props)))
    original_centroid = ctrl_nuclear_props[i].centroid
    position = img.img_local_position(img_ctrl_nuclear_seg, original_centroid, local_size)
    nuclear_seg = img.img_local_seg(img_ctrl_nuclear_seg, position, i+1)
    centroid = regionprops(label(nuclear_seg))[0].centroid
    nuclear = img_ctrl_nuclear[position[0]:position[1], position[2]:position[3]]
    FISH = img_ctrl_FISH[position[0]:position[1], position[2]:position[3]]
    mean_int_FISH = regionprops(label(nuclear_seg), FISH)[0].mean_intensity
    mean_int_nucleus = regionprops(label(nuclear_seg), nuclear)[0].mean_intensity
    if ctrl_nuclear_props[i].area*1.0/ctrl_nuclear_props[i].convex_area > 0.95:
        nuclear_edge = img.edge_from_seg(nuclear_seg)
        shape = nuclear_seg.shape
        for m in range(shape[0]):
            for n in range(shape[1]):
                if nuclear_seg[m, n] == 1:
                    distance_temp = math.sqrt((m - centroid[0]) ** 2 + (n - centroid[1]) ** 2)
                    k_temp = (n - centroid[1])/(m-centroid[0])
                    border_temp = 1 if nuclear_edge[m, n] == 1 else 0
                    intensity_FISH_temp = FISH[m][n]
                    intensity_nucleus_temp = nuclear[m][n]
                    data.loc[len(data.index)] = [i+1, distance_temp, k_temp, border_temp, intensity_FISH_temp,
                                                 intensity_nucleus_temp, mean_int_FISH, mean_int_nucleus]
    else:
        print("Nuclear does not pass convex filter.")

data.to_csv('%ssummary.txt' % master_folder, index=False, sep='\t')

data_r = pd.DataFrame()
for i in range(int(max(data['nuclear']))):
    print("Start relative r calculation %s/%s..." % (i+1, max(data['nuclear'])))
    r = []
    data_temp = data[data['nuclear'] == i+1]
    if len(data_temp) != 0:
        data_border = data_temp[data_temp['border'] == 1].sort_values(by=['k'])
        for j in range(len(data_temp)):
            pos = dat.find_pos(data_temp['k'].tolist()[j], data_border['k'].tolist())
            if pos == len(data_border):
                pos = pos-1
            r.append(data_border['distance'].tolist()[pos])
    data_temp = data_temp.assign(r=r)
    data_temp['relative_r'] = data_temp['distance'] * 1.0 / data_temp['r']
    data_r = pd.concat([data_r, data_temp])

data_r.to_csv('%ssummary_r.txt' % master_folder, index=False, sep='\t')

print("END!")

"""viewer = napari.Viewer()
viewer.add_image(img_ctrl_nuclear)
viewer.add_image(img_ctrl_FISH)
viewer.add_image(img_ctrl_nuclear_seg)
viewer.add_image(img_ctrl_nuclear_edge)
napari.run()"""
