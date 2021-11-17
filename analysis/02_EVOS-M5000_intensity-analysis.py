import os
import pandas as pd
import matplotlib.pyplot as plt
from skimage.measure import label
from datetime import datetime
import shared.segmentation as seg
import shared.objects as obj
import shared.dataframe as dat
import napari

# input parameters
exp_name = '20211022_EVOS-M5000_ColoDM-Cas9_nucleofectionTest'
master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20211022_ColoDM-Cas9_nucleofectionTest/"
data_source = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20211022_ColoDM-Cas9_nucleofectionTest/EVOS_M5000/"
ch1 = 'GFP'  # channel for intensity measurement
ch2 = 'TRANS'  # channel for segmentation, generally TRANS or BFP

# write log file
f = open("%slog.txt" % master_folder, "a+")
now = datetime.now()
dt_string = now.strftime("%d/%m/%Y %H:%M:%S")
f.write("exp_name: %s\n" % exp_name)
f.write("datetime: %s\n" % dt_string)
f.write("code_running: %s\n" % __file__)
f.write("master_folder: %s\n" % master_folder)
f.write("data_source: %s\n" % data_source)
f.write("Notes:\n")

# get all the files to be analyzed
files = [x for x in os.listdir("%s%s/" % (data_source, ch1))]

# lists for creating data dataframe
data_sample = list()
data_int = list()

f.write("going through data to extract intensity information...\n")
count = 0

for i in files:
    print('# Analyzing %s ... (%d/%d)' % (i[:-(len(ch1)+5)], count+1, len(files)))

    # load images
    img_ch1 = plt.imread('%s%s/%s' % (data_source, ch1, i))[:, :, 1]
    img_ch2 = plt.imread('%s%s/%s%s.tif' % (data_source, ch2, i[:-(len(ch1)+4)], ch2))
    # segment cells and get intensity information
    cell = seg.cell_seg_trans(img_ch2)
    cell_int = dat.list_unwrap(obj.get_intensity(label(cell), [img_ch1]))
    # add data into lists
    data_sample.append(i[:-(len(ch1)+5)])
    data_int.append(cell_int)

    count = count+1

data = pd.DataFrame({'sample': data_sample, 'intensity': data_int})
save_path = master_folder
data.to_csv('%ssummary_intensity.txt' % save_path, index=False, sep='\t')

print("DONE!")
f.write("DONE!\n\n")
f.close()

"""
viewer = napari.Viewer()
viewer.add_image(img_GFP)
viewer.add_image(img_trans)
napari.run()
"""