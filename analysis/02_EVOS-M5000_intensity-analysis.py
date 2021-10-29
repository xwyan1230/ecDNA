import matplotlib.pyplot as plt
import napari

data_path = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20211022_ColoDM-Cas9_nucleofectionTest/"
img = plt.imread('%s/SF-A1_0129_GFP.tif' % data_path, format=None)

viewer = napari.Viewer()
viewer.add_image(img)
napari.run()