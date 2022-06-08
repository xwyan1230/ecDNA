import skimage.io as skio
import tifffile as tif
import numpy as np
import napari
import skimage as ski
import matplotlib.pyplot as plt
from skimage.filters import sobel, threshold_otsu
from scipy import ndimage as ndi
from skimage.segmentation import watershed

master_folder = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20220325_Natasha_THZ1/"
sample = 'DMSO'
fov = 3

nuclei = skio.imread("%s%s/%s_RAW_ch00_fov%s.tif" % (master_folder, sample, sample, fov), plugin="tifffile")
print("shape: {}".format(nuclei.shape))
print("dtype: {}".format(nuclei.dtype))
print("range: ({}, {})".format(np.min(nuclei), np.max(nuclei)))

# The microscope reports the following spacing (in µm)
original_spacing = np.array([0.50, 0.065, 0.065])
# We downsampled each slice 4x to make the data smaller
rescaled_spacing = original_spacing * [1, 4, 4]
# Normalize the spacing so that pixels are a distance of 1 apart
spacing = rescaled_spacing / rescaled_spacing[2]


# Helper function for plotting histograms.
def plot_hist(ax, data, title=None):
    ax.hist(data.ravel(), bins=256)
    ax.ticklabel_format(axis="y", style="scientific", scilimits=(0, 0))

    if title:
        ax.set_title(title)

"""equalized = ski.exposure.equalize_hist(nuclei)

fig, ((a, b), (c, d)) = plt.subplots(nrows=2, ncols=2)

plot_hist(a, nuclei, title="Original")
plot_hist(b, equalized, title="Histogram equalization")

cdf, bins = ski.exposure.cumulative_distribution(nuclei.ravel())
c.plot(bins, cdf, "r")
c.set_title("Original CDF")

cdf, bins = ski.exposure.cumulative_distribution(equalized.ravel())
d.plot(bins, cdf, "r")
d.set_title("Histogram equalization CDF");

fig.tight_layout()
# plt.show()

vmin, vmax = np.quantile(nuclei, q=(0.005, 0.995))

stretched = ski.exposure.rescale_intensity(
    nuclei,
    in_range=(vmin, vmax),
    out_range=np.float32
)"""

edges = sobel(nuclei)

viewer = napari.view_image(nuclei, blending='additive', colormap='green', name='nuclei')
#viewer.add_image(edges, blending='additive', colormap='magenta', name='edges')

thresholded = nuclei > threshold_otsu(nuclei)
# viewer.add_image(thresholded, name='thresholded', opacity=0.3)

transformed = ndi.distance_transform_edt(thresholded, sampling=spacing)

maxima = ski.morphology.local_maxima(transformed)
viewer.add_points(np.transpose(np.nonzero(maxima)), name='bad points')
viewer.camera.angles

viewer.layers['bad points'].visible = False
points = viewer.add_points(name='interactive points', ndim=3)
points.mode = 'add'
napari.run()

marker_locations = points.data

markers = np.zeros(nuclei.shape, dtype=np.uint32)
marker_indices = tuple(np.round(marker_locations).astype(int).T)
markers[marker_indices] = np.arange(len(marker_locations)) + 1
markers_big = ski.morphology.dilation(markers, ski.morphology.ball(5))

segmented = watershed(
    edges,
    markers_big,
    mask=thresholded
)

viewer = napari.view_image(nuclei, contrast_limits=[0, 3000], scale=spacing)
viewer.add_labels(segmented, name='segmented')
napari.run()