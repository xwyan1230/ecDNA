import scipy.io
import shared.math as mat
import matplotlib.pyplot as plt
import numpy as np
import shared.image as ima
import napari

data_path = "/Users/xwyan/Dropbox/LAB/ChangLab/Projects/Data/20220315_auto-correlation_test/"

# step 1 test if homologous image gives out 1

"""im_FISH = scipy.io.loadmat('%stest-FISH_N250-R50-homo5.mat' % data_path)['signal']
im_seg = scipy.io.loadmat('%stest-mask_N250-R50.mat' % data_path)['out']
im_FISH1 = im_seg * 50

_, r, g, _ = mat.auto_correlation(im_FISH, im_seg, 100)
_, r1, g1, _ = mat.auto_correlation(im_FISH1, im_seg, 100)

plt.subplots(figsize=(6, 4))
plt.plot(r[:60], g[:60], color='#FFD700', label='avgInt=5')
plt.plot(r1[:60], g1[:60], alpha=0.4, color='#40E0D0', label='avgInt=50')
plt.xlabel('r')
plt.ylabel('g')
plt.legend(loc=2, bbox_to_anchor=(0.02, 0.99))
plt.savefig('%s/auto_correlation_test.pdf' % data_path)
plt.close()"""

# step 2 test different shapes

"""im_FISH = np.ones((250, 250)) * 5
im_test = np.zeros((250, 250))
im_seg = scipy.io.loadmat('%stest-mask_N250-R50.mat' % data_path)['out']
im_seg1 = ima.logical_ellipse(im_test, 125, 125, 60, 40)
im_seg2 = ima.logical_ellipse(im_test, 125, 125, 70, 30)

_, r, g, _ = mat.auto_correlation(im_FISH, im_seg, 100)
_, r1, g1, _ = mat.auto_correlation(im_FISH, im_seg1, 100)
_, r2, g2, _ = mat.auto_correlation(im_FISH, im_seg2, 100)
print(g)
print(g1)
print(g2)

plt.subplots(figsize=(6, 4))
plt.plot(r[:60], g[:60], color='#FFD700', label='circle R50')
plt.plot(r1[:60], g1[:60], alpha=0.4, color='#40E0D0', label='ellipse R60-40')
plt.plot(r2[:60], g2[:60], alpha=0.4, color='#DA70D6', label='ellipse R70-30')
plt.ylim([0.95, 1.05])
plt.xlabel('r')
plt.ylabel('g')
plt.legend(loc=2, bbox_to_anchor=(0.02, 0.99))
plt.savefig('%s/auto_correlation_test.pdf' % data_path)
plt.close()"""

# step 3: single cluster, different intensity
"""im_test = np.zeros((250, 250))
im_seg = ima.logical_ellipse(im_test, 125, 125, 50, 50)
im_FISH = ima.logical_ellipse(im_test, 125, 125, 10, 10)
im_FISH1 = ima.logical_ellipse(im_test, 125, 125, 10, 10, 5)
im_FISH2 = ima.logical_ellipse(im_test, 125, 125, 10, 10, 50)

_, r, g, _ = mat.auto_correlation(im_FISH, im_seg, 100)
_, r1, g1, _ = mat.auto_correlation(im_FISH1, im_seg, 100)
_, r2, g2, _ = mat.auto_correlation(im_FISH2, im_seg, 100)
print(g)
print(g1)
print(g2)

plt.subplots(figsize=(6, 4))
plt.plot(r[:60], g[:60], color='#FFD700', label='avgInt=1')
plt.plot(r1[:60], g1[:60], alpha=0.4, color='#40E0D0', label='avgInt=5')
plt.plot(r2[:60], g2[:60], alpha=0.4, color='#DA70D6', label='avgInt=50')
#plt.ylim([0.95, 1.05])
plt.xlabel('r')
plt.ylabel('g')
plt.legend(loc=2, bbox_to_anchor=(0.02, 0.99))
plt.savefig('%s/auto_correlation_test.pdf' % data_path)
plt.close()"""

# step 4: single cluster, different background
"""im_test = np.zeros((250, 250))
im_seg = ima.logical_ellipse(im_test, 125, 125, 50, 50)
im_FISH = ima.logical_ellipse(im_test, 125, 125, 10, 10, 50)
im_test1 = np.ones((250, 250))
im_FISH1 = ima.logical_ellipse(im_test1*5, 125, 125, 10, 10, 50)
im_FISH2 = ima.logical_ellipse(im_test1*25, 125, 125, 10, 10, 50)

_, r, g, _ = mat.auto_correlation(im_FISH, im_seg, 100)
_, r1, g1, _ = mat.auto_correlation(im_FISH1, im_seg, 100)
_, r2, g2, _ = mat.auto_correlation(im_FISH2, im_seg, 100)
print(g)
print(g1)
print(g2)

plt.subplots(figsize=(6, 4))
plt.plot(r[:60], g[:60], color='#FFD700', label='BgInt=0')
plt.plot(r1[:60], g1[:60], alpha=0.4, color='#40E0D0', label='BgInt=5')
plt.plot(r2[:60], g2[:60], alpha=0.4, color='#DA70D6', label='BgInt=25')
#plt.ylim([-0.5, 5])
plt.xlabel('r')
plt.ylabel('g')
plt.legend(loc=2, bbox_to_anchor=(0.02, 0.99))
plt.savefig('%s/auto_correlation_test.pdf' % data_path)
plt.close()
"""
# step 5: single cluster, same intensity ratio Bg/Int
"""im_test = np.zeros((250, 250))
im_seg = ima.logical_ellipse(im_test, 125, 125, 50, 50)
im_test1 = np.ones((250, 250))
im_FISH = ima.logical_ellipse(im_test1, 125, 125, 10, 10, 10)
im_FISH1 = ima.logical_ellipse(im_test1*2, 125, 125, 10, 10, 20)
im_FISH2 = ima.logical_ellipse(im_test1*5, 125, 125, 10, 10, 50)

_, r, g, _ = mat.auto_correlation(im_FISH, im_seg, 100)
_, r1, g1, _ = mat.auto_correlation(im_FISH1, im_seg, 100)
_, r2, g2, _ = mat.auto_correlation(im_FISH2, im_seg, 100)
print(g)
print(g1)
print(g2)

plt.subplots(figsize=(6, 4))
plt.plot(r[:60], g[:60], color='#FFD700', label='BgInt=1, AvgInt=10')
plt.plot(r1[:60], g1[:60], alpha=0.4, color='#40E0D0', label='BgInt=2, AvgInt=20')
plt.plot(r2[:60], g2[:60], alpha=0.4, color='#DA70D6', label='BgInt=5, AvgInt=50')
#plt.ylim([-0.5, 5])
plt.xlabel('r')
plt.ylabel('g')
plt.legend(loc=2, bbox_to_anchor=(0.02, 0.99))
plt.savefig('%s/auto_correlation_test.pdf' % data_path)
plt.close()"""

# step 6: 2,3,4 cluster, same size, binary image
"""im_test = np.zeros((250, 250))
im_seg = ima.logical_ellipse(im_test, 125, 125, 50, 50)
im_FISH = ima.logical_ellipse(im_test, 125, 125, 5, 5)
im_FISH1 = ima.logical_ellipse(im_FISH, 90, 100, 5, 5)
im_FISH2 = ima.logical_ellipse(im_FISH1, 130, 150, 5, 5)

_, r, g, _ = mat.auto_correlation(im_FISH, im_seg, 100)
_, r1, g1, _ = mat.auto_correlation(im_FISH1, im_seg, 100)
_, r2, g2, _ = mat.auto_correlation(im_FISH2, im_seg, 100)
print(g)
print(g1)
print(g2)

plt.subplots(figsize=(6, 4))
plt.plot(r[:60], g[:60], color='#FFD700', label='1 cluster')
plt.plot(r1[:60], g1[:60], alpha=0.4, color='#40E0D0', label='2 clusters')
plt.plot(r2[:60], g2[:60], alpha=0.4, color='#DA70D6', label='3 clusters')
#plt.ylim([-0.5, 5])
plt.xlabel('r')
plt.ylabel('g')
plt.legend(loc=2, bbox_to_anchor=(0.02, 0.99))
plt.savefig('%s/auto_correlation_test.pdf' % data_path)
plt.close()"""

# step 7 & 8: DNA FISH object mimic (random)
"""im_test = np.zeros((250, 250))
im_seg = ima.logical_ellipse(im_test, 125, 125, 50, 50)
im_FISH = ima.logical_dot_sample(im_test, im_seg, 100)
im_FISH1 = ima.logical_dot_sample(im_test, im_seg, 100)
im_FISH2 = ima.logical_dot_sample(im_test, im_seg, 100)

_, r, g, _ = mat.auto_correlation(im_FISH, im_seg, 100)
_, r1, g1, _ = mat.auto_correlation(im_FISH1, im_seg, 100)
_, r2, g2, _ = mat.auto_correlation(im_FISH2, im_seg, 100)
print(g)
print(g1)
print(g2)

plt.subplots(figsize=(6, 4))
plt.plot(r[:60], g[:60], color='#FFD700', label='100 dots rep1')
plt.plot(r1[:60], g1[:60], alpha=0.4, color='#40E0D0', label='100 dots rep2')
plt.plot(r2[:60], g2[:60], alpha=0.4, color='#DA70D6', label='100 dots rep3')
#plt.ylim([-0.5, 5])
plt.xlabel('r')
plt.ylabel('g')
plt.legend(loc=2, bbox_to_anchor=(0.02, 0.99))
plt.savefig('%s/auto_correlation_test.pdf' % data_path)
plt.close()"""

# step 9: DNA FISH object mimic (cluster)
"""im_test = np.zeros((250, 250))
im_seg = ima.logical_ellipse(im_test, 125, 125, 50, 50)
im_FISH = ima.logical_dot_sample(im_test, im_seg, 100)
im_test1 = ima.logical_ellipse(im_test, 125, 125, 35, 35)
im_FISH1 = ima.logical_dot_sample(im_test, im_test1, 100)
im_test2 = ima.logical_ellipse(im_test, 125, 125, 20, 20)
im_FISH2 = ima.logical_dot_sample(im_test, im_test2, 100)

_, r, g, _ = mat.auto_correlation(im_FISH, im_seg, 100)
_, r1, g1, _ = mat.auto_correlation(im_FISH1, im_seg, 100)
_, r2, g2, _ = mat.auto_correlation(im_FISH2, im_seg, 100)
print(g)
print(g1)
print(g2)"""

# step 10: DNA FISH object mimic (cluster, different number)
"""im_test = np.zeros((250, 250))
im_seg = ima.logical_ellipse(im_test, 125, 125, 50, 50)
im_test1 = ima.logical_ellipse(im_test, 125, 125, 20, 20)
im_FISH = ima.logical_dot_sample(im_test, im_test1, 100)
im_FISH1 = ima.logical_dot_sample(im_test, im_test1, 50)
im_FISH2 = ima.logical_dot_sample(im_test, im_test1, 25)

_, r, g, _ = mat.auto_correlation(im_FISH, im_seg, 100)
_, r1, g1, _ = mat.auto_correlation(im_FISH1, im_seg, 100)
_, r2, g2, _ = mat.auto_correlation(im_FISH2, im_seg, 100)
print(g)
print(g1)
print(g2)

plt.subplots(figsize=(6, 4))
plt.plot(r[:60], g[:60], color='#FFD700', label='100 dots R20')
plt.plot(r1[:60], g1[:60], alpha=0.4, color='#40E0D0', label='50 dots R20')
plt.plot(r2[:60], g2[:60], alpha=0.4, color='#DA70D6', label='25 dots R20')
#plt.ylim([-0.5, 5])
plt.xlabel('r')
plt.ylabel('g')
plt.legend(loc=2, bbox_to_anchor=(0.02, 0.99))
plt.savefig('%s/auto_correlation_test.pdf' % data_path)
plt.ylim([-0.5, 10.5])
plt.savefig('%s/auto_correlation_test_part.pdf' % data_path)
plt.close()"""

# step 11: 2,3,4 cluster, same size, binary image, random sampling
"""im_test = np.zeros((250, 250))
im_seg = ima.logical_ellipse(im_test, 125, 125, 50, 50)
im_FISH_template = ima.logical_ellipse(im_test, 125, 125, 8, 8)
im_FISH = ima.logical_dot_sample(im_test, im_FISH_template, 100)
im_FISH1_template = ima.logical_ellipse(im_test, 125, 125, 6, 5)
im_FISH1_template = ima.logical_ellipse(im_FISH1_template, 90, 100, 6, 5)
im_FISH1 = ima.logical_dot_sample(im_test, im_FISH1_template, 100)
im_FISH2_template = ima.logical_ellipse(im_test, 125, 125, 4, 4)
im_FISH2_template = ima.logical_ellipse(im_FISH2_template, 90, 100, 4, 5)
im_FISH2_template = ima.logical_ellipse(im_FISH2_template, 130, 150, 5, 5)
im_FISH2 = ima.logical_dot_sample(im_test, im_FISH2_template, 100)

print(np.sum(im_FISH_template))
print(np.sum(im_FISH1_template))
print(np.sum(im_FISH2_template))

_, r, g, _ = mat.auto_correlation(im_FISH, im_seg, 100)
_, r1, g1, _ = mat.auto_correlation(im_FISH1, im_seg, 100)
_, r2, g2, _ = mat.auto_correlation(im_FISH2, im_seg, 100)
print(g)
print(g1)
print(g2)

plt.subplots(figsize=(6, 4))
plt.plot(r[:60], g[:60], color='#FFD700', label='1 cluster')
plt.plot(r1[:60], g1[:60], alpha=0.4, color='#40E0D0', label='2 clusters')
plt.plot(r2[:60], g2[:60], alpha=0.4, color='#DA70D6', label='3 clusters')
#plt.ylim([-0.5, 5])
plt.xlabel('r')
plt.ylabel('g')
plt.legend(loc=2, bbox_to_anchor=(0.02, 0.99))
plt.savefig('%s/auto_correlation_test.pdf' % data_path)
plt.close()"""

# step 12: effect of size of cell
"""im_test = np.zeros((250, 250))
im_seg = ima.logical_ellipse(im_test, 125, 125, 50, 50)
im_FISH_template = ima.logical_ellipse(im_test, 125, 125, 8, 8)
im_FISH = ima.logical_dot_sample(im_test, im_FISH_template, 100)
im_seg1 = ima.logical_ellipse(im_test, 125, 125, 60, 60)
im_seg2 = ima.logical_ellipse(im_test, 125, 125, 70, 70)

_, r, g, _ = mat.auto_correlation(im_FISH, im_seg, 100)
_, r1, g1, _ = mat.auto_correlation(im_FISH, im_seg1, 100)
_, r2, g2, _ = mat.auto_correlation(im_FISH, im_seg2, 100)
print(g)
print(g1)
print(g2)

plt.subplots(figsize=(6, 4))
plt.plot(r[:60], g[:60], color='#FFD700', label='R50')
plt.plot(r1[:60], g1[:60], alpha=0.4, color='#40E0D0', label='R60')
plt.plot(r2[:60], g2[:60], alpha=0.4, color='#DA70D6', label='R70')
#plt.ylim([-0.5, 5])
plt.xlabel('r')
plt.ylabel('g')
plt.legend(loc=2, bbox_to_anchor=(0.02, 0.99))
plt.savefig('%s/auto_correlation_test.pdf' % data_path)
plt.close()"""

# test13 different ecDNA copy number
im_test = np.zeros((250, 250))
im_seg = ima.logical_ellipse(im_test, 125, 125, 50, 50)
im_FISH = ima.logical_ellipse(im_test, 125, 125, 5, 5)
im_FISH1 = ima.logical_ellipse(im_FISH, 90, 100, 5, 5)
im_FISH2 = ima.logical_ellipse(im_FISH1, 130, 150, 5, 5)

_, r, g, _ = mat.auto_correlation(im_FISH, im_seg, 100)
_, r1, g1, _ = mat.auto_correlation(im_FISH1, im_seg, 100)
_, r2, g2, _ = mat.auto_correlation(im_FISH2, im_seg, 100)
print(g)
print(g1)
print(g2)

plt.subplots(figsize=(6, 4))
plt.plot(r[:60], g[:60], color='#FFD700', label='1 cluster')
plt.plot(r1[:60], 2*np.array(g1[:60]), alpha=0.4, color='#40E0D0', label='2 clusters')
plt.plot(r2[:60], 3*np.array(g2[:60]), alpha=0.4, color='#DA70D6', label='3 clusters')
#plt.ylim([-0.5, 5])
plt.xlabel('r')
plt.ylabel('g')
plt.legend(loc=2, bbox_to_anchor=(0.02, 0.99))
plt.savefig('%s/auto_correlation_test.pdf' % data_path)
plt.close()

viewer = napari.Viewer()
viewer.add_image(im_seg, colormap='red')
napari.run()

viewer = napari.Viewer()
viewer.add_image(im_FISH, colormap='green', contrast_limits=[0, 1])
napari.run()

viewer = napari.Viewer()
viewer.add_image(im_seg, colormap='red', blending='additive')
viewer.add_image(im_FISH, colormap='green', blending='additive', contrast_limits=[0, 1])
napari.run()

viewer = napari.Viewer()
viewer.add_image(im_seg1, colormap='red', contrast_limits=[0, 1])
napari.run()

viewer = napari.Viewer()
viewer.add_image(im_seg1, colormap='red', blending='additive')
viewer.add_image(im_FISH, colormap='green', blending='additive', contrast_limits=[0, 1])
napari.run()

viewer = napari.Viewer()
viewer.add_image(im_seg2, colormap='red', contrast_limits=[0, 1])
napari.run()

viewer = napari.Viewer()
viewer.add_image(im_seg2, colormap='red', blending='additive')
viewer.add_image(im_FISH, colormap='green', blending='additive', contrast_limits=[0, 1])
napari.run()