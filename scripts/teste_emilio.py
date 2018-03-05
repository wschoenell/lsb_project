import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from matplotlib.colors import LogNorm

hdr = fits.getheader("/Users/william/data/ngc3115/work/coadd_10:07:38.23_-7:07:07.8.fits.fz", ext=1)
img = fits.getdata("/Users/william/data/ngc3115/work/coadd_10:07:38.23_-7:07:07.8.fits.fz")
wcs = WCS(hdr)

datapoints = np.loadtxt("/Users/william/Downloads/rej_cat3115.dat")
pixels = wcs.all_world2pix(datapoints[:, 0], datapoints[:, 1], 0)

x0, x1, y0, y1 = 4784, 7710, 21500, 24481
img = img[x0:x1, y0:y1]


if img.min() < 0:
    img -= img.min()

median = np.median(img[1400:1800,400:800])
std = np.std(img[1400:1800,400:800])
# plt.figure(2)
# plt.clf()
# plt.imshow(np.log10(img[1400:1800,400:800]))

plt.clf()
plt.imshow(img.T, norm=LogNorm(), cmap=plt.cm.gray_r, vmin=median-std, vmax=median+std)
plt.scatter(pixels[1] - x0, pixels[0] - y0, c=datapoints[:, 2], alpha=0.5)
plt.colorbar()
