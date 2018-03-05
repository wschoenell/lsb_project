from __future__ import print_function
import os

import sys
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy import units as u
from astropy.io.ascii import Csv
import numpy as np
from astropy.wcs import WCS
from scipy import ndimage

from matplotlib.colors import LogNorm

try:
    from lsb.image_slice import make_stamp
except ImportError:
    sys.path.append('/Users/william/PyCharmProjects/lsb_project/')
    from lsb.image_slice import make_stamp

files_dir = os.path.expanduser("~/data/lsb_work/")

import matplotlib.pyplot as plt

plt.ion()

pointings = []

for dirpath, dirnames, filenames in os.walk(files_dir):
    for filename in filenames:
        if not 'weight' in filename and filename.endswith("fits.fz"):
            f = os.path.join(dirpath, filename)
            hdr = fits.getheader(f, ext=1)
            pointings.append([os.path.join(dirpath, filename), WCS(hdr), hdr])
            # print " ".join(pointings[-1])

# img_coords = SkyCoord(["%s %s" % (p[1], p[2]) for p in pointings], unit=(u.deg, u.deg))
img_coords = SkyCoord([SkyCoord.from_pixel(p[2]['NAXIS1'] / 2, p[2]['NAXIS2'] / 2, p[1]) for p in pointings])
c = Csv()
# catalog_lsb = c.read(os.path.expanduser("~/workspace/lsb_project/docs/Muller_candidates.csv"))
catalog_lsb = c.read(os.path.expanduser("~/PyCharmProjects/lsb_project/docs/Muller_candidates.csv"))
lsb_coords = SkyCoord(["%s %s" % (p['ra'], p['dec']) for p in catalog_lsb], unit=(u.hour, u.deg))

idxsearcharound, idxself, sep2d, dist3d = lsb_coords.search_around_sky(img_coords, seplimit=1.0 * u.deg)

# grab the nearest neighbor
mask = np.zeros_like(idxself, dtype=bool)
mask[np.ravel(
    [np.argwhere(idxself == x)[np.argmin(sep2d[np.argwhere(idxself == x)])] for x in np.unique(idxself)])] = True
# idxsearcharound, idxself, sep2d, dist3d = idxsearcharound[mask], idxself[mask], sep2d[mask], dist3d[mask]

print("Found %i matches" % len(idxself))

for i_img in np.unique(idxsearcharound):
    for i_lsb in idxself[idxsearcharound == i_img]:
        xycoord = lsb_coords[i_lsb].to_pixel(pointings[i_img][1])
        xycoord = [int(i) for i in xycoord]
        # print(xycoord)
        # print(catalog_lsb[i_lsb], pointings[i_img])
        img = fits.getdata(pointings[i_img][0])
        wei = fits.getdata(pointings[i_img][0].replace('.fits.fz', '.weight.fits.fz'))

        img_file, wei_file, seg_file = make_stamp(img, wei, xycoord[1] - 150, xycoord[1] + 150, xycoord[0] - 150,
                                                  xycoord[0] + 150)

        img_stamp = fits.getdata(img_file)
        wei_stamp = fits.getdata(wei_file)
        seg_stamp = fits.getdata(seg_file)

        os.rename(img_file, "%s/%s.img.fits" % (os.path.dirname(img_file), catalog_lsb[i_lsb]["Name"]))
        os.rename(wei_file, "%s/%s.wei.fits" % (os.path.dirname(wei_file), catalog_lsb[i_lsb]["Name"]))
        os.rename(seg_file, "%s/%s.seg.fits" % (os.path.dirname(seg_file), catalog_lsb[i_lsb]["Name"]))

        mask = np.bitwise_or(seg_stamp > 0, wei_stamp < 1e-2)

        if img_stamp.min() < 0:
            img_stamp -= img_stamp.min()
        median = np.median(img_stamp[~mask])
        stdev = np.std(img_stamp[~mask])

        # img_stamp = np.ma.masked_array(img_stamp, mask=mask)
        # img_stamp[mask] = 0
        # wei_stamp[mask] = 0
        # wei_stamp = np.ma.masked_array(wei_stamp, mask=mask)

        fig = plt.figure(1)
        plt.clf()
        # plt.title(catalog_lsb[i_lsb]["Name"])

        ax = fig.add_axes([0.1, 0.5, 0.4, 0.4])
        ax.axes.get_xaxis().set_visible(False)
        ax.axes.get_yaxis().set_visible(False)

        vmin = 0 if median < stdev else median - stdev
        ax.imshow(img_stamp, cmap=plt.cm.gray_r, norm=LogNorm(), vmin=vmin, vmax=(median + stdev))

        ax.set_title(catalog_lsb[i_lsb]["Name"])

        img_stamp[mask] = np.random.normal(loc=median, scale=stdev, size=len(img_stamp[mask]))
        if img_stamp.min() < 0:
            img_stamp -= img_stamp.min()

        ax = fig.add_axes([0.5, 0.5, 0.4, 0.4], sharey=fig.axes[-1])
        ax.axes.get_xaxis().set_visible(False)
        ax.axes.get_yaxis().set_visible(False)
        img_g = ndimage.gaussian_filter(img_stamp, sigma=(2, 2), order=0)
        ax.imshow(np.ma.masked_array(img_g, mask=mask), cmap=plt.cm.hot_r, vmin=vmin, vmax=(median + stdev))

        ax = fig.add_axes([0.1, 0.1, 0.4, 0.4])

        ax.axes.get_xaxis().set_visible(False)
        ax.axes.get_yaxis().set_visible(False)
        ax.imshow(img_stamp, cmap=plt.cm.gray_r, norm=LogNorm(), vmin=vmin, vmax=(median + stdev))

        ax = fig.add_axes([0.5, 0.1, 0.4, 0.4], sharey=fig.axes[-1])
        ax.axes.get_xaxis().set_visible(False)
        ax.axes.get_yaxis().set_visible(False)

        ax.imshow(np.ma.masked_array(img_stamp, mask=mask), cmap=plt.cm.gray_r, norm=LogNorm(), vmin=vmin,
                  vmax=(median + stdev))

        try:
            plt.savefig("%s.png" % catalog_lsb[i_lsb]["Name"])
        except ValueError:
            print("Skipped %s" % catalog_lsb[i_lsb]["Name"])

        # plt.figure(5)
        # mask += img_stamp > (1.1 * sigma)
        # img_stamp[mask] = np.random.normal(loc=sigma, scale=rms, size=len(img_stamp[mask]))
        # img_g = ndimage.gaussian_filter(img_stamp, sigma=(4, 4), order=0)
        # plt.imshow(img_g, cmap=plt.cm.hot_r)  # , vmin=np.log10(sigma-2*rms), vmax=np.log10(sigma+2*rms))
        # plt.imshow(np.ma.masked_array(img_stamp, mask=mask), cmap=plt.cm.hot_r)  # , vmin=np.log10(sigma - 20 * rms),

        # plt.figure(2)
        # plt.clf()
        # plt.imshow(np.log10(img_stamp), cmap=plt.cm.hot_r, vmin=sigma*1, vmax=sigma*.1)
        # plt.title(catalog_lsb[i_lsb]["Name"])
        # plt.show()
        # plt.savefig("%s.png" % catalog_lsb[i_lsb]["Name"])

        # raw_input("next...")
