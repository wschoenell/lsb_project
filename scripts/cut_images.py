import matplotlib
from astropy.wcs import WCS

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import os

import sys
from astropy.io import fits
import numpy as np
from matplotlib.colors import LogNorm
from scipy import ndimage

try:
    from lsb.image_slice import make_stamp, remove_borders
except ImportError:
    sys.path.append('/Users/william/PyCharmProjects/lsb_project/')
    from lsb.image_slice import make_stamp, remove_borders

work_dir = os.path.expanduser("~/data/ngc3115/work/")
bins_x = 50
bins_y = 50

for dirpath, dirnames, filenames in os.walk(work_dir):
    # if "2015A-0616" not in dirpath:
    for filename in filenames:
        if filename.startswith('coadd') and filename.endswith("fz") and not 'weight' in filename:

            print("Slicing file %s" % filename)

            image_file = os.path.join(dirpath, filename)
            weight_file = image_file.replace("fits.fz", "weight.fits.fz")

            img = fits.getdata(image_file)
            wcs = WCS(fits.getheader(image_file, ext=1))
            wei = fits.getdata(weight_file)

            # img = remove_borders(img)

            x0 = 0
            y0 = 0
            i_bin = 0
            bins = []

            for yy in np.linspace(0, img.shape[1], bins_y, dtype=int)[1:-1]:
                for xx in np.linspace(0, img.shape[0], bins_x, dtype=int)[1:-1]:
                    bins.append([x0, yy, y0, yy])
                    print(x0, xx, y0, yy)
                    img_file, wei_file, seg_file = make_stamp(img, wei, x0, xx, y0, yy)

                    img_stamp = fits.getdata(img_file)
                    wei_stamp = fits.getdata(wei_file)

                    mask = wei_stamp <= 1e-4

                    invalid = mask.sum() / mask.size
                    if invalid > 0.90 or len(mask) < 1:
                        print("Skipping %.3f..." % invalid)
                    else:

                        seg_stamp = fits.getdata(seg_file)

                        # Write a fits file with the stamp
                        imageHeader = fits.Header()
                        hdu = fits.CompImageHDU(img_stamp, imageHeader)
                        hdu.writeto("%s/stamps/stamp_%04i_img.fits.fz" % (work_dir, i_bin), overwrite=True)

                        title = filename.replace('.fits.fz', '').replace('coadd_', '%i  ' % i_bin)
                        # Figure 1: Normal image w/o scaling
                        fsize = img_stamp.shape[1] / 100., img_stamp.shape[0] / 100.
                        fig = plt.figure(8, figsize=fsize, dpi=100)
                        plt.clf()
                        ax = fig.add_axes([0, 0, 1, 1])
                        ax.axis('off')
                        plt.imshow(np.ma.masked_array(np.log10(img_stamp[::-1]), mask=mask[::-1]), cmap=plt.cm.gray_r)
                        # plt.title(title)
                        plt.savefig("%s/stamps/stamp_%04i_1.png" % (work_dir, i_bin))
                        # with open("%s/stamps/stamp_%04i_1.png" % (work_dir, i_bin), 'w') as outfile:
                        #     fig.canvas.print_png(outfile)

                        # Figure 2: Image with agressive scaling
                        mask += seg_stamp > 0

                        fig = plt.figure(2, figsize=fsize, dpi=100)
                        plt.clf()
                        ax = fig.add_axes([0, 0, 1, 1])
                        ax.axis('off')
                        aux = np.ma.masked_array(img_stamp, mask=mask)
                        if aux.min() < 0:
                            aux -= aux.min()
                        median = np.ma.median(aux)
                        stdev = np.ma.std(aux)
                        vmin = 0 if median < stdev else median - stdev
                        plt.imshow(img_stamp[::-1], cmap=plt.cm.gray_r, vmin=vmin, vmax=(median + stdev),
                                   norm=LogNorm())  #
                        plt.savefig("%s/stamps/stamp_%04i_2.png" % (work_dir, i_bin))

                        # Figure 3: Image with a gaussian filter applied on it
                        fig = plt.figure(3, figsize=fsize, dpi=100)
                        plt.clf()
                        ax = fig.add_axes([0, 0, 1, 1])
                        ax.axis('off')
                        img_stamp[mask] = np.random.normal(loc=median, scale=stdev)

                        img_g = ndimage.gaussian_filter(img_stamp, sigma=(3, 3), order=0)
                        median = np.median(img_g)
                        stdev = np.std(img_g)
                        vmin = 0 if median < stdev else median - stdev

                        plt.imshow(np.ma.masked_array(img_g, mask=mask)[::-1], cmap=plt.cm.hot_r, vmin=vmin,
                                   vmax=(median + stdev))
                        plt.savefig("%s/stamps/stamp_%04i_3.png" % (work_dir, i_bin))

                        i_bin += 1

                        # raw_input('next...')

                    x0 = xx
                y0 = yy

            # plt.clf()
            # print 'n_bins', n_bins
            # img[wei == 0] = 0
            # plt.imshow(np.log10(img), cmap=plt.cm.hot, origin='lower')
            # plt.hlines(np.linspace(0, img.shape[0], bins_x, dtype=int)[1:], *plt.xlim())
            # plt.vlines(np.linspace(0, img.shape[1], bins_y, dtype=int)[1:], *plt.ylim())
            # plt.show()
