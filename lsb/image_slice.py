# cut image in boxes only where there is data
from __future__ import print_function
from __future__ import division
import os
from tempfile import mkstemp

from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt

plt.ion()

import astromatic_wrapper as aw


def remove_borders(img, wei=None, box_size=3, tresh=0.1):
    """
    Remove image borders by a running blank detector on the four corners

    box_size: Box size to calculate percentage of 0's
    tresh: Stop when, i.e., 10% is non-zero
    """

    x0, y0 = 0, 0
    x1, y1 = img.shape
    # x axis: from left to right
    for i in np.arange(0, img.shape[0], box_size):
        if (img[i:i + box_size] > 0).sum() / (img.shape[1] * box_size) > tresh:
            x0 = i
            break
    # x axis: right to left
    for i in np.arange(0, img.shape[0], box_size)[::-1]:
        if (img[i - box_size:i] > 0).sum() / (img.shape[1] * box_size) > tresh:
            x1 = i
            break
    # y axis: bottom to up
    for i in np.arange(0, img.shape[1], box_size):
        if (img[x0:x1, i:i + box_size] > 0).sum() / (img.shape[0] * box_size) > tresh:
            y0 = i
            break
    # y axis: up to bottom
    for i in np.arange(0, img.shape[1], box_size)[::-1]:
        if (img[x0:x1, i - box_size:i] > 0).sum() / (img.shape[0] * box_size) > tresh:
            y1 = i
            break
    print("Removing image borders: (%i,%i:%i,%i)" % (x0, x1, y0, y1))

    if wei is None:
        return img[x0:x1, y0:y1]
    else:
        return img[x0:x1, y0: y1], wei[x0:x1, y0:y1]


def make_stamp(img, wei, x0, x1, y0, y1):
    # if ((img[x0:x1, y0:y1] > 0).sum() / len(img[x0:x1, y0:y1])) > 0.99 and (wei[x0:x1, y0:y1] > 0).sum() / len(
    #         wei[x0:x1, y0:y1]) > 0.99:
    _, img_file = mkstemp(dir=os.path.expanduser("~/data/tmp_stamps/"))

    # Temporary file names
    wei_file = img_file + ".wei.fits"
    img_file = img_file + ".fits"
    seg_file = wei_file.replace("wei", "seg")

    img_stamp = img[x0:x1, y0:y1]
    wei_stamp = wei[x0:x1, y0:y1]

    img_stamp, wei_stamp = remove_borders(img_stamp, wei_stamp)

    try:
        fits.writeto(img_file, img_stamp)
        fits.writeto(wei_file, wei_stamp)
    except IOError:
        pass

    # Run SExtractor

    if not os.path.exists(seg_file):
        kwargs = {'code': 'SExtractor',
                  'config_file': os.path.expanduser('~/PyCharmProjects/lsb_project/lsb/config.sex'),
                  'cmd': '/usr/local/bin/sex',
                  'temp_path': '.',
                  'config': {'WEIGHT_TYPE': 'MAP_WEIGHT',
                             'WEIGHT_IMAGE': wei_file,
                             "FILTER_NAME": "/usr/local/share/sextractor/default.conv",
                             "CHECKIMAGE_TYPE": "SEGMENTATION",
                             "CHECKIMAGE_NAME": seg_file
                             },
                  "params": ["NUMBER", "ALPHA_J2000", "DELTA_J2000", "XMIN_IMAGE", "YMIN_IMAGE", "FLUX_AUTO",
                             "FLUXERR_AUTO",
                             "ELLIPTICITY", "FLAGS"]}
        sex = aw.api.Astromatic(**kwargs)
        sex.run(img_file)

    return img_file, wei_file, seg_file


if __name__ == '__main__':

    # image_file = os.path.expanduser("~/data/lsb_work/coadd_12:42:13.18_-34:25:39.3.fits.fz")
    image_file = os.path.expanduser("/Users/william/data/lsb/coadd_10:02:50.03_-8:19:04.0.fits.fz")
    weight_file = image_file.replace("fits.fz", "weight.fits.fz")

    bins_x = 50
    bins_y = 50

    img = fits.getdata(image_file)
    wei = fits.getdata(weight_file)

    x0 = 0
    y0 = 0
    n_bins = 0
    bins = []

    for yy in np.linspace(0, img.shape[1], bins_y, dtype=int)[1:-1]:
        for xx in np.linspace(0, img.shape[0], bins_x, dtype=int)[1:-1]:
            bins.append([x0, yy, y0, yy])
            print(x0, xx, y0, yy)
            img_file, wei_file, seg_file = make_stamp(img, wei, x0, xx, y0, yy)

            img_stamp = fits.getdata(img_file)
            # seg_stamp = fits.getdata(seg_file)
            wei_stamp = fits.getdata(wei_file)

            mask = wei_stamp <= 1e-4

            invalid = mask.sum() / mask.size
            if invalid > 0.90:
                print("Skipping %.3f..." % invalid)
            else:

                plt.figure(1)
                plt.clf()
                plt.imshow(np.ma.masked_array(np.log10(img_stamp), mask=mask),
                           cmap=plt.cm.gray_r)  # , vmin=np.log10(0.1*sigma), vmax=np.log10(3*sigma))
                # plt.imshow(np.log10(img_stamp), cmap=plt.cm.gray)
                plt.draw()
                plt.show()
                raw_input('next...')

            x0 = xx
        y0 = yy

    # plt.clf()
    # print 'n_bins', n_bins
    # img[wei == 0] = 0
    # plt.imshow(np.log10(img), cmap=plt.cm.hot, origin='lower')
    # plt.hlines(np.linspace(0, img.shape[0], bins_x, dtype=int)[1:], *plt.xlim())
    # plt.vlines(np.linspace(0, img.shape[1], bins_y, dtype=int)[1:], *plt.ylim())
    # plt.show()
