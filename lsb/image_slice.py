# cut image in boxes only where there is data
from __future__ import print_function
from __future__ import division
import os
import pdb
from tempfile import mkstemp

from astropy import units
from astropy.coordinates import SkyCoord
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from copy import deepcopy

from astropy.wcs import WCS

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


def make_stamp(img, wei, x0, x1, y0, y1, wcs=None, filter_images=None, color_image=None, color_image_wcs=None):
    # if ((img[x0:x1, y0:y1] > 0).sum() / len(img[x0:x1, y0:y1])) > 0.99 and (wei[x0:x1, y0:y1] > 0).sum() / len(
    #         wei[x0:x1, y0:y1]) > 0.99:

    # Check limits
    if x0 >= x1 or y0 >= y1:
        raise ValueError("must: x0 < x1 and y0 < y1")

    fd, img_file = mkstemp(dir=os.path.expanduser("~/data/tmp_stamps/"))
    os.close(fd)

    # Temporary file names
    os.unlink(img_file)
    wei_file = img_file + ".wei.fits"
    img_file = img_file + ".fits"
    seg_file = wei_file.replace("wei", "seg")

    img_stamp = img[x0:x1, y0:y1]
    wei_stamp = wei[x0:x1, y0:y1]

    # Test if sum(wei < thesh) is less than 90%
    mask = wei_stamp <= 1e-4
    invalid = mask.sum() / mask.size
    if invalid > 0.90 or len(mask) < 1:
        print("Skipping %.3f..." % invalid)
        return False

    if wcs is not None:
        wcs = deepcopy(wcs)
        crval = wcs.all_pix2world(y0, x0, 0)
        wcs.wcs.crpix = [0, 0]
        wcs.wcs.crval = crval
        hdr = wcs.to_header()
    else:
        hdr = None

    # img_stamp, wei_stamp = remove_borders(img_stamp, wei_stamp)

    # try:
    fits.writeto(img_file, img_stamp, header=hdr)
    fits.writeto(wei_file, wei_stamp, header=hdr)
    # except IOError:
    #     pass

    # Run SExtractor

    if not os.path.exists(seg_file):
        kwargs = {'code': 'SExtractor',
                  'config_file': os.path.expanduser('~/workspace/lsb/lsb_project/lsb/config.sex'),
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

    stamp_coord = wcs.all_pix2world([[(y1 - y0) / 2 + 0.5, (x1 - x0) / 2 + 0.5]], 0)  # Central coordinate
    if wcs is not None and filter_images is not None:
        print("coord: ", SkyCoord(stamp_coord[0][0] * units.deg, stamp_coord[0][1] * units.deg).to_string('hmsdms'))
        filter_files = dict()
        filter_files_weight = dict()
        for filter_name in filter_images:
            exts = fits.open(filter_images[filter_name])
            exts_wei = fits.open(
                filter_images[filter_name].replace('.fits.fz', '.weight.fits.fz'))  # FIXME: too hardocoded
            for i_ext, e in enumerate(exts[1:]):
                w = WCS(e.header)
                # print([stamp_coord[0][0]], [stamp_coord[0][1]], 0)
                pix = w.all_world2pix([stamp_coord[0][0]], [stamp_coord[0][1]], 0)
                pix = pix[1][0], pix[0][0]
                if 0 < pix[0] < e.shape[0] and 0 < pix[1] < e.shape[1]:
                    dx, dy = (x1 - x0) / 2, (y1 - y0) / 2
                    xx0 = int(pix[0] - dx) if pix[0] - dx > 0 else 0
                    yy0 = int(pix[1] - dy) if pix[1] - dy > 0 else 0
                    xx1 = int(pix[0] + dx) if pix[0] + dx < img.shape[0] else img.shape[0]
                    yy1 = int(pix[1] + dy) if pix[1] + dy < img.shape[1] else img.shape[1]
                    out_fname = img_file.replace('.fits', '%s.fits' % filter_name)

                    # slice
                    print(e.header['EXTNAME'], xx0, xx1, yy0, yy1)
                    filter_stamp = e.data[xx0:xx1, yy0:yy1]
                    filter_stamp_weight = exts_wei[i_ext + 1].data[xx0:xx1, yy0:yy1]

                    # Do the wcs
                    w = deepcopy(w)
                    crval = w.all_pix2world(yy0, xx0, 0)
                    w.wcs.crpix = [0, 0]
                    w.wcs.crval = crval
                    hdr = w.to_header()

                    # write image
                    fits.writeto(out_fname, filter_stamp, header=hdr)
                    filter_files.update({filter_name: out_fname})
                    # write weight
                    out_fname = out_fname.replace(".fits", "wei.fits")
                    fits.writeto(out_fname, filter_stamp_weight, header=hdr)
                    filter_files_weight.update({filter_name: out_fname})

                    # break
            exts.close()
            exts_wei.close()

        ret = img_file, wei_file, seg_file, filter_files, filter_files_weight
    else:

        ret = img_file, wei_file, seg_file

    if color_image is not None and color_image_wcs is not None:
        stamp_coord = wcs.all_pix2world([[0, 0]], 0)  # Stamp coordinate
        pix = color_image_wcs.all_world2pix([stamp_coord[0][0]], [stamp_coord[0][1]], 0)
        print(stamp_coord, pix)
        xx0, yy0 = pix[1][0], pix[0][0]
        xx0 = xx0 if xx0 > 0 else 0
        yy0 = yy0 if yy0 > 0 else 0
        xx1 = xx0 + img_stamp.shape[0] if xx0 + img_stamp.shape[0] < color_image.shape[0] else color_image.shape[0]
        yy1 = yy0 + img_stamp.shape[1] if yy0 + img_stamp.shape[1] < color_image.shape[1] else color_image.shape[1]
        # print(xx0,xx1,yy0,yy1)
        # pdb.set_trace()
        color_stamp = color_image[int(xx0):int(xx1), int(yy0):int(yy1)]

        ret = list(ret)
        ret += [color_stamp]

    return ret


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
