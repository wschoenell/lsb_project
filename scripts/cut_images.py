from __future__ import print_function
import os

import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.image as mpimg

from plot_lib import make_plots

import pickle

try:
    from lsb.image_slice import make_stamp, remove_borders
    from lsb.util import label_pc
except ImportError:
    sys.path.append('/Users/wschoenell/home_william/workspace/lsb/lsb_project/')
    from lsb.image_slice import make_stamp, remove_borders


from astropy.wcs import WCS

from astropy.io import fits
import numpy as np


work_dir = os.path.expanduser("~/data/ngc3115/work/")


stamp_size = [800, 800]  # 10 * 2 * 500pc x 500pc radius ~ 800 px
stamp_margin = [40, 40]  # Overlap between stamps


# pdb.set_trace()

i_bin = 0
bins = dict()
for dirpath, dirnames, filenames in os.walk(work_dir):
    # if "2015A-0616" not in dirpath:
    for filename in filenames:
        if filename.startswith('coadd') and filename.endswith(
                "fz") and not 'weight' in filename and not 'g' in filename \
                and not 'r' in filename:

            print("Slicing file %s" % filename)

            image_file = os.path.join(dirpath, filename)
            weight_file = image_file.replace("fits.fz", "weight.fits.fz")
            filter_images = {f: image_file.replace('fits.fz', '%s.fits.fz' % f) for f in ['g', 'r']}
            color_image = image_file.replace('fits.fz', 'png')
            color_image_wcs = image_file.replace('fits.fz', 'pngwcs.fits')

            img = fits.getdata(image_file)
            wcs = WCS(fits.getheader(image_file, ext=1))
            wei = fits.getdata(weight_file)

            color = np.flip(mpimg.imread(color_image), 0)  # PNG should be flipped along axis=0 to match the WCS
            color_wcs = WCS(fits.getheader(color_image_wcs))

            y0 = 0
            for yy in np.arange(0, img.shape[1], stamp_size[1], dtype=int)[1:-1]:
                # # For testing:
                # if i_bin > 100:
                #     break
                x0 = 0
                for xx in np.arange(0, img.shape[0], stamp_size[0], dtype=int)[1:-1]:
                    file_prefix = "%s/stamps/%04i/stamp_%04i" % (work_dir, i_bin, i_bin)
                    bins[i_bin] = [filename, [x0, yy, y0, yy], None, False,
                                   file_prefix]  # filename, (x0,x1,y0,y1), wcs, exists?, file_prefix

                    aux_ox0 = 0 if x0 < stamp_margin[0] else x0 - stamp_margin[0]
                    aux_ox1 = img.shape[0] if xx + stamp_margin[0] > img.shape[0] else xx + stamp_margin[0]

                    aux_oy0 = 0 if y0 < stamp_margin[1] else y0 - stamp_margin[1]
                    aux_oy1 = img.shape[1] if yy + stamp_margin[1] > img.shape[1] else yy + stamp_margin[1]

                    stamp = make_stamp(img, wei, aux_ox0, aux_ox1, aux_oy0, aux_oy1, wcs=wcs,
                                       filter_images=filter_images, color_image=color, color_image_wcs=color_wcs)
                    if not stamp:
                        print("Skipped...")
                    else:

                        img_file, wei_file, seg_file, filter_files, filter_files_weight, color_stamp = stamp

                        img_stamp = fits.getdata(img_file)
                        wei_stamp = fits.getdata(wei_file)

                        # Create stamp dir if not exists
                        if not os.path.exists("%s/stamps/%04i/" % (work_dir, i_bin)):
                            os.mkdir("%s/stamps/%04i/" % (work_dir, i_bin))

                        # Write a fits file with the stamp
                        imageHeader = fits.getheader(img_file)
                        hdu = fits.CompImageHDU(img_stamp, imageHeader)
                        hdu.writeto("%s_img.fits.fz" % file_prefix, overwrite=True)

                        # Store stamp WCS on output file list
                        image_wcs = WCS(imageHeader)
                        bins[i_bin][2] = image_wcs
                        # Set exists to True.
                        bins[i_bin][3] = True

                        # Write a fits file with the single-filter images
                        for fname in filter_files:
                            imageHeader = fits.getheader(filter_files[fname])
                            imageData = fits.getdata(filter_files[fname])
                            hdu = fits.CompImageHDU(imageData, imageHeader)
                            hdu.writeto("%s_%s_img.fits.fz" % (file_prefix, fname),
                                        overwrite=True)

                        for fname in [filename]:
                            # title = fname.replace('.fits.fz', '').replace('coadd_', '%04i  ' % i_bin).strip()
                            make_plots(img_file, wei_file=wei_file, seg_file=seg_file, wcs=image_wcs,
                                       color_stamp=color_stamp, file_prefix=file_prefix, stamp_margin=stamp_margin)

                        # For the individual filters
                        # for filter_name in filter_files:
                        #     print(filter_files[filter_name])
                        #     make_plots(filter_files[filter_name], wei_file=filter_files_weight[filter_name],
                        #                suffix="_%s" % filter_name)

                        # raw_input('next...')
                        # Delete temporary files
                        for fname in [img_file, wei_file, seg_file] + list(filter_files.values()) + list(
                                filter_files_weight.values()):
                            os.unlink(fname)

                    i_bin += 1

                    x0 = xx
                y0 = yy

            # plt.clf()
            # print 'n_bins', n_bins
            # img[wei == 0] = 0
            # plt.imshow(np.log10(img), cmap=plt.cm.hot, origin='lower')
            # plt.hlines(np.linspace(0, img.shape[0], bins_x, dtype=int)[1:], *plt.xlim())
            # plt.vlines(np.linspace(0, img.shape[1], bins_y, dtype=int)[1:], *plt.ylim())
            # plt.show()

pickle.dump(bins, open("%s/stamps/stamp_list.pkl" % work_dir, "wb"))
