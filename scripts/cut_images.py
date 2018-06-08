from __future__ import print_function
import os
import pdb

import sys
import time

import matplotlib
import matplotlib.image as mpimg
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.wcs.utils import skycoord_to_pixel
from astroquery.simbad import Simbad

Simbad = Simbad()
Simbad.add_votable_fields('otype')

matplotlib.use("Agg")
import pickle

try:
    from lsb.image_slice import make_stamp, remove_borders
    from lsb.util import label_pc
except ImportError:
    sys.path.append('/Users/william/workspace/lsb/lsb_project/')
    from lsb.image_slice import make_stamp, remove_borders

from lsb.util import label_pc
from lsb.patches import circles

from astropy.wcs import WCS
from matplotlib.patches import Circle

import matplotlib.pyplot as plt
from astropy.io import fits
import numpy as np
from matplotlib.colors import LogNorm
from scipy import ndimage

import astromatic_wrapper as aw

work_dir = os.path.expanduser("~/data/ngc3115/work/")
z = 0.002212  # NGC3115 redshift
plate_scale = 0.263  # DECam

typical_radius = {i: label_pc(z, i).value / plate_scale for i in [70, 500, 1000]}  # in pixels

stamp_size = [800, 800]  # 10 * 2 * 500pc x 500pc radius ~ 800 px
stamp_margin = [40, 40]  # Overlap between stamps


# pdb.set_trace()

def make_plots(img_file, wei_file=None, seg_file=None, suffix='', wcs=None, color_stamp=None, sextractor_plot=False,
               plot_catalog_color=True):
    img_stamp, header_stamp = fits.getdata(img_file, header=True)

    # Load weight file
    if wei_file is not None:
        wei_stamp = fits.getdata(wei_file)
        mask = wei_stamp <= 1e-4
    else:
        mask = np.zeros_like(img_stamp, dtype=bool)

    # Figure 1: Normal image w/o scaling
    # FIXME: fsize is not working anymore!!!
    fsize = img_stamp.shape[1] / 100., img_stamp.shape[0] / 100.
    fig = plt.figure(1, figsize=fsize, dpi=100)
    plt.clf()
    ax = fig.add_axes([0, 0, 1, 1])
    ax.axis('off')
    aux = np.ma.masked_array(np.log10(img_stamp[::-1]), mask=mask[::-1]).compressed()
    # vmin, vmax = np.percentile(aux, [15, 95])
    plt.imshow(np.ma.masked_array(np.log10(img_stamp[::-1]), mask=mask[::-1]),
               cmap=plt.cm.gray_r)  # , vmin=vmin, vmax=vmax)

    aux_r = 2 * min(typical_radius.values()) + 10
    r_max = max(typical_radius.values())

    for radius in np.sort(typical_radius.keys()):
        # print(aux_r, r_max + 10, typical_radius[radius])
        circ = Circle((aux_r + typical_radius[radius], r_max + 10), typical_radius[radius], ec="black", alpha=0.8,
                      fill=False)

        ax = plt.gca()
        ax.add_patch(circ)
        ax.text(aux_r + typical_radius[radius], r_max + typical_radius[radius] + 20, "%i" % radius, color='red',
                ha="center", va="center")

        aux_r += 2 * typical_radius[radius] + 10

    ### ADD KNOWN OBJECTS FROM Simbad ###
    catalog = None
    if wcs:
        coordinate = wcs.all_pix2world([[0, 0]], 0)
        while True:
            try:
                catalog = Simbad.query_region(SkyCoord(coordinate, unit='deg'),
                                              10 * u.arcmin)  # FIXME: radius is set to a constant
                break
            except:
                time.sleep(2)
                print("Simbad error... Retrying...")
        print("Catalog length" + str(len(catalog)))
        print(img_stamp.shape)
        for obj in catalog:
            # pdb.set_trace()
            x, y = skycoord_to_pixel(SkyCoord([[obj["RA"], obj["DEC"]]], unit=('hour', 'deg')), wcs)
            x, y = float(x), float(y)
            if (0 < x < img_stamp.shape[0]) and (0 < y < img_stamp.shape[1]):
                print(float(x), float(img_stamp.shape[1] - y), obj["OTYPE"])
                ax.text(x, img_stamp.shape[1] - y, obj["OTYPE"], color='blue')
                ax.plot([x], [img_stamp.shape[1] - y], 'x', color='red')
            plt.xlim(0, img_stamp.shape[0])
            plt.ylim(img_stamp.shape[1], 0)

    # plt.vlines([stamp_margin[0], img_stamp.shape[0] - stamp_margin[0]], 1, img_stamp.shape[1] - 1)
    # plt.hlines([stamp_margin[1], img_stamp.shape[1] - stamp_margin[1]], 1, img_stamp.shape[0] - 1)

    plt.savefig("%s_1%s.png" % (file_prefix, suffix))

    # Figure 2: Image with agressive scaling

    # Load segmentation map file
    if seg_file is not None:
        print("Loading seg_stamp...")
        seg_stamp = fits.getdata(seg_file)
        mask += seg_stamp > 0

    print(fsize)
    fig = plt.figure(2, figsize=fsize, dpi=100)
    plt.clf()
    ax = fig.add_axes([0, 0, 1, 1])
    ax.axis('off')
    aux = np.ma.masked_array(img_stamp, mask=mask)
    if aux.min() < 0:
        aux -= np.min(aux.compressed())
    median = np.median(aux.compressed())
    stdev = np.std(aux.compressed())
    vmin = 0.1 if median <= stdev else median - stdev
    plt.imshow(img_stamp[::-1], cmap=plt.cm.gray_r, vmin=vmin, vmax=(median + stdev), norm=LogNorm())  #
    plt.savefig("%s_2%s.png" % (file_prefix, suffix))

    # Figure 3: Image with a gaussian filter applied on it
    fig = plt.figure(3, figsize=fsize, dpi=100)
    plt.clf()
    ax = fig.add_axes([0, 0, 1, 1])
    ax.axis('off')
    img_stamp[mask] = np.random.normal(loc=median, scale=stdev)

    img_g = ndimage.gaussian_filter(img_stamp, sigma=(3, 3), order=0)
    # img_g /= img_stamp
    median = np.median(img_g)
    stdev = np.std(img_g)
    vmin = 0.1 if median <= stdev else median - stdev
    # plt.imshow(np.ma.masked_array(img_g, mask=mask)[::-1], cmap=plt.cm.gray_r, norm=LogNorm()) #, vmin=vmin, vmax=(median + stdev))
    plt.imshow(np.ma.masked_array(img_g, mask=mask)[::-1], cmap=plt.cm.hot_r, vmin=vmin, vmax=(median + stdev))

    plt.savefig("%s_3%s.png" % (file_prefix, suffix))

    # Figure 4: same as Fig. 3, with SEXtractor
    ## Plot same as Fig. 3:
    fig = plt.figure(4, figsize=fsize, dpi=100)
    plt.clf()
    ax = fig.add_axes([0, 0, 1, 1])
    ax.axis('off')
    plt.imshow(np.ma.masked_array(img_g, mask=mask), origin='lower', cmap=plt.cm.hot_r, vmin=vmin,
               vmax=(median + stdev))

    ## SExtractor
    if sextractor_plot:
        tmp_file = os.path.split(img_file)
        tmp_file = os.path.join(tmp_file[0], "fig4_" + tmp_file[1])
        tmp_wei = tmp_file.replace('fig4', 'fig4wei')
        tmp_cat = tmp_file.replace('fig4_', 'catalog_')
        fits.writeto(tmp_file, img_g)
        kwargs = {'code': 'SExtractor',
                  'config_file': os.path.expanduser('~/workspace/lsb/lsb_project/lsb/config.sex'),
                  'cmd': '/usr/local/bin/sex',
                  'temp_path': '.',
                  'config': {"FILTER_NAME": "/usr/local/share/sextractor/default.conv",
                             "CATALOG_NAME": tmp_cat,
                             "CATALOG_TYPE": "FITS_LDAC",
                             },
                  "params": ["XMIN_IMAGE", "YMIN_IMAGE", "FWHM_IMAGE", "A_IMAGE", "B_IMAGE", "THETA_IMAGE", "FLAGS"]
                  }

        if wei_file is not None:
            wei = np.copy(wei_stamp)
            fits.writeto(tmp_wei, img_g)
            wei[mask] = 0
            kwargs['config'].update({
                'WEIGHT_TYPE': 'MAP_WEIGHT',
                'WEIGHT_IMAGE': wei_file,
            })
        sex = aw.api.Astromatic(**kwargs)
        sex.run(tmp_file)
        catalog = fits.getdata(tmp_cat, ext=2)
        catalog = catalog[catalog["FLAGS"] == 0]
        # remove edge detections
        # pdb.set_trace()
        catalog = catalog[np.bitwise_and(catalog["XMIN_IMAGE"] > stamp_margin[0],
                                         catalog["XMIN_IMAGE"] < img_g.shape[0] - stamp_margin[0])]
        catalog = catalog[np.bitwise_and(catalog["YMIN_IMAGE"] > stamp_margin[1],
                                         catalog["YMIN_IMAGE"] < img_g.shape[1] - stamp_margin[1])]

        ## Plot circles on SEXtractor detections
        mask = (catalog["FWHM_IMAGE"] / 2) >= typical_radius[70]
        circles(catalog["XMIN_IMAGE"][mask], catalog["YMIN_IMAGE"][mask], catalog["FWHM_IMAGE"][mask] / 2,
                edgecolor="blue",
                facecolor="none", alpha=1)
        # ellipses(catalog["XMIN_IMAGE"], catalog["YMIN_IMAGE"], w=catalog["A_IMAGE"], h=catalog["A_IMAGE"], rot=catalog["THETA_IMAGE"], edgecolor="blue", facecolor="none",
        #         alpha=1)
        plt.savefig("%s_4%s.png" % (file_prefix, suffix))
        os.unlink(tmp_file)
        os.unlink(tmp_cat)

    # Figure 5: Color stamp
    if color_stamp is not None:
        fig = plt.figure(5, figsize=fsize, dpi=100)
        plt.clf()
        ax = fig.add_axes([0, 0, 1, 1])
        ax.axis('off')
        plt.imshow(color_stamp, origin='lower')
        plt.savefig("%s_5%s.png" % (file_prefix, suffix))
        # if catalog is not None and plot_catalog_color:
        # for obj in catalog:
        #     # pdb.set_trace()
        #     x, y = skycoord_to_pixel(SkyCoord([[obj["RA"], obj["DEC"]]], unit=('hour', 'deg')), wcs)
        #     x, y = float(x), float(y)
        #     if (0 < x < img_stamp.shape[0]) and (0 < y < img_stamp.shape[1]):
        #         print(float(x), float(img_stamp.shape[1] - y), obj["OTYPE"])
        #         ax.text(x, img_stamp.shape[1] - y, obj["OTYPE"], color='blue')
        #         ax.plot([x], [img_stamp.shape[1] - y], 'x', color='red')


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
                                       color_stamp=color_stamp)

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
