import os
import pickle
import shutil
import sys

import astromatic_wrapper as aw
from astromatic_wrapper.api import AstromaticError
from astropy.io import fits, ascii
from astropy.wcs import WCS
import numpy as np
import matplotlib.pyplot as plt

import matplotlib.image as mpimg
from matplotlib.colors import LogNorm

config = {'min_matches': 2, "work_dir": os.path.expanduser("~/data/ngc3115/work/"), "candidate_size": 150,
          'psf_dir': '/Users/wschoenell/workspace/lsb/lsb_project/scripts/psfex',
          'candidates_dir': '/Users/wschoenell/data/ngc3115/work/candidates'}

# sync cmd: rsync -av --delete /Users/wschoenell/data/ngc3115/work/candidates minerva:public_html/candidates/

with open('match_list.pkl', 'r') as fp:
    match_list = pickle.load(fp)
    match_list = np.unique(match_list)

mask = [len(x) >= config['min_matches'] for x in match_list]

coadd_data = []
for dirpath, dirnames, filenames in os.walk(config["work_dir"]):
    # if "2015A-0616" not in dirpath:
    for filename in filenames:
        if filename.startswith('coadd') and filename.endswith(
                "fz") and not 'weight' in filename and not 'g' in filename and not 'r' in filename:
            # filename = filename.replace('fits.fz', 'g.fits.fz')
            img_file = fits.open(os.path.join(dirpath, filename))
            wei_file = fits.open(os.path.join(config['work_dir'], filename.replace('fits.fz', 'weight.fits.fz')))
            hdr = [x.header for x in img_file]
            coadd_data.append({'filename': filename,
                               'hdr': hdr,
                               'wcs': [WCS(h) for h in hdr],
                               'wei': [ext.data for ext in wei_file],
                               'img': [ext.data for ext in img_file]
                               })

# print('deu')
# sys.exit(0)

for i_match, match in enumerate(np.array(match_list)[mask]):
    # matches shall be for more than min_matches within different users.
    if len(np.unique(np.array(match)[:, 0])) < config['min_matches']:
        continue

    i_bin = match[0][1]
    file_prefix = "%s/stamps/%04i/stamp_%04i" % (config['work_dir'], i_bin, i_bin)

    color_image = "%s_5.png" % file_prefix
    color_stamp = np.flip(mpimg.imread(color_image), 0)
    fig = plt.figure(1)
    fig.clf()
    ax = fig.add_axes([0, 0, 1, 1])
    ax.axis('off')

    # make pixel red
    for m in match:
        color_stamp[m[3], m[2]] = [0, 0, 0, 0]

    ymin, ymax, xmin, xmax = match[0][3] - config["candidate_size"], match[0][3] + config["candidate_size"], match[0][
        2] - config["candidate_size"], match[0][2] + config["candidate_size"]
    xmin = 0 if xmin < 0 else xmin
    ymin = 0 if ymin < 0 else ymin

    # plot paths
    plots_path = "%s/candidate_%03d" % (config["candidates_dir"], i_match)
    if os.path.exists(plots_path):
        shutil.rmtree(plots_path)
    os.mkdir(plots_path)

    c = color_stamp[ymin:ymax, xmin:xmax]
    ax.imshow(c, origin='lower')
    plt.draw_all()
    # plt.show()
    plt.ion()
    plt.savefig("%s/candidate_%03d_1.png" % (plots_path, i_match))

    # plot image
    img_stamp = fits.getdata("%s_img.fits.fz" % file_prefix)[ymin:ymax, xmin:xmax][::-1]
    fig = plt.figure(2, dpi=100)
    plt.clf()
    ax = fig.add_axes([0, 0, 1, 1])
    ax.axis('off')
    aux = np.ma.masked_array(img_stamp)
    if aux.min() < 0:
        aux -= np.min(aux.compressed())
    median = np.median(aux.compressed())
    stdev = np.std(aux.compressed())
    vmin = 0.1 if median <= stdev else median - stdev
    plt.imshow(img_stamp, cmap=plt.cm.gray_r, norm=LogNorm())  # vmin=vmin, vmax=(median + stdev)
    plt.draw_all()
    # plt.show()
    plt.savefig("%s/candidate_%03d_2.png" % (plots_path, i_match))

    # SEXtractor
    x, y, ra, dec = np.array(match[0][2:6])
    aux_coords = np.array([[w.wcs_world2pix(ra, dec, 0) for w in coadd['wcs']] for coadd in coadd_data])

    aux_found = False
    for i in range(aux_coords.shape[0]):
        for j in range(aux_coords.shape[1])[1:]:
            x, y = aux_coords[i][j]
            if not aux_found and \
                    (0 < x < coadd_data[i]['hdr'][j]['NAXIS1']) and \
                    (0 < y < coadd_data[i]['hdr'][j]['NAXIS2']):
                i_file = i
                i_ext = j
                aux_found = True

    # FIXME: hardcoded 30
    if i_match in (30, 31, 32, 33):  # (0, 1, 2, 3):  #
        i_file = 1

    x, y = aux_coords[i_file][i_ext]

    xmin = x - config['candidate_size']
    ymin = y - config['candidate_size']
    xmin = 0 if xmin < 0 else int(xmin)
    ymin = 0 if ymin < 0 else int(ymin)
    xmax = x + config['candidate_size']
    ymax = y + config['candidate_size']
    xmax = int(xmax) if xmax < coadd_data[i_file]['hdr'][i_ext]['NAXIS1'] else coadd_data[i_file]['hdr'][i_ext][
        'NAXIS1']
    ymax = int(ymax) if ymax < coadd_data[i_file]['hdr'][i_ext]['NAXIS2'] else coadd_data[i_file]['hdr'][i_ext][
        'NAXIS2']

    img_stamp = coadd_data[i_file]['img'][i_ext][ymin:ymax, xmin:xmax][::-1]
    wei_stamp = coadd_data[i_file]['wei'][i_ext][ymin:ymax, xmin:xmax][::-1]

    fits.writeto('tmp_img.fits', img_stamp, overwrite=True)
    fits.writeto('tmp_wei.fits', wei_stamp, overwrite=True)

    zp = coadd_data[i_file]['hdr'][0]['MAGZERO'] if "MAGZERO" in coadd_data[i_file]['hdr'][0] else 25
    saturate = coadd_data[i_file]['hdr'][0]["SATURATE"] if "SATURATE" in coadd_data[i_file]['hdr'][0] else 70000
    kwargs = {'code': 'SExtractor',
              'config_file': os.path.expanduser('~/workspace/lsb/lsb_project/lsb/config.sex'),
              'cmd': '/usr/local/bin/sex',
              'temp_path': '.',
              'config': {'WEIGHT_TYPE': 'MAP_WEIGHT',
                         'WEIGHT_IMAGE': 'tmp_wei.fits',
                         "FILTER_NAME": "/usr/local/share/sextractor/default.conv",
                         "CHECKIMAGE_TYPE": "SEGMENTATION",
                         "MAG_ZEROPOINT": str(zp),
                         "CHECKIMAGE_NAME": 'tmp_seg.fits',
                         # "ANALYSIS_THRESH": '0.5',
                         # "DETECT_THRESH": '26,%s' % str(zp),
                         "SATUR_LEVEL": str(saturate),
                         "CLEAN": "Y",
                         # "DEBLEND_MINCONT": "1",   # disables de-blending
                         "DEBLEND_MINCONT": "0.01",
                         "DETECT_MINAREA": "20"
                         },
              "params": ["NUMBER", "ALPHA_J2000", "DELTA_J2000", "XWIN_IMAGE", "YWIN_IMAGE",
                         "FLUX_AUTO", "FLUXERR_AUTO", "MAG_AUTO", "MAGERR_AUTO",
                         "MAG_ISO", "MAGERR_ISO", "ISOAREA_IMAGE",
                         "MU_MAX",
                         "ELLIPTICITY", "FLAGS"]}
    sex = aw.api.Astromatic(**kwargs)
    sex.run('tmp_img.fits')

    try:
        catalog = ascii.SExtractor().read("test.cat")
        # catalog = catalog[catalog["FLAGS"] == 0]
    except IndexError:
        catalog = None

    fig = plt.figure(3, dpi=100)
    plt.clf()
    ax = fig.add_axes([0, 0, 1, 1])
    ax.axis('off')
    plt.imshow(img_stamp, cmap=plt.cm.gray_r, norm=LogNorm())  # vmin=vmin, vmax=(median + stdev)
    plt.plot([x - xmin], [y - ymin], 'o', color='blue')

    if catalog is not None:
        # Draw the magnitudes
        for pt in catalog:
            plt.text(pt['XWIN_IMAGE'], pt['YWIN_IMAGE'],
                     "%.2f\n%.2f" % (pt['MAG_AUTO'], pt['MU_MAX']), color='red')

        seg = fits.getdata("tmp_seg.fits")
        plt.imshow(seg, alpha=0.5, cmap=plt.cm.Paired)
        plt.draw_all()
        plt.show()

        ## Run SEXtractor in dual mode.
        filename_r = os.path.join(config['work_dir'], coadd_data[i_file]["filename"].replace("fits.fz", "r.fits.fz"))
        if not os.path.exists(filename_r):
            print("error", filename_r)
            sys.exit(1)
        file_r = fits.open(filename_r)

        ### GET STAMP ####
        aux_coords = np.array([WCS(ext.header).wcs_world2pix(ra, dec, 0) for ext in file_r])
        aux_found = False
        for j in range(aux_coords.shape[0]):
            x, y = aux_coords[j]
            if not aux_found and \
                    (0 < x < file_r[j].header['NAXIS1']) and \
                    (0 < y < file_r[j].header['NAXIS2']):
                i_ext = j
                aux_found = True

        x, y = aux_coords[i_ext]

        xmin = x - config['candidate_size']
        ymin = y - config['candidate_size']
        xmin = 0 if xmin < 0 else int(xmin)
        ymin = 0 if ymin < 0 else int(ymin)
        xmax = x + config['candidate_size']
        ymax = y + config['candidate_size']
        xmax = int(xmax) if xmax < file_r[i_ext].header['NAXIS1'] else file_r[i_ext].header['NAXIS1']
        ymax = int(ymax) if ymax < file_r[i_ext].header['NAXIS2'] else file_r[i_ext].header['NAXIS2']

        img_stamp = file_r[i_ext].data[ymin:ymax, xmin:xmax][::-1]
        filename_r_wei = filename_r.replace("fits.fz", "weight.fits.fz")
        file_r_wei = fits.open(filename_r_wei)
        wei_stamp = file_r_wei[i_ext].data[ymin:ymax, xmin:xmax][::-1]

        fits.writeto('tmp_r_img.fits', img_stamp, overwrite=True)
        fits.writeto('tmp_r_wei.fits', wei_stamp, overwrite=True)

        ## SAVE FITS STAMP FOR GALFIT ##
        zp = file_r[0].header['MAGZERO'] if "MAGZERO" in file_r[0].header else 25
        saturate = file_r[0].header["SATURATE"] if "SATURATE" in file_r[0].header else 70000

        fits.writeto("%s/candidate_%03d_r_img.fits" % (plots_path, i_match), img_stamp,
                     header=fits.Header({"MAGZERO": zp, "SATURATE": saturate, "EXPTIME": 1}))  #file_r[0].header["EXPTIME"]
        fits.writeto("%s/candidate_%03d_r_wei.fits" % (plots_path, i_match), 1/wei_stamp)  # FIXME: check wei

        ### PSF
        # filename_r_psf = os.path.join(config["psf_dir"],
        #                               "snap_" + os.path.basename(filename_r.replace("fits.fz", "catalog.fits")))
        # fits.writeto(os.path.join(plots_path, "candidate_%03d_r_psf.fits" % i_match),
        #              fits.getdata(filename_r_psf, ext=i_ext))

        ##################

        kwargs["config"].update({'WEIGHT_TYPE': "MAP_WEIGHT,MAP_WEIGHT",
                                 'WEIGHT_IMAGE': "tmp_wei.fits,tmp_r_wei.fits",
                                 "MAG_ZEROPOINT": str(zp),
                                 "CHECKIMAGE_NAME": 'tmp_r_seg.fits',
                                 "SATUR_LEVEL": str(saturate)})

        sex = aw.api.Astromatic(**kwargs)
        sex_err = False
        try:
            sex.run(['tmp_img.fits,tmp_r_img.fits'])
        except AstromaticError:
            print("Problem running SEXtractor")
            sex_err = True

        fig = plt.figure(4, dpi=100)
        plt.clf()
        ax = fig.add_axes([0, 0, 1, 1])
        ax.axis('off')
        plt.imshow(img_stamp, cmap=plt.cm.gray_r, norm=LogNorm())  # vmin=vmin, vmax=(median + stdev)
        # plt.plot([x - xmin], [y - ymin], 'o', color='blue')
        if not sex_err:
            try:
                catalog = ascii.SExtractor().read("test.cat")
                # catalog = catalog[catalog["FLAGS"] == 0]
            except IndexError:
                catalog = None

            if catalog is not None:
                # Draw the magnitudes
                for pt in catalog:
                    plt.text(pt['XWIN_IMAGE'], pt['YWIN_IMAGE'],
                             "%.2f\n%.2f" % (pt['MAG_AUTO'], pt['MU_MAX']), color='red')
                # Draw segmentation map
                seg = fits.getdata("tmp_seg.fits")
                plt.imshow(seg, alpha=0.5, cmap=plt.cm.Paired)

    plt.figure(3)
    plt.savefig("%s/candidate_%03d_3.png" % (plots_path, i_match))
    plt.figure(4)
    plt.savefig("%s/candidate_%03d_4.png" % (plots_path, i_match))

    # save stamps for GALFIT

    ### FOR G FILTER ####
    filename_g = os.path.join(config['work_dir'], coadd_data[i_file]["filename"].replace("fits.fz", "r.fits.fz"))
    file_g = fits.open(filename_g)
    aux_coords = np.array([WCS(ext.header).wcs_world2pix(ra, dec, 0) for ext in file_g])
    aux_found = False
    for j in range(aux_coords.shape[0]):
        x, y = aux_coords[j]
        if not aux_found and \
                (0 < x < file_g[j].header['NAXIS1']) and \
                (0 < y < file_g[j].header['NAXIS2']):
            i_ext = j
            aux_found = True

    x, y = aux_coords[i_ext]

    xmin = x - config['candidate_size']
    ymin = y - config['candidate_size']
    xmin = 0 if xmin < 0 else int(xmin)
    ymin = 0 if ymin < 0 else int(ymin)
    xmax = x + config['candidate_size']
    ymax = y + config['candidate_size']
    xmax = int(xmax) if xmax < file_g[i_ext].header['NAXIS1'] else file_g[i_ext].header['NAXIS1']
    ymax = int(ymax) if ymax < file_g[i_ext].header['NAXIS2'] else file_g[i_ext].header['NAXIS2']

    img_stamp_g = file_g[i_ext].data[ymin:ymax, xmin:xmax][::-1]
    filename_g_wei = filename_g.replace("fits.fz", "weight.fits.fz")
    file_g_wei = fits.open(filename_g_wei)
    wei_stamp_g = file_g_wei[i_ext].data[ymin:ymax, xmin:xmax][::-1]
    zp = file_g[0].header['MAGZERO'] if "MAGZERO" in file_g[0].header else 25
    saturate = file_g[0].header["SATURATE"] if "SATURATE" in file_g[0].header else 70000

    fits.writeto("%s/candidate_%03d_g_img.fits" % (plots_path, i_match), img_stamp_g,
                 header=fits.Header({"MAGZERO": zp, "SATURATE": saturate, "EXPTIME": 1}))  # file_g[0].header["EXPTIME"]
    fits.writeto("%s/candidate_%03d_g_wei.fits" % (plots_path, i_match), 1/wei_stamp_g)  # FIXME: check wei

    ### PSF
    # filename_g_psf = os.path.join(config["psf_dir"],
    #                               "snap_" + os.path.basename(filename_g.replace("fits.fz", "catalog.fits")))
    # fits.writeto(os.path.join(plots_path, "candidate_%03d_g_psf.fits" % i_match),
    #              fits.getdata(filename_g_psf, ext=i_ext))

    print(match, coadd_data[i_file]["filename"])
    # raw_input("next...")

    # break

    #
