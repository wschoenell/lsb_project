import os
import astromatic_wrapper as aw
from astropy.coordinates import Angle, SkyCoord
from astropy.io import fits
from astropy import units as u
import matplotlib.pyplot as plt

from astroquery.vizier import Vizier

import numpy as np
from scipy.optimize import least_squares

Vizier.ROW_LIMIT = -1

filenames_img = ["/Users/william/data/ngc3115/c4d_170216_050619_osi_g_v2.fits",
                 "/Users/william/data/ngc3115/c4d_170216_061516_osi_r_v2.fits",
                 "/Users/william/data/ngc3115/c4d_170217_075805_osi_g_v2.fits",
                 "/Users/william/data/ngc3115/c4d_170218_051026_osi_r_v2.fits",
                 ]

colors = ['blue', 'red', 'magenta', 'black']

i_color = 0
plt.clf()
for fname in filenames_img:
    hdr = fits.getheader(fname, ext=0)
    filename_catalog = "catalog_%s" % os.path.basename(fname)
    filename_checkimg = "seg_%s" % os.path.basename(fname)
    kwargs = {'code': 'SExtractor',
              'config_file': os.path.expanduser('~/workspace/lsb/lsb_project/lsb/config.sex'),
              'cmd': '/usr/local/bin/sex',
              'temp_path': '.',
              'config': {"FILTER_NAME": "/usr/local/share/sextractor/default.conv",
                         "CATALOG_NAME": filename_catalog,
                         "CATALOG_TYPE": "FITS_LDAC",
                         'WEIGHT_TYPE': 'MAP_WEIGHT',
                         'WEIGHT_IMAGE': fname.replace("osi", "osw"),
                         # 'MAG_ZEROPOINT': float(hdr['MAGZPT']) + 2.5 * np.log10(hdr["EXPTIME"])  # MAGZERO
                         'MAG_ZEROPOINT': str(hdr['MAGZERO']),
                         'CHECKIMAGE_TYPE': "SEGMENTATION",
                         'CHECKIMAGE_NAME': filename_checkimg
                         },
              "params": ["XMIN_IMAGE", "YMIN_IMAGE", "ALPHA_J2000", "DELTA_J2000",
                         "FWHM_IMAGE", "A_IMAGE", "B_IMAGE", "THETA_IMAGE", "FLAGS",
                         "MAG_AUTO", "CLASS_STAR"],
              }

    if not os.path.exists(filename_catalog):
        sex = aw.api.Astromatic(**kwargs)
        sex.run(fname)

    catalog = fits.getdata(filename_catalog, ext=2)
    # catalog = catalog[catalog["FLAGS"] == 0]
    # stars = catalog[catalog["CLASS_STAR"] > 0.9]
    stars = catalog[catalog["MAG_AUTO"] < 90]

    hdr = fits.getheader(fname, ext=1)

    result = Vizier.query_region(SkyCoord(ra=hdr["CRVAL1"], dec=hdr["CRVAL2"], unit=(u.deg, u.deg)),
                                 radius=Angle(3.0, "deg"),
                                 catalog='II/336/apass9')
    values = result.values()[0]
    apass_coordinates = SkyCoord(values["RAJ2000"], values["DEJ2000"], unit=(u.deg, u.deg))
    stars_coordinates = SkyCoord(stars["ALPHA_J2000"], stars["DELTA_J2000"], unit=(u.deg, u.deg))

    idx, sep2d, _ = apass_coordinates.match_to_catalog_sky(stars_coordinates)

    if '_g_' in filename_catalog:
        magname = 'g_mag'
    elif '_r_' in filename_catalog:
        magname = 'r_mag'

    aux_label = fname.split('/')[-1].split('_v')[0].replace('_osi', '')
    mask = np.bitwise_and(sep2d.value < 0.0005, values['e_' + magname] < 0.1)
    plt.plot(values[magname][mask], stars[idx]['MAG_AUTO'][mask], '.', label=aux_label,
             color=colors[i_color])
    # plt.errorbar(values[magname][sep2d.value < 0.0003], stars[idx]['MAG_AUTO'][sep2d.value < 0.0003], )

    # mask = np.bitwise_and(mask, values[magname] > -17)

    ####
    # def residual(p, position, fwhm):
    #     model = [p[0] * x ** 2 + p[1] * x + p[2] for x in position]
    #     return np.array(fwhm) - np.array(model)
    # A, B, C = np.polyfit(values[magname][mask], stars[idx]['MAG_AUTO'][mask], 2)
    # A=0
    # robust_fit = least_squares(residual, (A, B, C), args=(values[magname][mask], stars[idx]['MAG_AUTO'][mask]), loss='soft_l1', f_scale=0.5)
    # poly = robust_fit.x
    ####
    poly = np.polyfit(values[magname][mask], stars[idx]['MAG_AUTO'][mask], 1)
    x = np.arange(plt.xlim()[0], plt.xlim()[1])
    y = np.polyval(poly, x)
    # plt.plot(x, y, color=colors[i_color])
    plt.plot(x, x, color=colors[i_color])
    # plt.plot([14], [18], color=colors[i_color])
    i_color += 1

plt.legend()
plt.xlim(12, 19)
plt.ylim(12, 19)
plt.xlabel("APASS mag")
plt.ylabel("SExtractor MAG\_AUTO")
plt.savefig("test.png")