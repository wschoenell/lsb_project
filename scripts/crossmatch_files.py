from __future__ import print_function
import os

import sys
from astropy.coordinates import SkyCoord
from astropy.io import fits
import numpy as np

from astropy import units as u

overwrite_coadd = False

# files_dir = os.path.expanduser("~/data/lsb/data_example/2015A-0616")
# files_dir = os.path.expanduser("~/data/lsb/data_example")
# object = None

files_dir = os.path.expanduser("/home/william/data/ngc3115")
workdir = "~/data/ngc3115/work/"
object = 'NGC 3115'

pointings = []

for dirpath, dirnames, filenames in os.walk(files_dir):
    # if "2015A-0616" not in dirpath:
    for filename in filenames:
        if 'osi' in filename and filename.endswith("fz"):
            f = os.path.join(dirpath, filename)
            hdr = fits.getheader(f)
            print("@@>", hdr["OBJECT"])
            if object is None or hdr["OBJECT"] == object:
                pointings.append([os.path.join(dirpath, filename), hdr["RA"], hdr["DEC"]])
                print(" ".join(pointings[-1]))

pointings = np.array(pointings)

np.savetxt(os.path.expanduser("%s/input_list.txt" % workdir), pointings, fmt="%s\t%s\t%s")

catalog = SkyCoord(["%s %s" % (p[1], p[2]) for p in pointings], unit=(u.hourangle, u.deg))

# Crossmatch catalog to itself
idxsearcharound, idxself, sep2d, dist3d = catalog.search_around_sky(catalog, seplimit=0.1 * u.deg)
## Remove the self matches
msk = idxsearcharound != idxself
idxsearcharound, idxself, sep2d, dist3d = idxsearcharound[msk], idxself[msk], sep2d[msk], dist3d[msk]
## Remove duplicates and join in groups - TODO: this can be more efficiently!
aux_list = list(np.unique(idxsearcharound))
final_match = []
for i, idx in enumerate(aux_list):
    i_match = np.argwhere(idxsearcharound == idx)[:, 0]
    final_match.append([idx] + list(idxself[i_match]))
final_match = [list(np.sort(m)) for m in final_match]
final_match = list(np.unique(final_match, axis=0))
# final_match = list(np.unique(final_match))

# Run SWARP
for m in final_match:

    # Co-added file names
    coadd_file = os.path.expanduser("%s/coadd_%s_%s.fits" % tuple(np.append([workdir], pointings[m[-1]][1:3])))
    coadd_weifile = os.path.expanduser(
        "%s/coadd_%s_%s.weight.fits" % tuple(np.append([workdir], pointings[m[-1]][1:3])))

    # if co-added already exists, skip...
    if os.path.exists(coadd_file + ".fz") and os.path.exists(coadd_weifile + ".fz") and not overwrite_coadd:
        print("Files %s and %s already exists. Skipping coadd." % (coadd_file, coadd_weifile))
    else:
        for i_ext in range(1, 10):
            for img_infile in [pointings[i_m][0] for i_m in m]:
                script = ""
                # Unpack images
                img_file = img_infile.replace("fits.fz", "%i.fits" % i_ext)
                script += "funpack -E %i -O %s %s \n" % (i_ext, img_file, img_infile)

                # Unpack weights
                wei_infile = img_infile.replace("osi", "osw")
                wei_file = "%s.%i.weight.fits" % (img_infile.replace(".fits.fz", ""), i_ext)
                script += "funpack -E %i -O %s %s \n" % (i_ext, wei_file, wei_infile)

                os.system(script)

                for f in img_infile, img_file, wei_infile, wei_file:
                    if not os.path.exists(f):
                        print("File not found", f)
                        sys.exit(1)

        # Run SWarp
        script = "swarp %s -c $HOME/workspace/lsb_project/conf/swarp.conf\n" % (" ".join(
            [" ".join([pointings[i_m][0].replace("fits.fz", "%i.fits" % i_ext) for i_ext in range(1, 10)]) for i_m in
             m]))

        os.system(script)

        for f in "coadd.fits", "coadd.weight.fits":
            if not os.path.exists(f):
                print("File not found", f)
                sys.exit(1)

        # Rename output files
        script = "mv coadd.fits %s\n " % coadd_file
        script += "mv coadd.weight.fits %s\n" % coadd_weifile

        os.system(script)

        for f in coadd_file, coadd_weifile:
            if not os.path.exists(f):
                print("File not found", f)
                sys.exit(1)
            else:
                os.system("fpack %s" % f)
                os.unlink(f)
                print("File %s fpacked and deleted" % f)

        # Remove unpacked files
        for i_ext in range(1, 10):
            for img_infile in [pointings[i_m][0] for i_m in m]:
                script = ""
                # images
                img_file = img_infile.replace("fits.fz", "%i.fits" % i_ext)
                if os.path.exists(img_file):
                    os.unlink(img_file)
                    print("deleted %s" % img_file)

                # weights
                wei_infile = img_infile.replace("osi", "osw")
                wei_file = "%s.%i.weight.fits" % (img_infile.replace(".fits.fz", ""), i_ext)
                if os.path.exists(wei_file):
                    os.unlink(wei_file)
                    print("deleted %s" % wei_file)

# with open("tmp.sh", 'w') as fp:
#     fp.write(swarp_script)
