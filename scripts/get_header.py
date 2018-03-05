import os
from astropy.io import fits

data = ""

for dirpath, dirnames, filenames in os.walk(os.path.expanduser("~/data/lsb_work/")):
    for filename in filenames:
        if filename.endswith("fits.fz") and not filename.endswith("weight.fits.fz"):
            hdr = fits.getheader(os.path.join(dirpath, filename), ext=1)
            data += "%s %s %s\n" % (filename, hdr['CRVAL1'], hdr['CRVAL2'])

with open('filelist.txt', 'w') as fp:
    fp.write(data)
