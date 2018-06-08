import os
import pickle

with open(os.path.expanduser("~/data/ngc3115/work/stamps/stamp_list.pkl"), 'r') as fp:
    stamp_data = pickle.load(fp)

stamp_size = [880, 880]

script = "ngc3115\n"
for bin_number in stamp_data:
    if stamp_data[bin_number][3]:
        wcs = stamp_data[bin_number][2]
        ra1, dec1 = wcs.all_pix2world([[0, 0]], 0)[0]
        ra2, dec2 = wcs.all_pix2world([stamp_size], 0)[0]
        # print('draw red box(%.2f, %.2f, %.2f, %.2f, 0, %s)' % (ra1, dec1, abs(-ra2 + ra1), abs(dec2 - dec1), bin_number))
        script += 'draw red box(%.2f, %.2f, %.2f, %.2f)\n' % (ra1, dec1, abs(-ra2 + ra1), abs(dec2 - dec1))
        # script += 'draw white string(%.2f, %.2f, %s)\n' % (ra1 - abs(ra1 - ra2)/2., dec2 + abs(dec2 - dec1)/2., bin_number)
        script += 'draw white string(%.2f, %.2f, %s)\n' % (ra1, dec1, bin_number)

os.unlink('aladin.ajs')
with open('aladin.ajs', 'w') as fp:
    fp.write(script)

