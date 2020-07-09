import pickle
import sqlite3

from astropy import units
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.wcs import WCS

import numpy as np

input_file = '/Users/william/workspace/lsb/lsb_find/backup/lsb_check_20180901_111156.db'
stamp_dir = '/Users/william/data/ngc3115/work/'
outfile = "matches.csv"
# fp = open(outfile, 'w')
# fp.write("username, stamp_number, x, y, ra, dec, comment\n")

bad_formatted_answers = [300]

conn = sqlite3.connect(input_file)

# Query stamps WCS
c = conn.cursor()
c.execute(
    """
select id, file_prefix from stamps
order by id
""")
stamps_wcs = dict()
for stamp in c.fetchall():
    aux = stamp_dir + stamp[1].replace("c:\\Users\\William Schoenell\\ownCloud\\ngc3115\\work\\\\", "").replace('\\',
                                                                                                                '/')
    # fitsfile = "stamp_%s_img.fits.fz" % aux[-4:]
    wcs = WCS(header=fits.getheader("%s_img.fits.fz" % aux, ext=1))
    stamps_wcs[stamp[0]] = wcs

coord_array = []
out_array = []

for username in ['william', 'apieres', 'Ana', 'crisf', 'basilio', 'Karla']:
    c = conn.cursor()
    c.execute(
        """
select stamps.id, coordinates, stamps.number from answers
  join users
    on users.id = answers.user_id
  join stamps
    on stamps.id = answers.stamp_id
 where answers.coordinates != ''
 and users.name = '%s'
 """ % username)
    # print(len(c.fetchall()))
    for answer in c.fetchall():
        for coordinate in answer[1].split('\n'):
            if ',' in coordinate and answer[0] not in bad_formatted_answers:
                coord = [int(x) if x != "" else None for x in coordinate.split('#')[0].split(',')]
                comment = coordinate.split("#", 1)
                comment = comment[1] if len(comment) > 1 else ""
                coord_wcs = stamps_wcs[answer[0]].all_pix2world([coord], 0)
                out_array.append([username, answer[2], coord[0], coord[1], coord_wcs.ravel()[0], coord_wcs.ravel()[1],
                                  comment.replace(',', ';')])
                coord_array.append([out_array[-1][4], out_array[-1][5]])
                out_str = "%s, %i, %i, %i, %.3f, %.3f, %s\n" % tuple(out_array[-1])
                # fp.write(out_str)
                print(out_str.replace('\n', ''))

# fp.close()


### Pythoninc crossmatch/clustering ####
# grp_array = []
# while len(out_array) > 0:
#     out_array.append(out_array.pop(0))
#     coord_array.pop(0)
#     for val in out_array:
#         distance = SkyCoord()
coord_array = SkyCoord(coord_array, unit=(units.deg, units.deg))
idx1, idx2, sep2d, dist3d = coord_array.search_around_sky(coord_array, seplimit=1 * units.arcmin)

match_list = []
for i, match in enumerate(idx1):
    if (idx1 == match).sum() > 1:
        #[out_array[match]] +
        m = [out_array[x] for x in idx2[idx1 == match]]
        if len(np.unique(np.array(m)[:, 0])) > 2:
            match_list.append(m)

# The final list of matches has many duplicate results. Save only the unique records.
match_list = np.unique(match_list)

with open('match_list.pkl', 'w') as fp:
    pickle.dump(match_list, fp)
