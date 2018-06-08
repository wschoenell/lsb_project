import sys

import numpy as np
import matplotlib.pyplot as plt
from astropy.modeling.functional_models import Sersic2D

try:
    from lsb.util import label_pc
except ImportError:
    sys.path.append('/Users/william/PyCharmProjects/lsb_project/')
    from lsb.util import label_pc

z = 0.002212  # NGC3115 redshift
plate_scale = 0.263  # DECam

typical_radius = {i: label_pc(z, i).value / plate_scale for i in [100, 1000, 5000]}  # in pixels

img_size = [800, 800]
x_0, y_0 = img_size[0] / 2, img_size[1] / 2

debug = False

for size in typical_radius:  # typical_radius[100]
    x, y = np.meshgrid(np.arange(img_size[0]), np.arange(img_size[1]))

    # mod = Sersic2D(amplitude=1, r_eff=typical_radius[size], n=4, x_0=img_size[0] / 2, y_0=img_size[1] / 2, ellip=.5,
    #                theta=-1)

    # Sersic indexes from fig 5 of eigenthaler etal 2018
    mod = Sersic2D(amplitude=1, r_eff=typical_radius[size], n=0.8, x_0=x_0, y_0=y_0, ellip=.5,
                   theta=-1)
    img = mod(x, y)
    # img += np.random.poisson(lam=img, size=img_size)
    img += np.random.poisson(lam=1, size=img_size) * 2
    log_img = np.log10(img)

    signal_slice = int(x_0 - typical_radius[size] / 2.), int(x_0 + typical_radius[size] / 2.), int(
        y_0 - typical_radius[size] / 2.), int(y_0 + typical_radius[size] / 2.)
    noise_slice = [0, 100, 0, 100]
    sn = np.average(img[signal_slice[0]:signal_slice[1], signal_slice[2]:signal_slice[3]]) / np.std(
        img[noise_slice[0]:noise_slice[1], noise_slice[2]:noise_slice[3]])
    print("S/N = %.2f" % sn)

    fsize = img_size[1] / 100., img_size[0] / 100.
    fig = plt.figure(1, figsize=fsize, dpi=100)
    plt.clf()
    ax = fig.add_axes([0, 0, 1, 1])
    ax.axis('off')

    median = np.median(img)
    stdev = np.std(img)

    # plt.imshow(log_img, origin='lower', interpolation='nearest', vmin=-1, vmax=2, cmap=plt.cm.gray_r)
    plt.imshow(img, origin='lower', interpolation='nearest', cmap=plt.cm.gray_r,
               vmin=median - stdev,
               vmax=median + stdev)
    if debug:
        plt.plot([signal_slice[0], signal_slice[1], signal_slice[1], signal_slice[0]],
                 [signal_slice[2], signal_slice[3], signal_slice[2], signal_slice[3]], '.', color='red')
        plt.plot([noise_slice[0], noise_slice[1], noise_slice[1], noise_slice[0]],
                 [noise_slice[2], noise_slice[3], noise_slice[2], noise_slice[3]], '.', color='blue')
    plt.show()
    plt.savefig("sersic_example_%i.png" % size)
