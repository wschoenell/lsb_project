import numpy as np
from astropy import units


def label_pc(z, sizes, cosmo=None):
    """
    Return the proper radius for each size in sizes in pc and z.
    """
    if cosmo is None:
        from astropy.cosmology import WMAP9 as cosmo
    pc_sec = cosmo.kpc_proper_per_arcmin(z).to("pc / arcsec")
    return np.array(sizes) * units.pc / pc_sec