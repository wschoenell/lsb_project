 * Mueller et al 2017. Data from 2014A-0624 and 2015A-0616
 * Eigenthaler et al 2018. Data from ???


LFTP download from archive::

    lftp -u anonymous,lftp archive.noao.edu
    cd user_2640
    mget -c *osi* *osw*

Weight Maps
-----------

The CP produces weight maps, or inverse variance maps, for many of the data products. These are derived from error propagation.
The values should be taken with some caution. When the values are zero this is equivalent to bad pixels and will also be reflected in the data quality map.
https://www.noao.edu/noao/staff/fvaldes/CPDocPrelim/PL201_3.html

    Considerations on the size of the UDGs
    --------------------------------------

    Fig 8. of Eigenthaler et al 2018 shows an estimate for the typical sizes of the dwarf galaxies that one should expect to find on data like the DECam.
    The observations are generally limited by the seeing (we use 1.5 arcsec in our case) and the upper limit of detection is of about log(r_eff) ~ 3.5.

    Using the distance to NGC3115 $z = 0.002222$, we estimate that 1.5 arcsec radius corresponds to 70 pc and 1000 pc corresponds to 21.5 arcsec of radius.


File name convention
--------------------

osi: Object Stacked Image
osw: Object Stacked Weight


Color Image
-----------
TODO: add weights and do the weighted-resampling!
a::


    fname="/Users/william/data/ngc3115/c4d_170216_050619_osi_g_v2.fits.fz" # 10:07:38.23 -7:07:07.8
    fout="coadd_g_1.fits"
    fname="/Users/william/data/ngc3115/c4d_170216_061516_osi_r_v2.fits.fz" # 10:07:38.08 -7:07:05.0
    fout="coadd_r_1.fits"

    fname="/Users/william/data/ngc3115/c4d_170217_075805_osi_g_v2.fits.fz" # 10:02:50.00 -8:19:08.0
    fout="coadd_g_2.fits"
    fname="/Users/william/data/ngc3115/c4d_170218_051026_osi_r_v2.fits.fz" # 10:02:50.03 -8:19:04.0
    fout="coadd_r_2.fits"
    for i_hdu in `seq 1 9`
    do
        echo $fname $i_hdu
        w=`echo $fname | sed s/osi/osw/`
        funpack -E $i_hdu -O $fname\_$i_hdu.fits $fname
        funpack -E $i_hdu -O $fname\_$i_hdu.weight.fits $w
    done
    swarp -c /Users/william/workspace/lsb/lsb_project/conf/swarp.conf $fname*[0-9].fits
    mv coadd.fits $fout
    rm c4d_*.fits coadd.weight.fits

    #field1 = coadd_10:07:38.23_-7:07:07.8.png
    #field2 = coadd_10\:02\:50.03_-8\:19\:04.0.png

    # Resampling...
    # Config for field 1
        CENTER         10:07:38.23 -7:07:07.8  # Coordinates of the image center
        PIXELSCALE_TYPE        MANUAL          # MANUAL,FIT,MIN,MAX or MEDIAN
        PIXEL_SCALE            0.263           # Pixel scale
        IMAGE_SIZE             30900,28000     # Image size (0 = AUTOMATIC)
    # Config for field 2
        CENTER         10:02:50.00 -8:19:08.0  # Coordinates of the image center
        PIXELSCALE_TYPE        MANUAL          # MANUAL,FIT,MIN,MAX or MEDIAN
        PIXEL_SCALE            0.263           # Pixel scale
        IMAGE_SIZE             30900,28000     # Image size (0 = AUTOMATIC)


    swarp -c ~/workspace/lsb/lsb_project/conf/swarp_resample_1.conf coadd_g_1.fits coadd_r_1.fits
    swarp -c ~/workspace/lsb/lsb_project/conf/swarp_resample_2.conf coadd_g_2.fits coadd_r_2.fits
    rm *.resamp.weight.fits *[0-9].fits

    # generate color image
    python ~/workspace/willy_scripts/fits_rgb/trilogy.py ~/workspace/lsb/lsb_project/conf/rgb.config

