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

File name convention
--------------------

osi: Object Stacked Image
osw: Object Stacked Weight