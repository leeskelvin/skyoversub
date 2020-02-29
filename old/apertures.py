#!/usr/bin/env python3

# imports
import os
import sys
import warnings
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import astropy.stats as stats
from collections import Counter as count
import pandas as pd
import pickle

# LSST stack imports
import lsst.daf.persistence as dafPersist
import lsst.afw.display as afwDisplay
from lsst.ip.isr import IsrTask
import lsst.afw.detection as afwDetection

# setup
butler = dafPersist.Butler("/datasets/hsc/repo/rerun/RC/w_2019_38/DM-21386/")
trac = 9813
pat = "2,2"
filt = "HSC-R"

# input butler data
dataid = {'tract': trac, 'patch': pat, 'filter': filt}
coadd = butler.get('deepCoadd', dataId=dataid)
wcs = coadd.getWcs()
table = butler.get('deepCoadd_meas', dataId=dataid)
skyobjects = table[table['merge_peak_sky']]

# # launch ds9
# !~/software/ds9/ds9 -mode region -regions shape projection -cmap sls -scale log -scale mode 99.5 -wcs degrees &


def makemaskdict(desired, current):
    keys = current.keys()
    output = {key: 'IGNORE' for key in keys}
    for key in keys:
        if key in desired.keys():
            output[key] = desired[key]
    return output


# display data
display1 = afwDisplay.getDisplay(backend='ds9', frame=1)
display2 = afwDisplay.getDisplay(backend='ds9', frame=2)
mpdict = {"DETECTED": "RED", "DETECTED_NEGATIVE": "GREEN", "BAD": "MAGENTA", "NO_DATA": "YELLOW"}
mpdictfull = makemaskdict(mpdict, coadd.mask.getMaskPlaneDict())
display1.setMaskPlaneColor(mpdictfull)
display1.mtv(coadd)
display1.setMaskTransparency(0)
with display1.Buffering():
    for s in skyobjects:
        # display1.dot("○", s.getX(), s.getY(), size=5, ctype='cyan')
        display1.dot("o", s.getX(), s.getY(), size=8, ctype='yellow')
        display1.dot("*", s.getX(), s.getY(), size=100, ctype='yellow')

# in config.py:
# config.measurement.plugins['base_CircularApertureFlux'].radii=[3.0, 4.5, 6.0, 9.0, 12.0, 17.0, 25.0, 35.0, 50.0, 70.0, 250.0, 500.0]
# on command line:
# measureCoaddSources.py /datasets/hsc/repo --rerun RC/w_2019_38/DM-21386:private/lskelvin/apertures_ntrue --id tract=9813 patch=2,2 filter=HSC-R --configfile config.py -c measurement.doReplaceWithNoise=True
# measureCoaddSources.py /datasets/hsc/repo --rerun RC/w_2019_38/DM-21386:private/lskelvin/apertures_nfalse --id tract=9813 patch=2,2 filter=HSC-R --configfile config.py -c measurement.doReplaceWithNoise=False


def getradflux(butler, dataid, radius):
    subtable = butler.get('deepCoadd_meas', dataId=dataid)
    subskyobjects = subtable[subtable['merge_peak_sky']]
    flux = subskyobjects['base_CircularApertureFlux_{}_0_instFlux'.format(radius)]
    flag = ~subskyobjects['base_CircularApertureFlux_{}_0_flag'.format(radius)]
    return flux, flag


# output butlers
nfalse = dafPersist.Butler("/project/lskelvin/hsc-rerun/apertures_nfalse")
ntrue = dafPersist.Butler("/project/lskelvin/hsc-rerun/apertures_ntrue")

# extract values
r500flux, r500flag = getradflux(ntrue, dataid, 500)
r9flux, r9flag = getradflux(ntrue, dataid, 9)



# output catalogues
radii = [9, 500]
for rad in radii:
    # read aperture fluxes
    flux_don, flag_don = getradflux(butler_donoise, dataid, rad)
    flux_dontn, flag_dontn = getradflux(butler_dontnoise, dataid, rad)
    # keep only sky objects that were kept in both cases
    joint_flag = flag_don*flag_dontn
    # Lee plot
    plt.plot(flux_dontn[joint_flag], flux_don[joint_flag], 'o')
    plt.xlabel('Aperture flux, keeping sources')
    plt.ylabel('"Sky"')
    plt.title('Aperture radius: {} ({} sky objects survived)'.format(rad, np.sum(joint_flag)))
    plt.show()
    plt.close()

