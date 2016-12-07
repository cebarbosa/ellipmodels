#. -*- coding: utf-8 -*-
"""

Created on 07/12/2016

@Author: Carlos Eduardo Barbosa

"""
import os

import numpy as np
import pyfits as pf
from pyraf import iraf

iraf.stsdas(_doprint=False)
iraf.analysis(_doprint=False)
iraf.isophote(_doprint=False)

def run_ellipse():
    wdir = os.path.join(os.getcwd(), "SFCG280")
    os.chdir(wdir)
    imgs = ["rband_SFCG280.wcs.fits"]
    for img in imgs:
        objid = img.split("_")[1].split(".")[0]
        pars = np.loadtxt("{}_apers.dat".format(objid))
        for i,p in enumerate(pars):
            x0, y0, a, b, pa = p
            pa0 = pa + 90
            pa0 = pa0 if pa0 <= 90. else pa0 - 180

            ellip = 1 - b / a
            outtable = "{0}_{1:02d}.tab".format(objid, i+1)
            if not os.path.exists(outtable):
                iraf.ellipse(input=img, output=outtable, x0=x0, y0=y0,
                             ellip0=ellip, pa0=pa0, hcenter="yes", hellip="yes",
                             hpa="yes", interactive="no")
            bmodel = outtable.replace(".tab", "_bmodel.fits")
            if not os.path.exists(bmodel):
                iraf.bmodel(table=outtable, output=bmodel, parent=img)
        data = pf.getdata(img)
        bmodel = np.array([pf.getdata(x) for x in os.listdir(".") if x.endswith("_bmodel.fits")])
        bmodel = bmodel.sum(axis=0)
        resid = data - bmodel
        sigimg = resid / np.sqrt(data)
        h = pf.getheader(img)
        pf.writeto("bmodel.fits", bmodel, h, clobber=True)
        pf.writeto("residuals.fits", resid, h, clobber=True)
        pf.writeto("sigmaimg.fits", sigimg, h, clobber=True)




if __name__ == "__main__":
    run_ellipse()
