#. -*- coding: utf-8 -*-
"""

Created on 07/12/2016

@Author: Carlos Eduardo Barbosa

"""
import os

import numpy as np
import pyfits as pf
from pyraf import iraf
from astropy.io.ascii import SExtractor
from astropy.modeling.models import Ellipse2D
import matplotlib.pyplot as plt

iraf.stsdas(_doprint=False)
iraf.analysis(_doprint=False)
iraf.isophote(_doprint=False)

def run_ellipse():
    wdir = os.path.join(os.getcwd(), "SFCG280")
    os.chdir(wdir)
    imgs = ["rband_SFCG280.wcs.fits"]
    for img in imgs:
        objid = img.split("_")[1].split(".")[0]
        cat = img.replace(".wcs.fits", ".2.5.cat")
        pars = np.loadtxt("{}_apers.dat".format(objid))
        table = SExtractor().read(cat)
        avoid = []
        for i, obj in enumerate(table):
            r = np.sqrt((pars[:,0] - obj["X_IMAGE"])**2 +
                        (pars[:,1] - obj["Y_IMAGE"])**2)
            if r.min() == 0.:
                avoid.append(i)
        # make_mask(img, cat, avoid_objs=avoid)
        iraf.imcopy(input="mask.fits[type=mask]", output=img + ".pl")
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
        bmodel = np.array([pf.getdata(x) for x in os.listdir(".")
                           if x.endswith("_bmodel.fits")])
        bmodel = bmodel.sum(axis=0)
        resid = data - bmodel
        sigimg = resid / np.sqrt(data)
        h = pf.getheader(img)
        pf.writeto("bmodel.fits", bmodel, h, clobber=True)
        pf.writeto("residuals.fits", resid, h, clobber=True)
        pf.writeto("sigmaimg.fits", sigimg, h, clobber=True)

def make_mask(im, cat, avoid_objs=None):
    """ Produces a mask for a given image from a SExtractor catalog. """
    table = SExtractor().read(cat)
    imdata = pf.getdata(im)
    h = pf.getheader(im)
    xis = [np.arange(x) for x in imdata.shape]
    xx, yy = np.meshgrid(*xis)
    mask = np.zeros_like(imdata)
    for i, obj in enumerate(table):
        if i in avoid_objs:
            continue
        x0, y0, a, b, theta, rkron= obj["X_IMAGE"], obj["Y_IMAGE"], \
                                    obj["A_IMAGE"], obj["B_IMAGE"], \
                                    obj["THETA_IMAGE"], obj["KRON_RADIUS"]
        e = Ellipse2D(amplitude=1., x_0=x0, y_0=y0, a=rkron * a,
                      b=rkron * b, theta=np.deg2rad(theta))
        mask += e(xx, yy)
    hdu = pf.PrimaryHDU(data=mask, header=h)
    hdu.writeto("mask.fits", clobber=True)






if __name__ == "__main__":
    run_ellipse()
