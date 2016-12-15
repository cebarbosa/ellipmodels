#. -*- coding: utf-8 -*-
"""

Created on 07/12/2016

@Author: Carlos Eduardo Barbosa

"""
import os
import shutil

import numpy as np
import pyfits as pf
from pyraf import iraf, iraffunctions
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.io.ascii import SExtractor
from astropy.modeling.models import Ellipse2D
import matplotlib.pyplot as plt

iraf.stsdas(_doprint=False)
iraf.analysis(_doprint=False)
iraf.isophote(_doprint=False)

def run_ellipse_fixed_geom(group):
    """ Run ellipse into a given group. """
    wdir = os.path.join(os.getcwd(), group)
    os.chdir(wdir)
    rband = "rband_{}.wcs.fits".format(group)
    rsub ="rband_{}.wcs.bcksub.fits".format(group)
    catalog = "rband_{}.2.5.cat".format(group)
    cat = SExtractor().read(catalog)
    data = np.genfromtxt("{}_apers.dat".format(group), dtype=None)
    names = []
    idxs = []
    for x0, y0, a, b, pa in data:
        # Finding object in catalog
        idx = np.argmin(np.sqrt((x0 - cat["X_IMAGE"])**2 +
                                (y0 - cat["Y_IMAGE"])**2))
        idxs.append(idx)
        center = SkyCoord(ra=cat[idx]["ALPHA_SKY"] * u.degree,
               dec=cat[idx]["DELTA_SKY"] * u.degree).to_string(style=u"hmsdms")
        # Make identification of object
        name = "".join([x.split(".")[0] for x in center.split()])
        for char in "hdms. ":
            name = name.replace(char, "")
        names.append(name)
        rundir = os.path.join(wdir, name)
        if not os.path.exists(rundir):
            os.mkdir(rundir)
        os.chdir(rundir)
        iraffunctions.chdir(rundir)
        shutil.copy(os.path.join(wdir, rsub), ".")
        shutil.copy(os.path.join(wdir, rband), ".")
        make_mask(rband, cat, avoid_objs=[idx])
        if not os.path.exists(rsub + ".pl"):
            iraf.imcopy(input="mask.fits[type=mask]", output=rsub + ".pl")
        pa0 = pa + 90
        pa0 = pa0 if pa0 <= 90. else pa0 - 180
        ellip = 1 - b / a
        outtable = "fixed_geom.tab"
        if not os.path.exists(outtable):
            iraf.ellipse(input=rsub, output=outtable, x0=x0, y0=y0,
                         ellip0=ellip, pa0=pa0, hcenter="yes", hellip="yes",
                         hpa="yes", interactive="no")
        bmodel = outtable.replace(".tab", "_bmodel.fits")
        if not os.path.exists(bmodel):
            iraf.bmodel(table=outtable, output=bmodel, parent=rsub)
    os.chdir(wdir)
    iraffunctions.chdir(wdir)
    make_mask(rband, cat, avoid_objs=idxs)
    mask = pf.getdata("mask.fits")
    r = pf.getdata(rband)
    sub = pf.getdata(rsub)
    sky = r - sub
    rms = np.sqrt(np.nanmean(sky**2))
    bmodel = np.array([pf.getdata(os.path.join(wdir, gal,
                                "fixed_geom_bmodel.fits")) for gal in names])
    bmodel = bmodel.sum(axis=0)
    resid = sub - bmodel
    sigimg = resid / np.sqrt(r)
    plt.imshow(np.ma.masked_where(mask!=0, sigimg), vmin=-3, vmax=3,
               interpolation="none", origin="bottom")
    plt.show(block=True)
    h = pf.getheader(rband)
    pf.writeto("bmodel.fits", bmodel, h, clobber=True)
    pf.writeto("residuals.fits", resid, h, clobber=True)
    pf.writeto("sigmaimg.fits", sigimg, h, clobber=True)
    return

def make_mask(im, cat, avoid_objs=None, redo=False):
    """ Produces a mask for a given image from a SExtractor catalog. """
    if os.path.exists("mask.fits") and not redo:
        return
    imdata = pf.getdata(im)
    h = pf.getheader(im)
    xis = [np.arange(x) for x in imdata.shape]
    xx, yy = np.meshgrid(*xis)
    mask = np.zeros_like(imdata)
    for i, obj in enumerate(cat):
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
    return






if __name__ == "__main__":
    run_ellipse_fixed_geom("SFCG280")
