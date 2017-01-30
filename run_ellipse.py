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
from astropy.stats import sigma_clip
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse

iraf.stsdas(_doprint=False)
iraf.analysis(_doprint=False)
iraf.isophote(_doprint=False)

def run_ellipse_fixed_geom(group, plot=True):
    """ Run ellipse into a given group. """
    wdir = os.path.join(data_dir, group)
    os.chdir(wdir)
    rband = "rband_{}.wcs.fits".format(group)
    rsub ="rband_{}.wcs.bcksub.fits".format(group)
    catalog = "rband_{}.2.5.cat".format(group)
    cat = SExtractor().read(catalog)
    data = np.genfromtxt("emission_pos_{}.dat".format(group), dtype=None)
    names = []
    idxs_all = []
    make_mask(rsub, cat, outfile="segmentation.fits", redo=False)
    pars = []
    for x0, y0, ra0, dec0 in data:
        # Finding object in catalog
        r = np.sqrt((x0 - cat["X_IMAGE"])**2 + (y0 - cat["Y_IMAGE"])**2)
        idx = np.argmin(r)
        center = SkyCoord(ra0, dec0, unit=(u.hourangle,
                                      u.deg)).to_string(style=u"hmsdms")
        # Make identification of object
        name = "".join([x.split(".")[0] for x in center.split()])
        for char in "hdms. ":
            name = name.replace(char, "")
        names.append(name)
        print name, x0, y0
        rundir = os.path.join(wdir, name)
        if not os.path.exists(rundir):
            os.mkdir(rundir)
        os.chdir(rundir)
        iraffunctions.chdir(rundir)
        shutil.copy(os.path.join(wdir, rsub), ".")
        shutil.copy(os.path.join(wdir, rband), ".")
        idxs = [idx]
        extra = "extra_masks.dat"
        if os.path.exists(extra):
            eidx =  np.atleast_1d(np.loadtxt(extra, dtype=int)).tolist()
            idxs += eidx
        idxs_all += idxs
        make_mask(rsub, cat, avoid_objs=idxs, redo=False)
        if not os.path.exists(rsub + ".pl"):
            iraf.imcopy(input="mask.fits[type=mask]", output=rsub + ".pl")
        pa = cat["THETA_IMAGE"][idx]
        a = cat["A_IMAGE"][idx]
        b = cat["B_IMAGE"][idx]
        rkron = cat["KRON_RADIUS"][idx]
        pa0 = pa + 90
        pa0 = pa0 if pa0 <= 90. else pa0 - 180
        ellip = np.maximum(0.05, 1 - b / a)
        outtable = "fixed_geom.tab"
        if not os.path.exists(outtable):
            iraf.ellipse(input=rsub, output=outtable, x0=x0, y0=y0,
                         ellip0=ellip, pa0=pa0, hcenter="yes", hellip="yes",
                         hpa="yes", interactive="no", sma0=rkron * a)
        bmodel = outtable.replace(".tab", "_bmodel.fits")
        if not os.path.exists(bmodel):
            iraf.bmodel(table=outtable, output=bmodel, parent=rsub)
            # Saving geometrical parameters for plotting
        pars.append([x0, y0, a, b, pa])
    ###########################################################################
    # Using input from ellipse runs to produce models and sigma images
    os.chdir(wdir)
    iraffunctions.chdir(wdir)
    idxs = np.unique(idxs_all)
    rimg = pf.getdata(rband)
    ###########################################################################
    # Reading image header
    h = pf.getheader(rband)
    ncombine = h["NCOMBINE"]
    gain = h["CCDSENS"]
    cd11 = h["CD1_1"]
    cd21 = h["CD2_1"]
    cd12 = h["CD1_2"]
    cd22 = h["CD2_2"]
    ###########################################################################
    # Calculating rotation angle
    m = np.matrix([[cd11, cd12], [cd21, cd22]])
    rad = 57.2957795131
    det = np.linalg.det(m)
    cdelt1 = np.sign(det) * np.sqrt(cd11 ** 2 + cd21 ** 2)
    cdelt2 = np.sqrt(cd12 ** 2 + cd22 ** 2);
    rot = np.arctan2(-cd21 /rad, np.sign(cdelt1) * cd11 / rad) * rad
    rot2 = np.arctan2(np.sign(cdelt1) * cd12 /rad,  cd22 /rad) * rad
    ###########################################################################
    # Loading images
    make_mask(rband, cat, avoid_objs=idxs, redo=False)
    mask =  pf.getdata("mask.fits")
    mask = np.clip(mask, 0, 1)
    r = pf.getdata(rband)
    r0 = pf.getdata(rsub) # Sky-subtracted image
    ###########################################################################
    # Producing sigma image
    shot_noise = np.sqrt(r * ncombine / gain)
    sky_noise = np.std(sigma_clip(r0))
    sigimg = np.sqrt(shot_noise**2 + sky_noise**2)
    ###########################################################################
    # Making bmodel image
    bmodel = np.array([pf.getdata(os.path.join(wdir, gal,
                                "fixed_geom_bmodel.fits")) for gal in names])
    bmodel = bmodel.sum(axis=0)
    resid = r0 - bmodel
    sigres = resid / sigimg
    ###########################################################################
    # Plotting residuals
    if plot:
        ax = plt.subplot(111)
        ax.minorticks_on()
        radius = [1, 4]
        im = ax.imshow(np.ma.array(sigres, mask=mask), origin="bottom",
                   vmin=-1.5, vmax=1.5, cmap="Spectral", interpolation="none")
        for p in pars:
            for r in radius:
                ell = Ellipse(xy=[p[0], p[1]], width= r * p[2],
                              height= r * p[3], angle=p[4], linestyle="--",
                              facecolor="none", lw=1.2, edgecolor="k")
                ax.add_artist(ell)
        plt.colorbar(im)
        plt.title(group)
        plt.show(block=True)
    ###########################################################################
    # Saving images
    pf.writeto("bmodel.fits", bmodel, h, clobber=True)
    pf.writeto("residuals.fits", resid, h, clobber=True)
    pf.writeto("sigmaimg.fits", sigimg, h, clobber=True)
    pf.writeto("sigres.fits", sigres, h, clobber=True)
    return

def make_mask(im, cat, avoid_objs=None, redo=False, outfile="mask.fits"):
    """ Produces a mask for a given image from a SExtractor catalog. """
    if os.path.exists(outfile) and not redo:
        return
    if avoid_objs is None:
        avoid_objs = []
    imdata = pf.getdata(im)
    h = pf.getheader(im)
    xis = [np.arange(x) for x in imdata.shape[::-1]]
    xx, yy = np.meshgrid(*xis)
    mask = np.zeros_like(xx, dtype="float64")
    for i, obj in enumerate(cat):
        if i in avoid_objs:
            continue
        x0, y0, a, b, theta, rkron= obj["X_IMAGE"], obj["Y_IMAGE"], \
                                    obj["A_IMAGE"], obj["B_IMAGE"], \
                                    obj["THETA_IMAGE"], obj["KRON_RADIUS"]
        e = Ellipse2D(amplitude=i, x_0=x0, y_0=y0, a=rkron * a,
                      b=rkron * b, theta=np.deg2rad(theta))
        mask += e(xx, yy)
    hdu = pf.PrimaryHDU(data=mask, header=h)
    hdu.writeto(outfile, clobber=True)
    return






if __name__ == "__main__":
    data_dir = "/home/kadu/Dropbox/halpha_groups/data"
    groups = sorted([x for x in os.listdir(data_dir)])
    # groups = ["SFCG280"]
    for group in groups:
        print group
        run_ellipse_fixed_geom(group)
