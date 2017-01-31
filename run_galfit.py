# -*- coding: utf-8 -*-
"""

Created on 30/01/2017

@Author: Carlos Eduardo Barbosa

"""
import os
from subprocess import call

import numpy as np
import pyfits as pf
import astropy.units as u
from astropy.io.ascii import SExtractor
from astropy.coordinates import SkyCoord
from astropy.stats import sigma_clip
import matplotlib.pyplot as plt

def prepare_input_first_interaction(group):
    """ Prepare input file for the first interaction of galfit"""
    rband = "rband_{}.wcs.fits".format(group)
    rsub ="rband_{}.wcs.bcksub.fits".format(group)
    catalog = "rband_{}.2.5.cat".format(group)
    cat = SExtractor().read(catalog)
    data = np.genfromtxt("emission_pos_{}.dat".format(group), dtype=None)
    ###########################################################################
    # Getting zero point
    zps = np.loadtxt(os.path.join(main_dir, "tables/con_to_AB.dat"), dtype="S")
    zps = dict([(x, float(y)) for x,y in zps])
    ###########################################################################
    # Filling header parameters
    fields = [galfit_guide(),
              "# IMAGE and GALFIT CONTROL PARAMETER",
              "A) {}".format(rband), "B) imgblock.fits", "C) none",
                 "D) none", "E) 1", "F) mask.fits", "G) none"]
    h = pf.getheader(rband)
    xdim = h["NAXIS1"]
    ydim = h["NAXIS2"]
    ps_x = np.sqrt(h['CD1_1'] ** 2. + h['CD1_2'] ** 2.) * 3600
    ps_y = np.sqrt(h['CD2_1'] ** 2. + h['CD2_2'] ** 2.) * 3600
    exptime = h["EXPTIME"]
    ncombine = h["NCOMBINE"]
    ###########################################################################
    # Update gain keyword
    if "GAIN" not in h.keys():
        d = pf.getdata(rband)
        h["GAIN"] = h["CCDSENS"]
        pf.writeto(rband, d, h, clobber=True)
    ###########################################################################
    fields.append("H) 1 {} 1 {} # xmin xmax ymin ymax".format(xdim, ydim))
    fields.append("I) 100 100")
    fields.append("J) {}".format(zps[group] + 25. - 2.5 * np.log10(exptime/ncombine)))
    fields.append("K {:.3f} {:.3f} # Plate scale (dx dy)   [arcsec per pixel]".format(
                   ps_x, ps_y))
    fields.append("O) regular # Display type (regular, curses, both)")
    fields.append("P) 0 # Options: 0=normal run; 1,2=make model/imgblock & quit")
    fields.append("U) 1 # Non-parametric component")
    fields.append("V) 0 # MultiNest")
    fields.append("W) default # output options")
    ###########################################################################
    # Filling body of input with Sersic components for the galaxies
    components = []
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
        pa = cat["THETA_IMAGE"][idx] + 90.
        a = cat["A_IMAGE"][idx]
        b = cat["B_IMAGE"][idx]
        rkron = cat["KRON_RADIUS"][idx]
        flux = cat["FLUX_AUTO"][idx]
        mag = -2.5 * np.log10(flux) + 25 + zps[group]
        comp = []
        comp.append("#" * 78)
        comp.append(" #Galaxy id: {} \n".format(name))
        comp.append(" 0) sersic     # Object type")
        comp.append(" 1) {} 1     # position x [pixel]".format(x0))
        comp.append(" 2) {} 1     # position y [pixel]".format(y0))
        comp.append(" 3) {} 1     # total magnitude in each band".format(mag))
        comp.append(" 4) {} 1     # R_e".format(rkron))
        comp.append(" 5) 2. 1     # Sersic exponent".format())
        comp.append(" 9) {} 1     # axis ration (b/a)".format(b/a))
        comp.append(" 10) {} 1     # position angle (PA)".format(pa))
        comp.append(" Z) 0   # Skip this model in output image? (yes=1, no=0)")
        comp.append("\n")
        components += comp
    ###########################################################################
    # Add simple sky
    sky = ["#" * 78, " # Sky component", " 0) sky"]
    data = pf.getdata(rband)
    skyval =  np.median(sigma_clip(data))
    sky.append(" 1) {:.1f} 1 # sky background [ADU counts]".format(skyval))
    sky.append(" 2) 0.000      0       # dsky/dx (sky gradient in x)")
    sky.append(" 3) 0.000      0       # dsky/dy (sky gradient in y)")
    sky.append(" Z) 0  #  Skip this model in output image?  (yes=1, no=0)")
    ###########################################################################
    # Saving to file
    with open("galfit00.txt", "w") as f:
        f.write("\n".join(fields) + "\n\n")
        f.write("\n".join(components))
        f.write("\n".join(sky))
    return

def galfit_guide():
    return """
================================================================================
# GUIDE TO INPUT FILE FOR GALFITM (a product of the MegaMorph project)
Including multi-band fitting, non-parametric component and MultiNest sampling.
CSL = comma separated list (must not contain any whitespace)
Where several lines are given for the same letter code, these are alternative
examples. The behaviour for multiple lines with the same letter is undefined.
================================================================================
"""



if __name__ == "__main__":
    main_dir = "/home/kadu/Dropbox/halpha_groups"
    data_dir = os.path.join(main_dir, "data")
    groups = sorted([x for x in os.listdir(data_dir)])
    # groups = ["SFCG002"]
    for group in groups:
        print group
        wdir = os.path.join(data_dir, group)
        os.chdir(wdir)
        prepare_input_first_interaction(group)
        call(["{}/bin/galfitm-1.2.1-linux-x86_64".format(main_dir),
              "galfit00.txt"])
        call(["ds9", "-scale", "mode", "99.", "-multiframe", "imgblock00.fits"])
        break