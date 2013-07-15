#!/usr/bin/env python

"""
Extract flux calib vectors for the first exposure from a set of spec files.
These are interpolated onto the coadd dloglam=1e-4 wavelength grid.

Writes a fits file with two image HDUs:
    HDU 0 : B calibration vectors [nfibers, nflux]
    HDU 1 : R calibration vectors [nfibers, nflux]

Stephen Bailey, LBL
Summer 2012
"""

import sys
import os
import numpy as N
from glob import glob
import fitsio

import optparse

parser = optparse.OptionParser(usage = "%prog [options]")
parser.add_option("-i", "--indir", type="string",  help="input data directory")
parser.add_option("-o", "--output", type="string",  help="output calib file")
opts, args = parser.parse_args()

#- Min and max ranges for loglam grid
bmin, bmax = (3.5500, 3.8000)
rmin, rmax = (3.7500, 4.0200)

#- log10(lambda) ranges for B and R
dll = 1e-4
bll = N.arange(bmin, bmax+dll/2, dll)
rll = N.arange(rmin, rmax+dll/2, dll)

nfibers = 1000  #- Number of BOSS fibers
calib_b = N.zeros( (nfibers, len(bll)) )
calib_r = N.zeros( (nfibers, len(rll)) )

for specfile in glob(opts.indir + 'spec*.fits'):
    hdr = fitsio.read_header(specfile, 0)
    fiber = hdr['FIBERID']
    print fiber
    
    for expid in range(1, 100):
        expkey = 'EXPID%02d' % expid
        if hdr[expkey].startswith('r'):
            break
    else:
        print >> sys.stderr, "ERROR: no R-exposures found"
        sys.exit(1)
    
    specb = fitsio.read(specfile, 4).view(N.recarray)
    specr = fitsio.read(specfile, hdr[expkey]).view(N.recarray)
    
    #--- DEBUG ---
    # import IPython
    # IPython.embed()
    # sys.exit(1)
    #--- DEBUG ---
    
    calib_b[fiber-1] = N.interp(bll, specb.loglam, specb.calib, left=0.0, right=0.0)
    calib_r[fiber-1] = N.interp(rll, specr.loglam, specr.calib, left=0.0, right=0.0)

#- Write output file
fitsio.write(opts.output, calib_b, clobber=True)    
fitsio.write(opts.output, calib_r)    

#- Update keywords
fx = fitsio.FITS(opts.output, 'rw')
fx[0].write_key('CRVAL1', bmin, 'log10(lambda) of first pixel')
fx[0].write_key('CDELT1', dll, 'delta(log10(lambda)) of wavelength grid')
fx[0].write_key('COMMENT', "Calibration vectors for B-channel")
fx[1].write_key('CRVAL1', rmin, 'log10(lambda) of first pixel')
fx[1].write_key('CDELT1', dll, 'delta(log10(lambda)) of wavelength grid')
fx[1].write_key('COMMENT', "Calibration vectors for R-channel")

