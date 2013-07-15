#!/usr/bin/env python

"""
Squash the x-spacing of spectral traces to mimic MaNGA fiber spacing.
For convenience in reformating the BOSS PIX-PSF format, increase the
number of fibers per bundle while keeping a constant 25 bundles.

Examples:
squash_psf.py -i spBasisPSF-PIX-r1-00131554.fits -o MaNGA-PSF-r1.fits
squash_psf.py -i spBasisPSF-PIX-b1-00131554.fits -o MaNGA-PSF-b1.fits

Stephen Bailey
Summer 2012
"""

import sys
import os
import numpy as N
import pyfits

import optparse
parser = optparse.OptionParser(usage = "%prog [options]")
parser.add_option("-i", "--input", type="string",  help="input BOSS PSF")
parser.add_option("-o", "--output", type="string",  help="output MaNGA PSF")
parser.add_option("-n", "--fibers_per_bundle", type="int", default=30,
    help="Number of fibers per bundle[%default]")

opts, args = parser.parse_args()
if opts.input is None or opts.output is None:
    print >> sys.stderr, "ERROR: you must specify both input and output"
    parser.print_help(sys.stderr)
    sys.exit(1)

nbundles = 25
nboss_fibers_per_bundle = 20
nmangafibers = nbundles * opts.fibers_per_bundle

#- BOSS X and loglam spacing from PSF file
xb = pyfits.getdata(opts.input, 0)
llb = pyfits.getdata(opts.input, 2)
nbossfibers, nflux = xb.shape

#- MaNGA X, Y, and loglam to fill in
xm = N.zeros( (nmangafibers, nflux) )
ym = N.tile( N.arange(nflux, dtype='float'), nmangafibers).reshape( (nmangafibers, nflux) )
llm = N.zeros( (nmangafibers, nflux) )

#- Loop over bundles, interpolating x positions
nfb = nboss_fibers_per_bundle
nfm = opts.fibers_per_bundle
ixb = N.arange(nfb)
ixm = N.linspace(0, nfb-1, opts.fibers_per_bundle)

for b in range(nbundles):
    ib0 = b*nfb
    im0 = b*nfm
    for j in range(nflux):
        xm[im0:im0+nfm, j] = N.interp(ixm, ixb, xb[ib0:ib0+nfb, j])
        llm[im0:im0+nfm, j] = N.interp(ixm, ixb, llb[ib0:ib0+nfb, j])
                
#- HDU 0 : X
hdr = pyfits.getheader(opts.input, 0)
hdr['NSPEC'] = nmangafibers
pyfits.writeto(opts.output, xm, header=hdr, clobber=True)

#- HDU 1 : Y
pyfits.append(opts.output, ym, header=pyfits.getheader(opts.input, 1))

#- HDU 2 : LogLam
pyfits.append(opts.output, llm, header=pyfits.getheader(opts.input, 2))

#- HDU 3 : Model exponents; copy through
data = pyfits.getdata(opts.input, 3)
header = pyfits.getheader(opts.input, 3)
pyfits.append(opts.output, data, header=header)

#- HDU 4 : X and Y scale factors
#- HARDCODE uniformity of bundles
d = pyfits.getdata(opts.input, 4)
header = pyfits.getheader(opts.input, 4)
data = d.copy()
data.resize( (nmangafibers, ) )
for i in range(nmangafibers):
    b = i / opts.fibers_per_bundle
    data[i] = d[b*nboss_fibers_per_bundle]
    data[i]['IGROUP'] = b

pyfits.append(opts.output, data, header=header)

#- HDU 5 : Model images; copy through
data = pyfits.getdata(opts.input, 5)
header = pyfits.getheader(opts.input, 5)
pyfits.append(opts.output, data, header=header)



#-------------------------------------------------------------------------
#- Scratch code for debugging

sys.exit(0)
from bbspec.spec2d.psf import load_psf
psf0 = load_psf('spBasisPSF-PIX-r1-00131554.fits')
psf1 = load_psf('MaNGA-PSF-r1.fits')
x, y = psf0.x(), psf0.y()
xx, yy = psf1.x(), psf0.y()

P.subplot(211)
for i in range(0,65): P.plot(x[i], y[i])

P.subplot(212)
for i in range(0,65): P.plot(xx[i], yy[i])
    