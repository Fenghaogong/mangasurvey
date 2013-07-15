#!/usr/bin/env python

"""
Convolve known input spectra with a PSF to generate an output image

Examples:
manga_spec2pix.py -i block_manga_v2.fits -p MaNGA-PSF-r1.fits -c calib.fits \
  -o MaNGA-image-r1.fits -n
"""

import sys
import optparse
import numpy as N
import fitsio
from bbspec.spec2d.psf import load_psf

parser = optparse.OptionParser(usage = "%prog [options]")
parser.add_option("-i", "--input",  type="string",  help="input spectra")
parser.add_option("-p", "--psf",    type="string",  help="input psf")
parser.add_option("-c", "--calib",  type="string",  help="input calibration vectors")
parser.add_option("-o", "--output", type="string",  help="output image")
parser.add_option('-n', "--noise",  action='store_true', help='Add noise')

opts, args = parser.parse_args()

if opts.input is None or opts.psf is None or \
   opts.output is None or opts.calib is None:
    parser.print_help()
    print >> sys.stderr, "\nERROR: input, psf, and output are all required parameters"
    sys.exit(1)

#- Load spectra and PSF
spectra = fitsio.read(opts.input)
psf = load_psf(opts.psf)

#- Strip wavelength off of spectra array
#- TODO : redo data format
loglam = N.log10(spectra[0])
spectra = spectra[1:]

#- Trim to just 1 CCD of spectra
spectra = spectra[0:psf.nspec]
nfiber = spectra.shape[0]

#- HACK ALERT : derive channel from psf filename...
if opts.psf.count('r1')>0 or opts.psf.count('r2')>0:
    ihdu = 1
else:
    ihdu = 0
    
#- Read calibration vectors
calib = fitsio.read(opts.calib, ihdu)
hdr = fitsio.read_header(opts.calib, ihdu)
loglam_calib = hdr['CRVAL1'] + hdr['CDELT1'] * N.arange(calib.shape[1])

#- calib==0 is meaningless.  Set to large number so flux/calib ~= 0
ii = N.where(calib == 0.0)
calib[ii] = 1e12

#- Trim spectra to just region where we have calibration
ii = N.where( (loglam_calib[0] < loglam) & (loglam < loglam_calib[-1]) )[0]
loglam = loglam[ii]
spectra = spectra[:, ii]

#- NOTE:
#- calib comes from 1000 BOSS spectra.  Use modulo (%) to wrap around
#- fiber counting when applying to MaNGA spectra.  For now, just start
#- counting from 0 for both sp1 and sp2.

#- Convert spectra from ergs into electrons
#- electrons = flux / calib
for i in range(nfiber):
    xcalib = N.interp(loglam, loglam_calib, calib[i%calib.shape[0]])
    ### xcalib[xcalib==0.0] = N.inf
    spectra[i] /= xcalib

#- TEST
### spectra = spectra[0:10]
### spectra = spectra[:, 1800:1900]
### loglam = loglam[1800:1900]

#- Project spectra -> image with PSF
psf.resample(loglam)
image = psf.spec2pix(spectra, xyrange=(0, psf.npix_x, 0, psf.npix_y), verbose=True)

#- Add noise
pure_image = image.copy()
if opts.noise:
    image = N.random.poisson(image.clip(0, N.inf)).astype(float)
    
    read_noise = 2.0
    image += N.random.normal(scale=read_noise, size=image.shape)

fitsio.write(opts.output, image, clobber=True)
fitsio.write(opts.output, pure_image)
    
        

