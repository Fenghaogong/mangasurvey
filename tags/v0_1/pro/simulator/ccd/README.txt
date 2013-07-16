This is a collection of routines that simulates MaNGA spectra on the BOSS CCDs.

It is split into two parts:
Part I runs in IDL and produces 1500 row stacked spectra and a wavelength solution
  in on-sky flux units.
Part II runs in Python and uses the output of part I to make the CCD images.

NOTE that the files
block_manga.fits
spBasisPSF-PIX-b1-00131554.fits
spBasisPSF-PIX-r1-00131554.fits
MaNGA-PSF-b1.fits
MaNGA-PSF-r1.fits
calib.fits

have been moved to the mangadb/simfiles/ccd/ directory and gzipped to cut down on file sizes.  
If using the CCD simulator they will
have to be put in the correct directory again on your local copy.

#####################################################

Part I: Use the routine manga_makeccdspec_part1.pro to produce the row stacked spectra
  in block_manga.fits
  All parameters are hard coded.

#####################################################

Part 2: Python routines from Steven Bailey

Bits and pieces for simulating MaNGA spectra on a CCD

Stephen Bailey, LBL
StephenBailey@lbl.gov
Summer 2012

### Quick Start Guide ###

1. Get the python requirements listed below
2. Run:
    
    manga_spec2pix.py -i block_manga.fits -c calib.fits -n \
        -p MaNGA-PSF-r1.fits -o MaNGA-specimage-r1.fits
        
    manga_spec2pix.py -i block_manga.fits -c calib.fits -n \
        -p MaNGA-PSF-b1.fits -o MaNGA-specimage-b1.fits
        
3. Update block_manga.fits with other spectra and rerun.
4. Use/Modify squash_psf.py to generate other MaNGA PSFs using the
   BOSS PSFs spBasisPSF-PIX-*-00131554.fits as inputs.
5. Iterate

### Contents ###

README.txt - This file
Data:
    spBasisPSF-PIX-*-00131554.fits  - BOSS PSFs
    MaNGA-PSF-*.fits                - MaNGA PSFs modified from the BOSS PSFs
    calib.fits                      - BOSS fiber flux = photons*calib
    block_manga.fits                - MaNGA input spectra from David Law
    MaNGA-image-r1.fits             - Output simulated MaNGA spectra image
Code:
    manga_spec2pix.py               - spectra+psf+calib -> image
    get_plate_calib.py              - Extract calib vectors from BOSS data
    squash_psf.py                   - Convert BOSS PSFs to MaNGA PSFs

### Python Code Requirements ###

bbspec product from BOSS svn repository
    Add $BBSPEC_DIR/python to $PYTHONPATH
External python packages:
    numpy
    scipy
    matplotlib
    pyfits [http://www.stsci.edu/institute/software_hardware/pyfits/Download]
    fitsio [https://github.com/esheldon/fitsio]
    
If you don't already have python + scientific packages installed, the
easist route is to get the Enthought Python Distribution [EPD]
from http://enthought.com/ [Free for academics] followed by
"easy_install fitsio".

### From Spectra to CCD Images ###

Simulated MaNGA CCD images where generated using manga_spec2pix.py

Usage: manga_spec2pix.py [options]

Options:
  -h, --help            show this help message and exit
  -i INPUT, --input=INPUT
                        input spectra
  -p PSF, --psf=PSF     input psf
  -c CALIB, --calib=CALIB
                        input calibration vectors
  -o OUTPUT, --output=OUTPUT
                        output image
  -n, --noise           Add noise

e.g.

manga_spec2pix.py -i block_manga_v2.fits -p MaNGA-PSF-r1.fits -c calib.fits \
  -o MaNGA-image-r1.fits -n

The input spectra from David Law are in block_manga.fits.  Row 0 is
the wavelength [angstroms]; additional rows are spectra in units of
ergs/s/cm2/A * 10^17.  The manga_spec2pix.py code uses the PSF file to
know how many spectra to project to where on a CCD.  It will only use
the first N spectra from the spectra file.  i.e. this example file
has 1500 spectra, but you will need to generate two separate spectra files
to separately simulate spectrographs 1 and 2.

### Generating MaNGA PSFs ###

This directory comes with MaNGA-PSF-*.fits PSF files with 750 fibers
split into 25 bundles of 30 fibers each.  These were generated with:

squash_psf.py -i spBasisPSF-PIX-r1-00131554.fits -o MaNGA-PSF-r1.fits
squash_psf.py -i spBasisPSF-PIX-b1-00131554.fits -o MaNGA-PSF-b1.fits

### Calibration Vectors ###

You probably won't need to to regenerate the calib.fits file, but for
the record...

calib.fits has the per-fiber calibration vectors for the B1 spectrograph
(HDU 0) and the R1 spectrograph (HDU 1).  NAXIS2 = 1000 = number of BOSS
fibers.  NAXIS1 = number of wavelength bins.
log10(wavelength) = CRVAL1 + CDELT1 * range(0..NAXIS2)

These were extracted from a directory of BOSS spectra using:

get_plate_calib.py -i $BOSS_SPECTRA/ -o calib.fits 

