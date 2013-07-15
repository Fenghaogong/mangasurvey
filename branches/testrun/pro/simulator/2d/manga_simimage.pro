; program manga_simimage
;
; Manga image simulator
;
; Program simulates image reconstruction at one wavelength slice given
; realistic observing scenarios.  Does not account for the spectral
; dimension, noise, or any detector characteristics, useful for
; exploring image reconstruction and observing strategy ONLY.
;
; Required input parameters:
; InpFile: String describing the location of the input FITS image
; pixscaleI: pixel scale of input image in arcsec/pixel
;
; Optional input parameters:
; xobj: X position of object of interest in input FITS file
; yobj: Y position of object of interest in input FITS file
; OutpFile: String name of output FITS image (default is 'SimImage.fits')
; sobfile: String name of simulated observing description file
;   Default is 'sample.sob'
; SimWave: Wavelength of simulation in Angstroms.
;   Default is 5500. (guide wavelength)
; /spiral: Set this flag to use spiral fiber layout. Default is serpentine.
;   Note, ONLY works with as-built bundles, not simulated bundles.
; /debug: Set this flag to run in debug mode and produce intermediary
;   files

; Output files:
; OutpFile: Output FITS image
; Assorted other intermediary files in the /working directory if /debug specified
;
; Examples:
; 1) Run a simulated observation of the sample image
;    'sampleinput.fits', assuming that it has pixel scale 0.1 arcsec/pixel.
;    Use a standard sobfile, but set the simulation wavelength to 9000 Angstroms
;    and run in debug mode to produce intermediary files.
;    The function call would look like
;    
;    manga_simimage,'sampleinput.fits',0.1,SimWave=9000.,/debug
;
; Author:
;   David R. Law (Dunlap Institute; drlaw@di.utoronto.ca)
;   First written: 04/29/2011
;   Last modified: 01/17/2013
;
; REVISION HISTORY:
;   v1: 17-Oct-2012  D. Law
;       Adapted previous code to this for-release.  Removed many
;       archaic options (9-pt dithers) and streamlined for use with
;       October 2012 vision of survey strategy.  All dialable
;       parameters have been moved out of command line into a
;       parameter file.
;   v1.2: 22-Dec-2012  D. Law
;       Rewrote formatting to interface with real MaNGA reduction pipeline,
;       call real bundle maps, etc.  Now includes dependencies on SDSS3
;       idlutils directory.
;   v1.3: 17-Jan-2013  D. Law
;       Rewrote to better handle real MaNGA metrology, change filepaths

; Function prototypes
forward_function mlbundlemap
forward_function mlcalifainterp

pro manga_simimage, InpFile, pixscaleI, xobj=xobj, yobj=yobj, OutpFile=OutpFile, sobfile=sobfile, SimWave=SimWave, spiral=spiral, debug=debug
  ; 32-bit integers
  Compile_Opt DEFINT32

  if (keyword_set(OutpFile)) then OutpFile=OutpFile $
  else OutpFile='SimImage.fits'

  if (keyword_set(sobfile)) then sobfile=sobfile $
  else sobfile='sample.sob'

  if (keyword_set(SimWave)) then SimWave=SimWave $
  else SimWave=5500.D

  ; Read simulated observations data description file
  ;
  ; This sets the single-valued parameters:
  ; pixscaleW: working pixel scale (arcsec/pixel) for the simulations
  ; pixscaleO: output pixel scale (arcsec/pixel) for the simulations
  ; boxxplier: Multiplier defining box size for the simulations as a fraction
  ;   of the bundle size
  ; ditherfrac: Fraction of the ideal dither offset that dithers move
  ;   (relevant for APOGEE co-observing tests)
  ; rseed: Seed for the random number generator
  ; guidelam: Guide wavelength in Angstroms
  ; tc: Outside temperature in degrees C
  ; rh: Relative humidity in percent
  ; atmp: Atmospheric pressure in mbar
  ; BType: Bundle map type to simulate (e.g., 7, 19, 37, 61, 91, 127, 169)
  ; BMfile: Actual bundle to use (e.g., cart01_127_1.bm)
  ; decl: Target declination in decimal degrees
  ; xyunc: UNKNOWN 1sigma uncertainty in fiber positions (in arcsec)
  ; xyerr: KNOWN 1sigma uncertainty in fiber positions (in arcsec)
  ; nexp: Total number of exposures (ndither*nnight-nbad)
  ;
  ; This also defines the dimension of other parameters, and populates
  ; the vectors appropriately:
  ; deadf: 3-element int vector containing numbers of fibers to call
  ;   dead (only relevant for simulated bundles; real bundles have dead
  ;   fibers handled during the call to mlbundlemap)
  ; PA: nexp-element vector (double type) defining position angle of bundle
  ; seeing: nexp-element vector (double type) defining seeing at guide
  ;   wavelength 5500 Ang
  ; exptime: nexp-element vector (double type) defining exposure time
  ;   in minutes
  ; ha: nexp-element vector (double type) defining hour angle in decimal hours
  ; dposn: nexp-element vector (integer type) defining dither position
  ;   (either position 1, 2, or 3)
  check=mlreadsob(sobfile,pixscaleW,pixscaleO,boxxplier,ditherfrac,rseed,crosstalk,guidelam,tc,rh,atmp,BType,BMfile,decl,knownxyerr,unknownxyerr,deadf,nexp,PA,seeing,exptime,ha,dposn)
  ; If mlreadsob failed, exit!
  if (check ne 0) then begin
     print,'mlreadsob failed to read .sob observing file properly.'
     return
  endif

  ; If debug keyword set, print out some key values here
  if (keyword_set(debug)) then begin
    print,'Working pixel scale: ',pixscaleW
    print,'Output pixel scale: ',pixscaleO
    print,'Bundle type: ',BType
    print,'Target declination: ',decl
    print,'Number of exposures: ',nexp
    print,'Dither positions: ',dposn
  endif

  ; Define the bundle map
  ; This establishes the baseline fiber locations in units of arcsec
  if keyword_set(spiral) then check=mlsimbmap(BType,nfiber,fiber_xcen,fiber_ycen,fiber_status,corediam,dither_rad,bundle_rmax,BMfile=BMfile,/spiral) $
  else check=mlsimbmap(BType,nfiber,fiber_xcen,fiber_ycen,fiber_status,corediam,dither_rad,bundle_rmax,BMfile=BMfile)
 
  ; If mlbundlemap failed, exit!
  if (check ne 0) then begin
     print,'mlbundlemap failed to define bundle properly.'
     return
  endif

  ; Now that all of the important parameters of the simulation have been
  ; defined, rescale input image to the proper working size
  InpImage=readfits(InpFile)
  xin=(size(InpImage))[1]; input x size in input pixels
  yin=(size(InpImage))[2]; input y size in input pixels
  typ=(size(InpImage))[(size(InpImage))[0]+1]
  xwin=xin*pixscaleI/pixscaleW; input x size in working pixels
  ywin=yin*pixscaleI/pixscaleW; input y size in working pixels
  fluxscale=(pixscaleI/pixscaleW)*(pixscaleI/pixscaleW)
  ; Resample the input image to the working pixel scale
  ; Use linear interpolation from pixel center, and conserve total flux
  ImageW=congrid(InpImage,xwin,ywin,/center,/interp)/fluxscale
  ; Now the target size in working pixels is
  Xsize=fix(bundle_rmax*2.*boxxplier/pixscaleW)
  Ysize=fix(bundle_rmax*2.*boxxplier/pixscaleW)

  ; Figure out what part of the image to use
  if (keyword_set(xobj)) then xobj=xobj*pixscaleI/pixscaleW $
  else xobj=xwin/2.
  if (keyword_set(yobj)) then yobj=yobj*pixscaleI/pixscaleW $
  else yobj=xwin/2.

  ; If necessary simulation box size is smaller than the input image, crop
  ; out the proper area from the center
  if (Xsize lt xwin) then begin
    print,'Simulation box is smaller than input image: cropping X dimension.'
    minx=xobj-Xsize/2. > 0; starting point is either desired position minus half box size
                          ; or zero (whichever is greater)
    maxx=minx+Xsize-1; End point is whatever the starting point is, plus box size
    ImageW=ImageW[minx:maxx,*]
  endif
  if (Ysize lt ywin) then begin
    print,'Simulation box is smaller than input image: cropping Y dimension.'
    miny=yobj-Ysize/2. > 0; starting point is either desired position minus half box size
                          ; or zero (whichever is greater)
    maxy=miny+Ysize-1; End point is whatever the starting point is, plus box size
    ImageW=ImageW[*,miny:maxy]
  endif
  ; If necessary simulation box size is larger than the input image, pad
  ; the input image with a zero buffer (and tell the user)
  if ((Xsize gt xwin)or(Ysize gt ywin)) then begin
    print,'Simulation box is bigger than input image: padding to fill.'
    deltax=Xsize-xwin > 0
    deltay=Ysize-ywin > 0

    temp=make_array(Xsize, Ysize, type=typ, value=0.)
    temp[deltax/2,deltay/2]=ImageW
    ImageW=temp
  endif
  if (keyword_set(debug)) then writefits,'working/imageworking.fits',ImageW

  ; Rescale fiber locations, radii, and fiducial dither distance to pixel units
  fiber_xcen=fiber_xcen/pixscaleW
  fiber_ycen=fiber_ycen/pixscaleW
  dither_rad=dither_rad/pixscaleW
  corediam=corediam/pixscaleW
  knownxyerr=knownxyerr/pixscaleW
  unknownxyerr=unknownxyerr/pixscaleW

  ; Assume fiber bundle centered on the center of the field
  ; (centering already done above)
  bundle_xcen=Xsize/2.
  bundle_ycen=Ysize/2.

  ; If certain fibers were specified to be 'dead' in the .sob file, make
  ; them so in fiber_status
  ndead=0
  for i=1,nfiber do begin
    if ((i eq deadf[0])or(i eq deadf[1])or(i eq deadf[2])) then fiber_status[i-1]=1
    if (fiber_status[i-1] ne 0) then ndead=ndead+1
  endfor
  if (keyword_set(debug)) then print,'There are ',ndead,' dead fibers.'

  ; Set up random number generator for xy offsets (known and unknown)
  ; Random draws for both x and y offset, where throw of each is sqrt(2)
  ; times the total throw.
  if (rseed eq -1) then rgenK=RANDOMN(systime_seed,nexp,2)*knownxyerr/sqrt(2.) $
  else rgenK=RANDOMN(rseed,nexp,2)*knownxyerr/sqrt(2.)
  if (rseed eq -1) then rgenU=RANDOMN(systime_seed+1,nexp,2)*unknownxyerr/sqrt(2.) $
  else rgenU=RANDOMN(rseed+1,nexp,2)*unknownxyerr/sqrt(2.)

  ; Define a matrix of fiber positions for all nexp observations
  xmatrix=fltarr(nfiber,nexp)
  ymatrix=fltarr(nfiber,nexp)
  for j=0,nexp-1 do begin
    ; Rotate by position angle of bundle
    xtemp=fiber_xcen*cos(pa[j]*3.1415926536/180.)+fiber_ycen*sin(pa[j]*3.1415926536/180.)
    ytemp=-fiber_xcen*sin(pa[j]*3.1415926536/180.)+fiber_ycen*cos(pa[j]*3.1415926536/180.)
    ; Add to base bundle center
    xtemp=xtemp+bundle_xcen
    ytemp=ytemp+bundle_ycen
    ; Dither offset to position 1, 2, or 3 (see diagram in obs strategy)
    if (dposn[j] eq 1) then dangle=240.*3.1415926536/180.
    if (dposn[j] eq 2) then dangle=0.*3.1415926536/180.
    if (dposn[j] eq 3) then dangle=120.*3.1415926536/180.
    xtemp=xtemp+ditherfrac*dither_rad*cos(dangle)
    ytemp=ytemp+ditherfrac*dither_rad*sin(dangle)
    ; DAR offset
    dr=mldar(ha[j],decl,SimWave,parang,xdar,ydar,waveREF=guidelam,TC=tc,RH=rh,P=atmp)
    xtemp=xtemp+xdar/pixscaleW
    ytemp=ytemp+ydar/pixscaleW
    ; Known XY bundle positioning error.  Random draws for each exposure
    ; from a Gaussian probability distribution with specified sigma.
    xtemp=xtemp+rgenK[j,0]
    ytemp=ytemp+rgenK[j,1]
    ; Set fiber positions
    xmatrix[*,j]=xtemp
    ymatrix[*,j]=ytemp
 endfor
 
  ; Define a matrix of 2d arrays describing the fiber maps
  ; fiber maps are -1 everywhere except where a fiber covers, where
  ; it is equal to the fiber number.
  fibermap=replicate(-1,Xsize,Ysize,nexp)
  coordX = rebin(dindgen(Xsize), [Xsize, Ysize])
  coordY = rebin(transpose(dindgen(Ysize)), [Xsize, Ysize])
  temp=intarr(Xsize,Ysize)
  for j=0,nexp-1 do begin
    temp[*,*]=-1
    for i=0,nfiber-1 do begin
      radmap=sqrt((coordX-xmatrix[i,j])^2+(coordY-ymatrix[i,j])^2)
      if (fiber_status[i] eq 0) then temp[where(radmap lt corediam/2.)]=i
    endfor
    fibermap[*,*,j]=temp
  endfor
  ; Determine fiber area in units of working pixels
  ; Use whatever the first non-dead fiber is to determine this
  fibareaW=(size(where(fibermap[*,*,0] eq (where(fiber_status eq 0))[0])))[1]
  if (keyword_set(debug)) then print,'fibareaW = ',fibareaW
  if (keyword_set(debug)) then writefits,'working/fibermap.fits',fibermap

  ; Generate a coverage map quantifying how many times each part of the
  ; field was observed
  CoverageMap=intarr(Xsize,Ysize)
  for j=0,nexp-1 do CoverageMap += fibermap[*,*,j] ne -1
  if (keyword_set(debug)) then writefits,'working/coveragemap.fits',CoverageMap

  ; Determine the amount of flux down each fiber, during each exposure
  fiberinput=dblarr(nfiber,nexp)
  for j=0,nexp-1 do begin
    ; Determine what the input image looks like when convolved with seeing
    blur1=mlfilterimg(ImageW,fwhm_gaussian=seeing[j]*((SimWave/5500.)^(-0.2))/pixscaleW,/all_pixels)*9./13.
    blur2=mlfilterimg(ImageW,fwhm_gaussian=2.*seeing[j]*((SimWave/5500.)^(-0.2))/pixscaleW,/all_pixels)*4./13.
    blur=blur1+blur2
    ; Figure out what is subtended by each fiber
    for i=0,nfiber-1 do begin
      fiberinput[i,j]=total(blur*(fibermap[*,*,j] eq i))
    endfor
  endfor

  ; In debug mode, make an 'image' for each exposure
  if (keyword_set(debug)) then begin
    IndivImages=fltarr(Xsize,Ysize,nexp)
    for j=0,nexp-1 do begin
      for i=0,nfiber-1 do begin
        IndivImages[*,*,j] += fiberinput[i,j]*(fibermap[*,*,j] eq i)/fibareaW; Flux conservation
      endfor
    endfor
    writefits,'working/IndivImages.fits',IndivImages
  endif

  ; Re-order info into 1-d vectors for interpolation routine
  ; Make sure to rescale x,y locations to output grid scale
  xout=fltarr(nexp*nfiber)
  yout=fltarr(nexp*nfiber)
  statout=fltarr(nexp*nfiber)
  fluxout=fltarr(nexp*nfiber)
  XsizeOUT=fix(Xsize*pixscaleW/pixscaleO)
  YsizeOUT=fix(Ysize*pixscaleW/pixscaleO)
  print,'XsizeOUT=',XsizeOUT
  for j=0,nexp-1 do begin
    for i=0,nfiber-1 do begin
      ; Rescale x and y locations
      ; If unknown bundle positioning errors were specified, add them here
      xout[i+j*nfiber]=(xmatrix[i,j]+rgenU[j,0])*pixscaleW/pixscaleO
      yout[i+j*nfiber]=(ymatrix[i,j]+rgenU[j,1])*pixscaleW/pixscaleO
      ; Adjust flux to different pixel scale
      fluxout[i+j*nfiber]=fiberinput[i,j]*(pixscaleW/pixscaleO)*(pixscaleW/pixscaleO)
      statout[i+j*nfiber]=fiber_status[i]
    endfor
  endfor

  ; Interpolate to a regular grid
  ; Emprically, using r=1.6'' and sigma=0.7 arcsec gives the best image
  ; quality tested so far, regardless of seeing (with califa algorithm)
  interpimage=mlcalifainterp(xout,yout,fluxout,statout,[XsizeOUT,YsizeOUT],1.6/pixscaleO,0.7/pixscaleO)
  writefits,OutpFile,interpimage
return
end

  
