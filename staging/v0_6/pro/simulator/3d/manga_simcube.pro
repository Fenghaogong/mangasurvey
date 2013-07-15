; program manga_simcube
;
; Manga data cube simulator
;
; Program simulates observation of a real MaNGA data cube.
; Takes an input flux-calibrated cube of 'real source' info
; and translates it to an 'as-recovered' cube.
;
; Required input parameters:
; InpCubeFile: String describing the location of the input FITS cube.
;   This cube must be flux calibrated in units of 1e-17 erg/s/cm^2/Angstrom
; pixscaleI: pixel scale of input image in arcsec/pixel
; InpWaveFile: String describing the location of a 1d FITS file
;   containing the wavelength information for the data cube.
;   Length must match spectral length of cube.
; wavestart: Starting wavelength for simulation in Angstroms
; wavestop: Ending wavelength for simulation in Angstroms
;
; Optional input parameters:
; xobj: X position of object of interest in input FITS cube
; yobj: Y position of object of interest in input FITS cube
; CCDfile: String name of output FITS spectra on CCD (default is 'CCDspec.fits')
; RSSfile: String name of output RSS FITS spectra (default is 'RSS.fits'),
;  this is row-stacked spectra, i.e. 'reduced' spectra that haven't been
;  assembled into a cube.
;  The error for each pixel in the RSSfile is written to err_RSSfile
; OutpFile: String name of output FITS cube (default is 'SimCube.fits')
; OutpWaveFile: String name of output FITS wavelength reference file
;   (default is 'SimWave.fits')
; sobfile: String name of simulated observing description file
;   Default is 'sample.sob'
; /ctflux: Set this flag to assume that CCD extraction cannot separate
;   flux introduced via crosstalk.  Default is to assume that it can.
; /nospecsmooth: Set this flag to NOT smooth input spectra to BOSS
;   resolution (assume already correct spec resolution)
; /spiral: Set this flag to use spiral fiber layout. Default is serpentine.
;   Note, ONLY works with as-built bundles, not simulated bundles.
; /debug: Set this flag to run in debug mode and produce intermediary
;   files.  This will slow down the program significantly, and write
;   quite a lot of large files to disk.
; /verbose: Alternative to /debug that prints progress messages to the terminal,
;   but doesn't take up extra time or write large files to disk
;
; Output files:
; CCDfile: Spectra of individual fibers as they would fall on CCD
;   (albeit on a different wavelength grid)
; OutpFile: Output FITS cube
; OutpWaveFile: Output FITS file listing wavelengths
; Assorted other intermediary files in the /working directory if /debug specified
;
; Examples:
; 1) Perform a simulated observation of the sample image
;     'samplecube.fits', assuming that it has pixel scale 1.0 arcsec/pixel and a
;     wavelength solution given by 'samplewave.fits'.
;     Note that these sample input files are the CALIFA early observations of
;     UGC00233 obtained from http://www.caha.es/CALIFA/public_html/?q=content/publications
;     Use a standard observing setup file, but set the simulation wavelength range to
;     bracket 6500-6800 Angstroms and run in verbose mode to print status messages:
;
;     manga_simcube,'samplecube.fits',1.0,'samplewave.fits',6500.,6800.,/verbose
;
; Author:
;   David R. Law (Dunlap Institute; drlaw@di.utoronto.ca)
;   First written: 04/29/2011
;   Last modified: 05/30/2013
;
; REVISION HISTORY:
;   v1: 24-Oct-2012  D. Law
;       Adapted previous code to this for-release.  Removed many
;       archaic options (9-pt dithers) and streamlined for use with
;       October 2012 vision of survey strategy.  All dialable
;       observing parameters have been moved out of command line into a
;       parameter file.
;   v1.1: 29-Oct-2012  D. Law
;       Tweaked random number generator, added info to output FITS headers
;   v1.2: 27-Nov-2012  D. Law
;       Added intermediate output of flux down fibers to facilitate
;       importing this into the CCD simulator.
;   v1.3: 22-Dec 2012  D. Law
;       Rewrote formatting to interface with real MaNGA reduction pipeline,
;       call real bundle maps, etc.  Now includes dependencies on SDSS3
;       idlutils directory.
;   v1.4: 22-Dec 2012  N. Drory & A. Weijmans
;       Added optional cross-talk between adjacent fibers on the slit.
;   v1.5: 23-Dec 2012  N. Drory & A. Weijmans
;       Added optional output of errors of RSS file
;   v1.6: 28-Dec 2012.  D. Law
;       Zeroed out unphysical negative flux values from object.
;   v1.7: 14-Jan-2013   D. Law
;       Tweaked crosstalk handling to read from sob file instead of
;       command line.  Default is now that crosstalk is removed from
;       flux and enters only in noise.  Change so that crosstalk is
;       NOT removed from flux using the /ctflux flag.
;       Changed filepath for database files.
;   v1.8: 17-Jan-2013  D. Law
;       Rewrote to better handle real MaNGA metrology, change filepaths.
;   v1.9: 22-Jan-2013  D. Law
;       Added treatment of noise in 3d cubes.
;   v1.10: 25-Jan-2013  D. Law
;       Added flag 'nospecsmooth' to NOT apply spectral smoothing to BOSS resolution
;   v1.11: 04-Mar-2013  D. Law
;       Tweaked flux calibration of output cubes
;   v1.12: 30-May-2013  D. Law
;       Added /flatsky option which uses a very simplified sky model
;       that has no lines and is just a tilted line in log(flux) space.
;
; Function prototypes
forward_function mlbundlemap
forward_function mlcalifainterp

pro manga_simcube, InpCubeFile, pixscaleI, InpWaveFile, wavestart, wavestop, xobj=xobj, yobj=yobj, CCDfile=CCDfile, RSSfile=RSSfile, OutpFile=OutpFile, OutpWaveFile=OutpWaveFile, sobfile=sobfile, spiral=spiral, debug=debug, verbose=verbose, ctflux=ctflux, nospecsmooth=nospecsmooth, flatsky=flatsky
  ; 32-bit integers
  Compile_Opt DEFINT32

  if (keyword_set(OutpFile)) then OutpFile=OutpFile $
  else OutpFile='SimCube.fits'

  if (keyword_set(CCDfile)) then CCDfile=CCDfile $
  else CCDfile='CCDspec.fits'

  if (keyword_set(RSSfile)) then RSSfile=RSSfile $
  else RSSfile='RSS.fits'
  RSSfile_err = 'err_'+RSSfile

  if (keyword_set(OutpWaveFile)) then OutpWaveFile=OutpWaveFile $
  else OutpWaveFile='SimWave.fits'

  if (keyword_set(sobfile)) then sobfile=sobfile $
  else sobfile='sample.sob'

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;; Set up various inputs and system parameters
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

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
  ; deadf: 3-element int vector containing address of fibers to call
  ;   dead (only relevant for simulated bundles; real bundles have dead
  ;   fibers handled during the call to mlbundlemap)
  ; PA: nexp-element vector (double type) defining position angle of bundle
  ; seeing: nexp-element vector (double type) defining seeing at guide
  ;   wavelength
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
  ; Convert exposure time from minutes to seconds
  exptime=exptime*60.

  ; If debug keyword set, print out some key values here
  if ((keyword_set(debug))or(keyword_set(verbose))) then begin
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

  ; Define spectrograph calibration information
  ; Calibration information is stored in the file boss_calib.dat, which
  ; was assembled from boss_sky.dat, boss_throughput.dat, and boss_nphotons.dat
  calibfile=getenv('MANGAROOT')+'/mangadb/' +getenv('MANGAVER') + '/bosscal/56280/calibmatrix.fits'
  calibmatrix=readfits(calibfile)
  nwave=(size(calibmatrix))[2]
  ; wave is the master wavelength vector for the simulation, log spacing
  ; from 3500 to 10000 Angstroms
  wave=dblarr(nwave)
  wave[*]=calibmatrix[0,*]
  ; Ensure appropriate wavelength boundaries set
  if (wavestart lt 3556.) then wavestart=3556.
  if (wavestop gt 10356.) then wavestop=10356.

  ; Make the flatsky spectrum for (possible) reference
  flatskyspectrum=10.^(findgen(nwave)*(0.3/nwave)-0.8)
  
  ; Isolate the relevant parts of the master wavelength vector
  temp=min(wave-wavestart,kstart,/absolute,/nan); Determinine corresponding kstart index
  temp=min(wave-wavestop,kstop,/absolute,/nan); Determinine corresponding kstop index
  ; Recrop wavelength and calibration arrays accordingly
  calibmatrix=calibmatrix[*,kstart:kstop]
  flatskyspectrum=flatskyspectrum[kstart:kstop]
  wave=wave[kstart:kstop]
  nwave=(size(calibmatrix))[2]
  ; Write wavelength file to disk
  writefits,OutpWaveFile,wave



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;; Put the input data cube in working form
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ; Read in data cube and wavelength reference file
  InpCube=readfits(InpCubeFile)
  xin=(size(InpCube))[1]; input x size in input pixels
  yin=(size(InpCube))[2]; input y size in input pixels
  zin=(size(InpCube))[3]; input z size in input pixels
  typ=(size(InpCube))[(size(InpCube))[0]+1]
  ; Input wavelength reference file
  InpWave=readfits(InpWaveFile)
  if (zin ne (size(InpWave))[1]) then begin
    print,'ERROR: Spectral length mismatch between input cube and'
    print,'specified wavelength solution'
    return
  endif

  ; Start by cropping out *roughly* the right spectral range
  temp=min(InpWave-wavestart,zstart,/absolute,/nan); Determinine corresponding zstart index
  zstart= zstart -2 > 0; Extra couple slices for later interpolation, but make sure valid index
  temp=min(InpWave-wavestop,zstop,/absolute,/nan); Determinine corresponding zstop index
  zstop= zstop +2 < zin-1; Extra couple slices for later interpolation, but make sure valid index
  InpCube=InpCube[*,*,zstart:zstop]
  InpWave=InpWave[zstart:zstop]
  zin=(size(InpCube))[3]

  ; Rescale image dimensions to correct pixel scale
  xwin=xin*pixscaleI/pixscaleW; input x size in working pixels
  ywin=yin*pixscaleI/pixscaleW; input y size in working pixels
  fluxscale=(pixscaleI/pixscaleW)*(pixscaleI/pixscaleW)
  ; Resample the input cube to the working pixel scale
  ; Use linear interpolation from pixel center, and conserve total flux
  tempcube=congrid(InpCube,xwin,ywin,zin,/center,/interp)/fluxscale
  ; Now the target size in working pixels is
  Xsize=fix(bundle_rmax*2.*boxxplier/pixscaleW)
  Ysize=fix(bundle_rmax*2.*boxxplier/pixscaleW)

  ; Figure out what part of the image to use
  if (keyword_set(xobj)) then xobj=xobj*pixscaleI/pixscaleW $
  else xobj=xwin/2.
  if (keyword_set(yobj)) then yobj=yobj*pixscaleI/pixscaleW $
  else yobj=ywin/2.

  ; If necessary simulation cube is smaller than the input cube, crop
  ; out the proper region
  if (Xsize lt xwin) then begin
    print,'Simulation box is smaller than input image: cropping X dimension.'
    minx=xobj-Xsize/2. > 0; starting point is either desired position minus half box size
                          ; or zero (whichever is greater)
    maxx=minx+Xsize-1; End point is whatever the starting point is, plus box size
    tempcube=tempcube[minx:maxx,*,*]
  endif
  if (Ysize lt ywin) then begin
    print,'Simulation box is smaller than input image: cropping Y dimension.'
    miny=yobj-Ysize/2. > 0; starting point is either desired position minus half box size
                          ; or zero (whichever is greater)
    maxy=miny+Ysize-1; End point is whatever the starting point is, plus box size
    tempcube=tempcube[*,miny:maxy,*]
  endif

  ; If necessary simulation box size is larger than the input image, pad
  ; the input image with a zero buffer (and tell the user)
  if ((Xsize gt xwin)or(Ysize gt ywin)) then begin
    print,'Simulation box is bigger than input image: padding to fill.'
    deltax=Xsize-xwin > 0
    deltay=Ysize-ywin > 0
    temp=make_array(Xsize, Ysize, zin,type=typ, value=0.D)
    for k=0,zin-1 do temp[deltax/2,deltay/2,k]=tempcube[*,*,k]
    tempcube=temp
  endif

  ; Now rescale/crop/pad appropriately in the spectral dimension
  ; to put the cube in a standard working format
  ; Loop over spatial pixels, smooth the spectra by the BOSS spectral
  ; resolution and interpolate the spectra
  ; the appropriate format.
  ; BOSS spectral resln is roughly R=1560 at 3700 Ang and R=2270
  ; 6000 Ang (blue channel), and R=1850 at 6000 Ang and R=2650 at 
  ; 9000 Ang (red channel).  This works out to FHWM in Angstroms
  ; varying from 2.37 Ang to 2.64 Ang in the blue channel and
  ; 3.24 to 3.40 Ang in the red.  These are pretty constant, so fix
  ; spectral resln to 2.51 Ang at lambda<6000 Ang, and to 3.32 Ang
  ; at lambda>6000 Ang.
  tempspec=dblarr(zin)
  CubeW=dblarr(Xsize,Ysize,nwave)
  for i=0,Xsize-1 do begin
  if (((keyword_set(debug))or(keyword_set(verbose)))and((i mod 10)eq 0)) then print,'Resampling input spectra: ',(i*100.)/Xsize,'% complete'
    for j=0,Ysize-1 do begin
      tempspec=tempcube[i,j,*]
      ; Unless /nospecsmooth set, smooth to BOSS resolution
      if not(keyword_set(nospecsmooth)) then begin
        if ((wavestart+wavestop)/2. le 6000.) then tempspec2=mlgaussfold(InpWave,tempspec,2.51) $
        else tempspec2=mlgaussfold(InpWave,tempspec,3.32)
        if ((i eq 0)and(j eq 0)) then print,'Smoothing to BOSS spectral resolution'
      ; If /nospecsmooth set, don't smooth spectra
      endif else begin
        tempspec2=tempspec
        if ((i eq 0)and(j eq 0)) then print,'No spectral smoothing selected'
      endelse

      tempspec3=interpol(tempspec2,InpWave,wave)
      CubeW[i,j,*]=tempspec3
    endfor
  endfor

  ; Zero out interpolated information beyond the bounds of original
  ; input wavelength range
  if ((keyword_set(debug))or(keyword_set(verbose))) then print,'Zeroing bad values'
  zerolow=where(wave lt InpWave[0])
  zerohigh=where(wave gt InpWave[zin-1])
  if (zerolow[0] ne -1) then CubeW[*,*,zerolow]=0.
  if (zerohigh[0] ne -1) then CubeW[*,*,zerohigh]=0.
  if (keyword_set(debug)) then writefits,'working/cubeworking.fits',CubeW,/compress
  if ((keyword_set(debug))or(keyword_set(verbose))) then print,'Got a working form cube.'

  ; Free memory by deleting unnecessary matrices
  mlundefine, InpCube
  mlundefine, tempcube

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;; Now we've got a functional cube in working form, use it!
;;;;; Make the fiber position matrices, and define other matrices
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ; Rescale fiber locations, radii, and fiducial dither distance to working pixel units
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
  if ((keyword_set(debug))or(keyword_set(verbose))) then print,'There are ',ndead,' dead fibers.'

  ; Set up random number generator for xy offsets (known and unknown)
   temp=randomn(systime_seed)
   if (rseed eq -1) then rseed=systime_seed
   seed=rseed; Reset RNG
   rgenK=RANDOMN(seed,nexp,2)*knownxyerr
   seed=rseed+1; Reset RNG
   rgenU=RANDOMN(seed,nexp,2)*unknownxyerr

  ; Define a matrix of fiber positions for all nexp observations
  ; This is also a function of wavelength
  xmatrix=fltarr(nfiber,nexp,nwave)
  ymatrix=fltarr(nfiber,nexp,nwave)

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
    ; Known XY bundle positioning error.  Random draws for each exposure
    ; from a Gaussian probability distribution with specified sigma.
    xtemp=xtemp+rgenK[j,0]
    ytemp=ytemp+rgenK[j,1]
    ; Add in the DAR and set fiber positions at each wavelength
    for k=0,nwave-1 do begin
      dr=mldar(ha[j],decl,wave[k],parang,xdar,ydar,waveREF=guidelam,TC=tc,RH=rh,P=atmp)
      xmatrix[*,j,k]=xtemp+xdar/pixscaleW
      ymatrix[*,j,k]=ytemp+ydar/pixscaleW
    endfor
  endfor

  ; Define a matrix of 2d arrays describing the fiber maps
  ; fiber maps are -1 everywhere except where a fiber covers, where
  ; it is equal to the fiber number.
  coordX = rebin(dindgen(Xsize), [Xsize, Ysize])
  coordY = rebin(transpose(dindgen(Ysize)), [Xsize, Ysize])
  radmap = replicate(0.,Xsize,Ysize)
  ; xmin,xmax,ymin,ymax say what region is important for each fiber
  xmin=fix(Xsize/2.-corediam)
  xmax=fix(Xsize/2.+corediam)
  ymin=fix(Ysize/2.-corediam)
  ymax=fix(Ysize/2.+corediam)
  workingregion=[n_elements(radmap[xmin:xmax]),n_elements(radmap[ymin:ymax])]
  ; Define lots of arrays that will be used within the upcoming loop over
  ; wavelength space.
  fibermap=replicate(-1,Xsize,Ysize,nexp)
  temp=intarr(Xsize,Ysize)
  ; Amount of flux down each fiber, during each exposure
  fiberinput=dblarr(nfiber,nexp)
  skyphot=dblarr(nfiber,nexp)
  recflux=dblarr(nfiber,nexp)
  recflux_err=dblarr(nfiber,nexp)
  ; Fiber area in units of working pixels
  fibareaW=0
  ; Images for each exposure (debug only)
  IndivImages=fltarr(Xsize,Ysize,nexp)
  ; Matrices for image interpolation
  xout=fltarr(nexp*nfiber)
  yout=fltarr(nexp*nfiber)
  statout=fltarr(nexp*nfiber)
  fluxout=fltarr(nexp*nfiber)
  errout=fltarr(nexp*nfiber)
  varout=fltarr(nexp*nfiber)
  ; Output sizes
  XsizeOUT=fix(Xsize*pixscaleW/pixscaleO)
  YsizeOUT=fix(Ysize*pixscaleW/pixscaleO)
  ; Final output cube
  FinalCube=fltarr(XsizeOUT,YsizeOUT,nwave)
  FinalErrCube=fltarr(XsizeOUT,YsizeOUT,nwave)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;; Loop through wavelength, doing the remaining caculations
;;;;; There are better ways to do this, but implement later
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; Make a big array that will be used to determine intermediate CCD
; spectra (for output file only, not used in calculations here)
CCDspec=dblarr(nwave,nexp*nfiber)
RSS=dblarr(nwave,nexp*nfiber)
RSS_err=dblarr(nwave,nexp*nfiber)

for k=0,nwave-1 do begin
  if (k eq 0) then T=systime()
  Tnow=systime()
  tperk=mltimer(T,Tnow)/(k > 1)
  if (((keyword_set(debug))or(keyword_set(verbose)))and(k eq 1)) then print,'Looping over wavelength space: ',(k*100.)/nwave,'% complete, ', tperk*(nwave-k),' seconds remaining'
  if (((keyword_set(debug))or(keyword_set(verbose)))and((k mod 10)eq 0)) then print,'Looping over wavelength space: ',(k*100.)/nwave,'% complete, ', tperk*(nwave-k),' seconds remaining'

  fibermap[*,*,*]=-1

  for j=0,nexp-1 do begin
    temp[*,*]=-1
    for i=0,nfiber-1 do begin
      xmin=xmatrix[i,j,k]-corediam > 0
      xmax=xmatrix[i,j,k]+corediam < Xsize-1
      ymin=ymatrix[i,j,k]-corediam > 0
      ymax=ymatrix[i,j,k]+corediam < Ysize-1
      radmap[*,*]=corediam; radmap is corediam by default (big enough to ignore)
      ; Calculate actual radius within a sensible box
      radmap[xmin:xmax,ymin:ymax]=sqrt((coordX[xmin:xmax,ymin:ymax]-Xmatrix[i,j,k])^2 +(coordY[xmin:xmax,ymin:ymax]-Ymatrix[i,j,k])^2)
      ; Assign fiber ID
      if (fiber_status[i] eq 0) then temp[where(radmap lt corediam/2.)]=i
    endfor
    fibermap[*,*,j]=temp
  endfor
  ; Determine fiber area in units of working pixels
  ; Use whatever the first non-dead fiber is to determine this
  ; This only has to be called once, not at all wavelengths
  if (k eq 0) then begin
    fibareaW=(size(where(fibermap[*,*,0] eq (where(fiber_status eq 0))[0])))[1]
    if (keyword_set(debug)) then print,'fibareaW = ',fibareaW
  endif

  ; Generate a coverage map quantifying how many times each part of the
  ; field was observed.  This is only done in debug mode
  if (keyword_set(debug)) then begin
    CoverageMap=intarr(Xsize,Ysize)
    for j=0,nexp-1 do CoverageMap += fibermap[*,*,j] ne -1
  endif

  ; If in debug mode, spit out the various maps produced so far
  if (keyword_set(debug)) then begin
    fmfile=strcompress('working/fibermap'+strtrim(k)+'.fits',/remove_all)
    writefits,fmfile,fibermap
    covfile=strcompress('working/coveragemap'+strtrim(k)+'.fits',/remove_all)
    writefits,covfile,CoverageMap
  endif

  ; Determine the amount of flux down each fiber, during each exposure
  for j=0,nexp-1 do begin
    ; Determine what the input image looks like when convolved with seeing
    blur1=mlfilterimg(CubeW[*,*,k],fwhm_gaussian=seeing[j]*((wave[k]/5500.)^(-0.2))/pixscaleW,/all_pixels)*9./13.
    blur2=mlfilterimg(CubeW[*,*,k],fwhm_gaussian=2.*seeing[j]*((wave[k]/5500.)^(-0.2))/pixscaleW,/all_pixels)*4./13.
    blur=blur1+blur2
    ; Figure out what is subtended by each fiber
    for i=0,nfiber-1 do begin
      fiberinput[i,j]=total(blur*(fibermap[*,*,j] eq i))
    endfor
  endfor

  ; Add crosstalk between adjacent fibers.
  ; Assume that they are ordered 1--->nfiber along the slit
  fiberinputCT=fiberinput
  if crosstalk gt 0.0 then begin
     if (keyword_set(verbose)and (k eq 0)) then print, "Adding slit/CCD crosstalk ", crosstalk 
     for fk=0,nfiber-1 do begin
        ; Special case for first fiber
        if (fk eq 0) then begin
           fiberinputCT[0,*]+=crosstalk*fiberinput[1,*]   ; add cross talk
           fiberinputCT[1,*]-=crosstalk*fiberinput[1,*]   ; subtract the flux that went into other fibers
        endif
        ; Fiber 2 to n-1
        if ((fk ne 0)and(fk ne nfiber-1)) then begin
           fiberinputCT[fk,*]+=crosstalk*(fiberinput[fk-1,*]+fiberinput[fk+1,*]) ; add cross-talk
           fiberinputCT[fk-1,*]-=crosstalk*(fiberinput[fk-1,*])  ; subtract the flux that went into other fibers
           fiberinputCT[fk+1,*]-=crosstalk*(fiberinput[fk+1,*])  ; subtract the flux that went into other fibers
        endif
        ; Special case for last fiber
        if (fk eq nfiber-1) then begin
           fiberinputCT[nfiber-1,*]+=crosstalk*fiberinput[nfiber-2,*] ; add cross-talk
           fiberinputCT[nfiber-2,*]-=crosstalk*fiberinput[nfiber-2,*] ; subtract from originating fiber
        endif
     endfor
  endif

  ; Now fiberinput is the input fiber spectra WITHOUT crosstalk,
  ; and fiberinputCT is the input fiber spectra WITH crosstalk.
  ; Carry both forward through simulation.

  ; Populate CCDspec with data for reference.  This gets values *with* crosstalk
  for j=0,nexp-1 do begin
    for i=0,nfiber-1 do CCDspec[k,j*nfiber+i]=fiberinput[i,j]
  endfor

  ; Convert from flux units to e-/channel/s at this wavelength slice
  ; using BOSS calibration matrix (dimension 2)
  objphot=fiberinput*calibmatrix[2,k]
  objphotCT=fiberinputCT*calibmatrix[2,k]
  ; Convert to total e-/channel by multiplying by exposure time
  for j=0,nexp-1 do begin
    objphot[*,j]=objphot[*,j]*exptime[j]
    objphotCT[*,j]=objphotCT[*,j]*exptime[j]
  endfor

  ; Make a sky slice in e-/channel
  ; If the /flatsky keyword is set then use an uber-simple flat spectrum
  if (keyword_set(flatsky)) then begin
    for j=0,nexp-1 do skyphot[*,j]=exptime[j]*calibmatrix[1,k]*flatskyspectrum[k]
  ; Otherwise, use the BOSS sky model
  endif else begin
    for j=0,nexp-1 do skyphot[*,j]=exptime[j]*calibmatrix[1,k]*calibmatrix[3,k]
  endelse

  ; Determine read noise
  rn=2.5 ; read noise per pixel
  if (k eq 0) then angppix=wave[k+1]-wave[k] $
  else angppix=wave[k]-wave[k-1]
  ; Multiply by sqrt(3.54)=1.88 for the effective 1-D size of a Gaussian spot
  ; with sigma=1 pix, which is 3.54 pix.
  rn=rn*1.88/sqrt(angppix); Convert to read noise per angstrom

  ; Generate a realization of the noise array
  ; using RANDOMN
  seed=rseed+2+k; Reset RNG
  ; Zero out non-physical negative object fluxes
  unphys=where(objphot lt 0.)
  if (size(unphys))[0] ne 0 then objphot[unphys]=0.
  unphys=where(objphotCT lt 0.)
  if (size(unphys))[0] ne 0 then objphotCT[unphys]=0.
  ; Compute noise array from the sum of sky noise, read noise,
  ; and the object flux WITH crosstalk
  var=sqrt(objphotCT+skyphot+rn*rn)
  noisearr=RANDOMN(seed,nfiber,nexp)*var

  ; Recovered flux per fiber is object+sky+noise, minus sky.
  ; Assume that the sky substracted off is known much more accurately than
  ; the sky added to the object so that it subtracts off perfectly leaving
  ; just object+noise.  (I.e., sky noise only enters once).
  ; Divide by exposure time to make it counts/sec (accomodates variable
  ; exposure times this way), and divide out the instrumental response
  ; to flux calibrate the output
  for j=0,nexp-1 do begin
     ; If flag CTFLUX is set, assume pipeline CANNOT remove crosstalk
     ; signal from either the flux vector or the noise vector
     if (keyword_set(ctflux)) then recflux[*,j]=(objphotCT[*,j]+noisearr[*,j])/exptime[j]/calibmatrix[2,k] $
     ; Otherwise, assume pipeline can remove crosstalk from flux vector
     ; but not the noise vector
     else recflux[*,j]=(objphot[*,j]+noisearr[*,j])/exptime[j]/calibmatrix[2,k]
     ; Error vector is just the noise sources
     var[*,j]=(var[*,j])/exptime[j]/calibmatrix[2,k]
     recflux_err[*,j]=(noisearr[*,j])/exptime[j]/calibmatrix[2,k]
  endfor
  ; Flux calibration looks pretty good in tests- output has about
  ; 51% the total flux in any given slice of the input, which corresponds
  ; to about the amount expected to lose to fill factor

  ; In debug mode, make an 'image' for each exposure
  if (keyword_set(debug)) then begin
    IndivImages[*,*,*]=0.
    for j=0,nexp-1 do begin
      for i=0,nfiber-1 do begin
        IndivImages[*,*,j] += recflux[i,j]*(fibermap[*,*,j] eq i)/fibareaW; Flux conservation
      endfor
    endfor
    imfile=strcompress('working/IndivImages'+strtrim(k)+'.fits',/remove_all)
    writefits,imfile,IndivImages
  endif

  ; Populate RSS with data for reference
  for j=0,nexp-1 do begin
    for i=0,nfiber-1 do begin
       RSS[k,j*nfiber+i]=recflux[i,j]
       RSS_err[k,j*nfiber+i]=recflux_err[i,j]
    endfor
  endfor

  ; Re-order info into 1-d vectors for interpolation routine
  ; Make sure to rescale x,y locations to output grid scale
  for j=0,nexp-1 do begin
    for i=0,nfiber-1 do begin
      ; Rescale x and y locations
      ; If unknown bundle positioning errors were specified, add them here
      xout[i+j*nfiber]=(xmatrix[i,j,k]+rgenU[j,0])*pixscaleW/pixscaleO
      yout[i+j*nfiber]=(ymatrix[i,j,k]+rgenU[j,1])*pixscaleW/pixscaleO
      ; Adjust flux to different pixel scale
      fluxout[i+j*nfiber]=recflux[i,j];*(pixscaleW/pixscaleO)*(pixscaleW/pixscaleO)
      errout[i+j*nfiber]=recflux_err[i,j];*(pixscaleW/pixscaleO)*(pixscaleW/pixscaleO)
      varout[i+j*nfiber]=var[i,j];*(pixscaleW/pixscaleO)*(pixscaleW/pixscaleO)
      statout[i+j*nfiber]=fiber_status[i]
    endfor
  endfor

  ; Interpolate to a regular grid
  ; Empirically, using rlim=1.6'' and sigma=0.7 arcsec gives
  ; the best image quality for CALIFA method thus far.
  interpimage=mlcalifainterp(xout,yout,fluxout,statout,[XsizeOUT,YsizeOUT],1.6/pixscaleO,0.7/pixscaleO)
  ; Multiply as necessary to take out scale differences from
  ; working vs output pixel scale
  FinalCube[*,*,k]=interpimage;*(pixscaleO/pixscaleW)
  interpimage=mlcalifainterp(xout,yout,varout,statout,[XsizeOUT,YsizeOUT],1.6/pixscaleO,0.7/pixscaleO)
  FinalErrCube[*,*,k]=interpimage

endfor

; Make information header
mkhdr,head,FinalCube

fxaddpar,head,'CDELT1',pixscaleO
fxaddpar,head,'CDELT2',pixscaleO
fxaddpar,head,'CUNIT1','arcsec'
fxaddpar,head,'CUNIT2','arcsec'
fxaddpar,head,'CTYPE1','LINEAR'
fxaddpar,head,'CTYPE2','LINEAR'

fxaddpar,head,'CUNIT3','Angstroms'
fxaddpar,head,'CRVAL3',alog10(wave[nwave-1])
fxaddpar,head,'CD3_3',alog10(wave[1])-alog10(wave[0])
fxaddpar,head,'CRPIX3',nwave
fxaddpar,head,'DC-FLAG',1
fxaddpar,head,'CTYPE3','LOG10'
fxaddpar,head,'DISPAXIS',3

fxaddpar,head,'UNITS','1e-17 erg/s/cm^2/Angstrom'

writefits,CCDfile,CCDspec
writefits,RSSfile,RSS
writefits,RSSfile_err,RSS_err
mwrfits,FinalCube,OutpFile,head,/create
mwrfits,FinalErrCube,OutpFile

return
end

  
