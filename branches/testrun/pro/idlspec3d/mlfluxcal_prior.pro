;+
; function mlfluxcal_prior
;
; As mlfluxcal, but doesn't construct new flux calibration vectors,
; uses previously determined vectors that have been archived.
; Somewhat kludgy for test-run data only.
;
;-

;------------------------------------------------------------------------------

function makelabel, hdr

   camera = strtrim(sxpar(hdr, 'CAMERAS'),2)
   expos =  strtrim(string(sxpar(hdr, 'EXPOSURE')),2)
   flat  =  strmid(sxpar(hdr, 'FLATFILE'),7,8)
   arc   =  strmid(sxpar(hdr, 'ARCFILE'),7,8)

   label = string(camera, expos, flat, arc, $
    format='(a2,"-",i8.8,"-",a8,"-",a8)')

   return, label
end

;------------------------------------------------------------------------------
pro add_iraf_keywords, hdr, wavemin, binsz

   sxaddpar, hdr, 'WAT0_001', 'system=linear'
   sxaddpar, hdr, 'WAT1_001', $
    'wtype=linear label=Wavelength units=Angstroms'
   sxaddpar, hdr, 'CRVAL1', wavemin, $
    ' Central wavelength (log10) of first pixel'
   sxaddpar, hdr, 'CD1_1', binsz, ' Log10 dispersion per pixel'
   sxaddpar, hdr, 'CRPIX1', 1, ' Starting pixel (1-indexed)'
   sxaddpar, hdr, 'CTYPE1', 'LINEAR'
   sxaddpar, hdr, 'DC-FLAG', 1, ' Log-linear flag'

   return
end


;------------------------------------------------------------------------------


function mlfluxcal_prior,fileb1,filer1,Cfile1,obsparam,slitmap,waveset,plotfile=plotfile

; Status starts off OK
status=0

; Set up plot file if specified, otherwise use display
if (keyword_set(plotfile)) then begin
  set_plot,'ps'
  device,filename=plotfile,/color
  splog,'Printing to ',plotfile
endif else set_plot,'x'

splog,'Full reduction set- assuming spSFrame inputs.'
splog,'Reading science + wavelength data from files:'
splog,fileb1
splog,filer1
splog,'NO good calibration info: using archived flux cal data'

adderr=0.03
minfracthresh=0.80
objname=strarr(2)
objname[0]=fileb1
objname[1]=filer1

; Populate some variables with header information
nfile=n_elements(objname)
plateid=lonarr(nfile)
mjd=lonarr(nfile)
camname=strarr(nfile); b1,r1,b2, or r2
expnum=lonarr(nfile)
spectroid=lonarr(nfile); Spectrograph id (1 or 2)
npixarr=lonarr(nfile); Number of spectral pixels

for ifile=0,nfile-1 do begin
  hdr=headfits(objname[ifile])
  plateid[ifile] = strtrim(sxpar(hdr, 'PLATEID'),2)
  mjd[ifile] = strtrim(sxpar(hdr, 'MJD'),2)
  camname[ifile] = strtrim(sxpar(hdr, 'CAMERAS'),2)
  spectroid[ifile] = strmid(camname[ifile],1,1)
  expnum[ifile] = sxpar(hdr, 'EXPOSURE')
  npixarr[ifile] = sxpar(hdr, 'NAXIS1')
endfor
maxmjd = max(mjd)

binsz = 1.0d-4
zeropoint = 3.5D
window = 100

npix = max(npixarr); Number of spectral pixels
nfiber = n_elements(slitmap) ; Number of fibers in one spectrograph

; Sort filenames such that this list contains first the blue then the red
spframes=objname
nfiles = n_elements(spframes)
filenames = spframes[sort(spframes)]

camnames=camname
ncam = N_elements(camnames)
nexpvec = lonarr(ncam)
exptimevec = fltarr(ncam)

; Loop through each 2D output and read in the data

; Start by determining the size of all the files
npixarr = lonarr(nfiles)
for ifile=0, nfiles-1 do begin
  mlframe_read, filenames[ifile], hdr=objhdr
  npixarr[ifile] = sxpar(objhdr,'NAXIS1')
endfor
npixmax = max(npixarr)
nobj = sxpar(objhdr,'NAXIS2') ; Number of fibers per spectrograph
; DRL unhappy with this...  Why?

for ifile=0, nfiles-1 do begin
  ; Read in all data from this input file.
  ; Reading the plug-map structure will fail if its structure is
  ; different between different files.

  splog, 'Reading file #', ifile, ': ', filenames[ifile]
  mlframe_read, filenames[ifile], objflux=tempflux, objivar=tempivar, $
    mask=temppixmask, wset=tempwset, dispset=tempdispset, $
    skyflux=tempsky, $
    hdr=hdr, adderr=adderr

  if (ifile EQ 0) then $
    hdrarr = ptr_new(hdr) $
  else $
    hdrarr = [hdrarr, ptr_new(hdr)]

  ; Read header info
  thismjd = sxpar(hdr, 'MJD')
  if (NOT keyword_set(mjdlist)) then mjdlist = thismjd $
  else mjdlist = [mjdlist, thismjd]
  cameras = strtrim(sxpar(hdr, 'CAMERAS'),2)
  expstr = string(sxpar(hdr, 'EXPOSURE'), format='(i8.8)')

  ; Solve for wavelength and lambda-dispersion at each pixel in the image
  ;need to reset junk since the array lengths change
  junk=0
  traceset2xy, tempwset, junk, tempwave
  traceset2xy, tempdispset, junk, tempdispersion

  ; Here is the correct conversion from pixels to log-lambda dispersion.
  ; We are converting from the dispersion in units of spFrame pixel sizes
  ; to the dispersion in units of the new rebinned pixel size, which is
  ; BINSZ in log-lambda units.
      
  ; this probably should be fixed elsewhere but limit where the fit range
  tempxmax=tempwset.xmax
  tempwset.xmax=(size(tempdispersion,/dimens))[0]-1
  correct_dlam, tempdispersion, 0, tempwset, dlam=binsz, /inverse
  tempwset.xmax=tempxmax

  dims = size(tempflux, /dimens)
  npix = dims[0]
  nfib = dims[1]

  ; Make a map of the size of each pixel in delta-(log10-Angstroms),
  ; and re-normalize the flux to electrons/(dloglam)
  correct_dlam, tempflux, tempivar, tempwset, dlam=binsz
  correct_dlam, tempsky, 0, tempwset, dlam=binsz

  ; Determine if this is a blue or red spectrum
  icam = (where(cameras EQ camnames))[0]
  if (icam EQ -1) then $
    message, 'Unknown camera ' + cameras
  nexpvec[icam] = nexpvec[icam] + 1
  exptimevec[icam] = exptimevec[icam] + sxpar(hdr, 'EXPTIME')

  ; Apply spectro-photometric calibration
  expnum = sxpar(hdr, 'EXPOSURE')

  if(camnames[icam] eq 'b1') then begin
    calibfile=getenv('MANGACORE_DIR')+'/'+getenv('MANGACORE_VER')+'/bosscal/56280/spFluxcalib-b1-00154443.fits.gz'
  endif
  if(camnames[icam] eq 'r1') then begin
    calibfile=getenv('MANGACORE_DIR')+'/'+getenv('MANGACORE_VER')+'/bosscal/56280/spFluxcalib-r1-00154443.fits.gz'
  endif

  if (keyword_set(calibfile)) then begin
    calibfac = mrdfits(calibfile, 0, calibhdr, /silent)
  endif else begin
    splog, 'WARNING: Reading default flux-calib vectors for camera=' $
       + camnames[icam]
    calibfac = fcalib_default(camnames[icam], tempwave, exptimevec[icam])
  endelse
  minval = 0.05 * mean(calibfac)

  divideflat, tempflux, invvar=tempivar, calibfac, minval=minval
  divideflat, tempsky, calibfac, minval=minval
  temppixmask = temppixmask $
    OR ((calibfac LE minval OR keyword_set(calibfile) EQ 0) $
    * pixelmask_bits('BADFLUXFACTOR'))

  ; This is ordinarily where BOSS would do a flux correction factor to 
  ; account for DAR, among other things.  I am NOT doing this.
  ; It's possible that I need to do part of it though...

  ; Apodize the errors
  ; Do this only for the dichroic overlap region, which are the first
  ; rows in both the blue and red CCD's.

  if (keyword_set(window)) then begin
    swin = window < npix
    indx = lindgen(swin)
    tempivar[indx,*] = $
    tempivar[indx,*] * (indx # replicate(1,nfib)) / swin
  endif

  ; Concatenate data from all images
  if (ifile EQ 0) then begin
    ; Construct the image arrays
    flux = make_array(npixmax,nobj*nfiles,type=size(tempflux,/type))
    fluxivar = make_array(npixmax,nobj*nfiles,type=size(tempivar,/type))
    wave = make_array(npixmax,nobj*nfiles,type=size(tempwave,/type))
    dispersion = make_array(npixmax,nobj*nfiles,type=size(tempdisp,/type))
    pixelmask = make_array(npixmax,nobj*nfiles,type=size(temppixmask,/type))
    skyflux = make_array(npixmax,nobj*nfiles,type=size(tempsky,/type))
    ;ximg = make_array(npixmax,nobj*nfiles,type=size(tempximg,/type))
    ;superflat = make_array(npixmax,nobj*nfiles,type=size(tempsuperflat,/type))

    ; Append as vectors...
    camerasvec = cameras
    label = makelabel(hdr)
    filenum = lonarr(nfib) + ifile
  endif else begin
    ; Append as vectors...
    camerasvec = [camerasvec, cameras]
    label = [label, makelabel(hdr)]
    filenum = [filenum, lonarr(nfib) + ifile]
  endelse

  flux[0:npixarr[ifile]-1,nobj*ifile:nobj*(ifile+1)-1] = tempflux
  fluxivar[0:npixarr[ifile]-1,nobj*ifile:nobj*(ifile+1)-1] = tempivar
  wave[0:npixarr[ifile]-1,nobj*ifile:nobj*(ifile+1)-1] = tempwave
  ; Pad the wavelengths with reasonable values
  if (npixarr[ifile] LT npixmax) then begin
    dwave = tempwave[npixarr[ifile]-1,*] - tempwave[npixarr[ifile]-2,*]
    addwave = tempwave[npixarr[ifile]-1,*] $
      ## (1+lonarr(npixmax-npixarr[ifile])) $
      + dwave ## (1+lindgen(npixmax-npixarr[ifile]))
    wave[npixarr[ifile]:npixmax-1,nobj*ifile:nobj*(ifile+1)-1] = addwave
  endif
  dispersion[0:npixarr[ifile]-1,nobj*ifile:nobj*(ifile+1)-1] = tempdispersion
  pixelmask[0:npixarr[ifile]-1,nobj*ifile:nobj*(ifile+1)-1] = temppixmask
  skyflux[0:npixarr[ifile]-1,nobj*ifile:nobj*(ifile+1)-1] = tempsky
  ;ximg[0:npixarr[ifile]-1,nobj*ifile:nobj*(ifile+1)-1] = tempximg
  ;superflat[0:npixarr[ifile]-1,nobj*ifile:nobj*(ifile+1)-1] = tempsuperflat
endfor

tempflux = 0
tempivar = 0
tempwave = 0
tempdispersion = 0
temppixmask = 0
tempsky = 0

; Check how many exposures we have in each of the (4) cameras

for icam=0, ncam-1 do begin
  junk = where(camerasvec EQ camnames[icam], nmatch)
  splog, 'Files for camera ' + camnames[icam] + ':', nmatch
  if (icam EQ 0) then nminfile = nmatch $
  else nminfile = nminfile < nmatch
endfor
; ??? Should make this routine robust to fewer files!!!
if (nminfile LT 1) then begin
  splog, 'ABORT: At least 1 file needed for each camera'
  status=-143L
  mlquitmanga3d,status
endif

; Construct output data structures, including the wavelength scale
totalpix = (size(flux, /dimens))[0]
nonzero = where(fluxivar GT 0.0)
minfullwave = min(wave[nonzero])
maxfullwave = max(wave[nonzero])

; Get max and min wavelength from good pixels
if (NOT keyword_set(wavemin)) then begin
  spotmin = long((minfullwave - zeropoint)/binsz) + 1L
  spotmax = long((maxfullwave - zeropoint)/binsz)
  wavemin = spotmin * binsz + zeropoint
  wavemax = spotmax * binsz + zeropoint
endif else begin
  spotmin = 0L
  if (NOT keyword_set(wavemax)) then begin
    spotmax = long((maxfullwave - wavemin)/binsz)
    wavemax = spotmax * binsz + wavemin
  endif else spotmax = long((wavemax - wavemin)/binsz)
endelse

; This is where output wavelengths are set
; kludge to be FIXED to BOSSCAL wavelengths
; so that everything is on same grid without
; requiring multiple resamplings
;nfinalpix = spotmax - spotmin + 1L
;finalwave = dindgen(nfinalpix) * binsz + wavemin
nfinalpix = (size(waveset))[1]
finalwave = dindgen(nfinalpix) * binsz + alog10(waveset[0])

; This seems to be BOSS scaling up from one detector to two??
;nfiber = 2 * nfib

finalflux = fltarr(nfinalpix, nfiber)
finalivar = fltarr(nfinalpix, nfiber)
finalandmask = lonarr(nfinalpix, nfiber)
finalormask = lonarr(nfinalpix, nfiber)
finaldispersion = fltarr(nfinalpix, nfiber)
finalsky = fltarr(nfinalpix, nfiber)

; Note that currently this is combining ALL fibers

; Combine each fiber, one at a time
for ifiber=0, nfiber-1 do begin
  ; Find the first occurance of fiber number IFIBER+1
  indx = (where(slitmap.fiberid EQ ifiber+1))[0]

  if (indx NE -1) then begin
    splog, 'Fiber', ifiber+1, ' ', slitmap[indx].ifuname, $
      slitmap[indx].mag, format = '(a, i5.4, a, a, f6.2, 5f6.2)'

    ; Identify all objects with the same XFOCAL,YFOCAL plate position, and
    ; combine all these objects into a single spectrum.
    ; If all pluggings are identical, then this will always be
    ; the same fiber ID.
    ; Also, insist that the object type is not 'NA', which would
    ; occur for unplugged fibers. <--- Disable this for BOSS ???

    ; DRL- comment this out, do it for ALL fibers.
    ; Ah- this is how we combine b+r, is by finding both locations
    ; where a fiber is in both b and r rows.
    ; WEIRD- why does BOSS do it this v awkward way???
    indx=[where(slitmap.fiberid eq ifiber+1),where(slitmap.fiberid eq ifiber+1)+nfiber]
    if (ifiber eq 0) then splog,indx
;    indx = where(abs(slitmap.xfocal - plugmap[indx].xfocal) LT 0.0001 $
 ;     AND abs(plugmap.yfocal - plugmap[indx].yfocal) LT 0.0001)
  endif

  if (indx[0] NE -1) then begin
    temppixmask = pixelmask[*,indx]
    combine1fiber, wave[*,indx], flux[*,indx], fluxivar[*,indx], $
      finalmask=temppixmask, indisp=dispersion[*,indx], $
      skyflux=skyflux[*,indx], $
      newloglam=finalwave, newflux=bestflux, newivar=bestivar, $
      andmask=bestandmask, ormask=bestormask, newdisp=bestdispersion, $
      newsky=bestsky, $
      nord=nord, binsz=binsz, bkptbin=bkptbin, maxsep=maxsep, $
      maxiter=50, upper=3.0, lower=3.0, maxrej=1

    finalflux[*,ifiber] = bestflux
    finalivar[*,ifiber] = bestivar
    finalandmask[*,ifiber] = bestandmask
    finalormask[*,ifiber] = bestormask
    finaldispersion[*,ifiber] = bestdispersion
    finalsky[*,ifiber] = bestsky

    ; DRL- a really stupid kludge to take out wacky values
    ; in the far blue.  Should do something better.
    for j=0,100 do begin
      if abs(finalflux[j,ifiber]) gt 100. then finalflux[j,ifiber]=0.
    endfor

    ; The following adds the COMBINEREJ bit to the input pixel masks
    pixelmask[*,indx] = temppixmask
  endif else begin
    splog, 'Fiber', ifiber+1, ' NO DATA'
    finalandmask[*,ifiber] = pixelmask_bits('NODATA')
    finalormask[*,ifiber] = pixelmask_bits('NODATA')
  endelse
endfor


   ; Modify the 1st file's header to use for the combined plate header.

   bighdr = *hdrarr[0]

; BOSS does a 'flux distortion' correction here, do I need to???


  ;Indiv calib frames?  DRL should remove


;   for ifile=0, nfiles-1 do begin
;      thisfile = fileandpath(filenames[ifile], path=thispath)
;      thisfile = djs_filepath(repstr(thisfile,'spSFrame','spTFrame'), $
;       root_dir=thispath)
;      splog, 'Writing file #', ifile, ': ', thisfile, prename=filenames[ifile]
;      indx = where(filenum EQ ifile, nthis)
;
;      hdr = *hdrarr[ifile]
;      sxaddpar, hdr, 'BUNIT', '1E-17 erg/cm^2/s/Ang'
;
;      ; Apply the flux-distortion image to each individual frame, by
;      ; interpolating off the full wavelength-scale distortion image
;      ; onto the wavelength mapping of each individual exposure+CCD.
;
;      for i=0L, nthis-1 do begin
;         thisflux1 = flux[*,indx[i]]
;         thisivar1 = fluxivar[*,indx[i]]
;         thissky1 = skyflux[*,indx[i]]
;         ;;j = slitmap[indx[i]].fiberid - 1
;         ;;thisicorr = interpol(invcorrimg[*,j], finalwave, wave[*,indx[i]])
;         ;;divideflat, thisflux1, invvar=thisivar1, thisicorr, minval=minicorrval
;         flux[*,indx[i]] = thisflux1
;         fluxivar[*,indx[i]] = thisivar1
;
;         ;divideflat, thissky1, thisicorr, minval=minicorrval
;         skyflux[*,indx[i]] = thissky1
;      endfor
;
;      mwrfits, flux[*,indx], thisfile, hdr, /create
;      mwrfits, fluxivar[*,indx], thisfile
;      mwrfits, pixelmask[*,indx], thisfile
;      mwrfits, wave[*,indx], thisfile
;      mwrfits, dispersion[*,indx], thisfile
;      mwrfits, slitmap[indx], thisfile
;      mwrfits, skyflux[*,indx], thisfile
;      mwrfits, ximg[*,indx], thisfile
;      mwrfits, superflat[*,indx], thisfile
;   endfor
;   splog, prename=''

   ;----------
   ; Clear memory

   wave = 0
   flux = 0
   fluxivar = 0
   temppixmask = 0
   dispersion = 0
   skyflux = 0
   superflat = 0



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Write the corrected spCFrame files
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   ;---------------------------------------------------------------------------
   ; Create the output header
   ;---------------------------------------------------------------------------

   ;----------
   ; Print roll call of bad fibers and bad pixels.

   fiber_rollcall, finalandmask, finalwave

   ;----------
   ; Remove header cards that were specific to this first exposure
   ; (where we got the header).

   ncoeff = sxpar(bighdr, 'NWORDER')
   for i=2, ncoeff-1 do sxdelpar, bighdr, 'COEFF'+strtrim(string(i),2)

   sxdelpar, bighdr, ['SPA', 'IPA', 'IPARATE']
   sxdelpar, bighdr, 'EXPOSURE'
   sxdelpar, bighdr, 'REQTIME'
   sxdelpar, bighdr, 'QUALITY'
   sxdelpar, bighdr, 'FILENAME'
   sxdelpar, bighdr, 'SEQID'
   sxdelpar, bighdr, 'DARKTIME'
   sxdelpar, bighdr, 'CAMERAS'
   sxdelpar, bighdr, 'PLUGMAPO'
   for i=1, 4 do sxdelpar, bighdr, 'GAIN'+strtrim(string(i),2)
   for i=1, 4 do sxdelpar, bighdr, 'RDNOISE'+strtrim(string(i),2)
   sxdelpar, bighdr, ['CAMCOL', 'CAMROW']
   sxdelpar, bighdr, ['AMPLL', 'AMPLR', 'AMPUL', 'AMPUR']
   sxdelpar, bighdr, ['FFS', 'FF', 'NE', 'HGCD']
   sxdelpar, bighdr, ['SPEC1', 'SPEC2']
   sxdelpar, bighdr, 'NBLEAD'
   sxdelpar, bighdr, 'PIXFLAT'
   sxdelpar, bighdr, 'PIXBIAS'
   sxdelpar, bighdr, 'FLATFILE'
   sxdelpar, bighdr, 'ARCFILE'
   sxdelpar, bighdr, 'OBJFILE'
   sxdelpar, bighdr, 'FRAMESN2'
   sxdelpar, bighdr, 'DEREDSN2'

   ;----------
   ; Average together some of the fields from the individual headers.

   cardname = [ 'AZ', 'ALT', 'TAI', 'WTIME', 'AIRTEMP', 'DEWPOINT', $
    'DEWDEP', 'DUSTA', 'DUSTB', 'DUSTC', 'DUSTD', 'GUSTS', 'HUMIDITY', $
    'HUMIDOUT', 'PRESSURE', 'WINDD', 'WINDS', 'TEMP01', 'TEMP02', $
    'TEMP03', 'TEMP04', 'HELIO_RV', 'SEEING20', 'SEEING50', 'SEEING80', $
    'RMSOFF20', 'RMSOFF50', 'RMSOFF80', 'XCHI2', 'SKYCHI2', $
    'WSIGMA', 'XSIGMA' ]
   sxcombinepar, hdrarr, cardname, bighdr, func='average'

   sxcombinepar, hdrarr, 'TAI-BEG', bighdr, func='min'
   sxcombinepar, hdrarr, 'TAI-END', bighdr, func='max'

   sxcombinepar, hdrarr, 'XCHI2', bighdr, func='max', outcard='XCHI2MAX', $
    after='XCHI2'
   sxcombinepar, hdrarr, 'XCHI2', bighdr, func='min', outcard='XCHI2MIN', $
    after='XCHI2'

   sxcombinepar, hdrarr, 'SKYCHI2', bighdr, func='max', outcard='SCHI2MAX', $
    after='SKYCHI2'
   sxcombinepar, hdrarr, 'SKYCHI2', bighdr, func='min', outcard='SCHI2MIN', $
    after='SKYCHI2'

   sxcombinepar, hdrarr, 'WSIGMA', bighdr, func='max', outcard='WSIGMAX', $
    after='WSIGMA'
   sxcombinepar, hdrarr, 'WSIGMA', bighdr, func='min', outcard='WSIGMIN', $
    after='WSIGMA'

   sxcombinepar, hdrarr, 'XSIGMA', bighdr, func='max', outcard='XSIGMAX', $
    after='XSIGMA'
   sxcombinepar, hdrarr, 'XSIGMA', bighdr, func='min', outcard='XSIGMIN', $
    after='XSIGMA'

   ; Add the NGUIDE keywords for all headers of one flavor of CAMERAS
   ; (e.g., for all the 'b1' exposures if the first frame is 'b1'.)
   cardname = 'NGUIDE'
   sxcombinepar, hdrarr[0], cardname, bighdr, func='total'
   cameras0 = sxpar(*(hdrarr[0]), 'CAMERAS')
   for ihdr=1, n_elements(hdrarr)-1 do begin
      if (sxpar(*(hdrarr[ihdr]), 'CAMERAS') EQ cameras0) then $
       sxcombinepar, hdrarr[ihdr], cardname, bighdr, func='total'
   endfor

   ;----------
   ; Use the MJD passed as a keyword, which will typically be for the most
   ; observation, and be consistent with the output file names

   ; DRL for some reason this was vector??? Kludge.
   if (keyword_set(mjd)) then $
    sxaddpar, bighdr, 'MJD', obsparam.mjd

   ; Get the list of MJD's used for these reductions, then convert to a string
   mjdlist = mjdlist[uniq(mjdlist, sort(mjdlist))]
   mjdlist = strtrim(strcompress(string(mjdlist,format='(99a)')),2)
   sxaddpar, bighdr, 'MJDLIST', mjdlist, after='MJD'

   ;----------
   ; Add new header cards

   sxaddpar, bighdr, 'VERSCOMB', idlspec2d_version(), $
    ' Version of idlspec2d for combining multiple spectra', after='VERS2D'
   sxaddpar, bighdr, 'NEXP', nfiles, $
    ' Number of exposures in this file', before='EXPTIME'
   for ifile=0,nfiles-1 do $
    sxaddpar, bighdr, string('EXPID',ifile+1, format='(a5,i2.2)'), label[ifile], $
     ' ID string for exposure '+strtrim(ifile+1,2), before='EXPTIME'
   if (keyword_set(bestexpnum)) then $
    sxaddpar, bighdr, 'BESTEXP', bestexpnum, before='EXPID01'

   sxaddpar, bighdr, 'EXPTIME', min(exptimevec), $
    ' Minimum of exposure times for all cameras'
   for icam=0, ncam-1 do $
    sxaddpar, bighdr, 'NEXP_'+camnames[icam], nexpvec[icam], $
     ' '+camnames[icam]+' camera number of exposures', before='EXPTIME'
   for icam=0, ncam-1 do $
    sxaddpar, bighdr, 'EXPT_'+camnames[icam], exptimevec[icam], $
     ' '+camnames[icam]+' camera exposure time (seconds)', before='EXPTIME'
   sxaddpar, bighdr, 'SPCOADD', systime(), $
    ' SPCOADD finished', after='EXPTIME'

   sxaddpar, bighdr, 'NWORDER', 2, ' Linear-log10 coefficients'
   sxaddpar, bighdr, 'NWORDER', 2, ' Linear-log10 coefficients'
   sxaddpar, bighdr, 'WFITTYPE', 'LOG-LINEAR', ' Linear-log10 dispersion'
   sxaddpar, bighdr, 'COEFF0', wavemin, $
    ' Central wavelength (log10) of first pixel'
   sxaddpar, bighdr, 'COEFF1', binsz, ' Log10 dispersion per pixel'

   sxaddpar, bighdr, 'NAXIS1', n_elements(bestflux)
   sxaddpar, bighdr, 'NAXIS2', nfiber

   spawn, 'uname -n', uname
   sxaddpar, bighdr, 'UNAME', uname[0]

   ;----------
   ; Compute the fraction of bad pixels in total, and on each spectrograph.
   ; Bad pixels are any with SKYMASK(INVVAR)=0, excluding those where
   ; the NODATA bit is set in the pixel mask.

   ifib1 = where(slitmap.spectrographid EQ 1, nfib1)
   ifib2 = where(slitmap.spectrographid EQ 2, nfib2)
   qbadpix = skymask(finalivar, finalandmask, finalormask) EQ 0 $
    AND (finalandmask AND pixelmask_bits('NODATA')) EQ 0
   if (nfib1 GT 0) then $
    fbadpix1 = total(qbadpix[*,ifib1]) / (nfib1 * nfinalpix) $
   else $
    fbadpix1 = 0
   if (nfib2 GT 0) then $
    fbadpix2 = total(qbadpix[*,ifib2]) / (nfib2 * nfinalpix) $
   else $
    fbadpix2 = 0
   if (nfib1 GT 0 AND nfib2 GT 0) then $
    fbadpix = total(qbadpix[*,[ifib1,ifib2]]) / ((nfib1+nfib2) * nfinalpix) $
   else if (nfib1 GT 0) then $
    fbadpix = fbadpix1 $
   else if (nfib2 GT 0) then $
    fbadpix = fbadpix1 $
   else $
    fbadpix = 0

   sxaddpar, bighdr, 'FBADPIX', fbadpix, ' Fraction of bad pixels'
   sxaddpar, bighdr, 'FBADPIX1', fbadpix1, ' Fraction of bad pixels on spectro-1'
   sxaddpar, bighdr, 'FBADPIX2', fbadpix2, ' Fraction of bad pixels on spectro-2'

   ;----------
   ; Add keywords for IRAF-compatability

   add_iraf_keywords, bighdr, wavemin, binsz

   mkhdr, hdrfloat, finalivar, /image, /extend
   add_iraf_keywords, hdrfloat, wavemin, binsz

   mkhdr, hdrlong, finalandmask, /image, /extend
   add_iraf_keywords, hdrlong, wavemin, binsz

   ; First write the file with the flux distortion vectors
   ; DRL- didn't use this...
   ;mwrfits, corrimg, distortfitsfile, bighdr, /create

   outname = Cfile1
   len=strlen(outname)
   suffix=strmid(outname,len-3,3)
   ; If filename ended in .gz suffix, remove it.
   ; Will be added back when we gzip later
   if (suffix eq '.gz') then outname=strmid(outname,0,len-3)

   ; HDU #0 is flux
   sxaddpar, bighdr, 'BUNIT', '1E-17 erg/cm^2/s/Ang'
   mwrfits, finalflux, outname, bighdr, /create

   ; HDU #1 is inverse variance
   sxaddpar, hdrfloat, 'BUNIT', '1/(1E-17 erg/cm^2/s/Ang)^2'
   mwrfits, finalivar, outname, hdrfloat

   ; HDU #2 is AND-pixelmask
   mwrfits, finalandmask, outname, hdrlong

   ; HDU #3 is OR-pixelmask
   mwrfits, finalormask, outname, hdrlong

   ; HDU #4 is dispersion map
   sxaddpar, hdrfloat, 'BUNIT', 'pixels'
   mwrfits, finaldispersion, outname, hdrfloat

   ; HDU #5 is slitmap
   mwrfits, slitmap, outname

   ; HDU #6 is the sky
   mwrfits, finalsky, outname

   ; HDU #7 is the wavelength solution
   mwrfits, finalwave, outname

; Calculate output flux in AB mag units
finalfluxAB=finalflux
for i=0,nfinalpix-1 do begin
  lam=10.^finalwave[i]
  finalfluxAB[i,*]=-2.5*alog10(abs(finalflux[i,*])*lam*lam/3.14)+40.09
endfor

   ; HDU #8 is final flux in AB units
   mwrfits, finalfluxAB, outname

; gzip output file
spawn, ['gzip', '-f', outname], /noshell

; Close out the plot file
if (keyword_set(plotfile)) then device,/close

return,status
end

