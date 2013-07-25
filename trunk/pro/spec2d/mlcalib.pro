;+
; NAME:
;   mlcalib
;
; PURPOSE:
;   Extract calibration frames.
;
; CALLING SEQUENCE:
;   mlcalib, flatname, arcname, fibermask=, cartid=, $
;    lampfile=, indir=, timesep=, ecalibfile=, plottitle=, $
;    minflat=, maxflat=, arcinfoname=, flatinfoname=, $
;    arcstruct=, flatstruct=, writeflatmodel=, /bbspec]
;
; INPUTS:
;   flatname   - Name(s) of flat-field SDSS image(s)
;   arcname    - Name(s) of arc SDSS image(s)
;   cartid     - Cartridge ID from plugmap
;
; OPTIONAL KEYWORDS:
;   fibermask  - Fiber status bits, set nonzero for bad status [NFIBER].
;                Note this is not modified, but modified copies appear
;                in the returned structures ARCSTRUCT and FLATSTRUCT.
;   lampfile   - Name of file describing arc lamp lines, which would
;                over-ride the default file read by FITARCIMAGE.
;   indir      - Input directory for FLATNAME, ARCNAME, OBJNAME;
;                default to '.'
;   timesep    - Maximum time separation between flats and arcs to pair them;
;                set to zero to disable this test; default to 7200 sec.
;   ecalibfile - opECalib file to pass to SDSSPROC
;   plottitle  - Prefix for titles in QA plots.
;   minflat    - Parameter for SDSSPROC for pixel flats; default to 0.8
;   maxflat    - Parameter for SDSSPROC for pixel flats; default to 1.2
;   arcinfoname- File name (with path) to output arc extraction and fitting
;                information
;   flatinfoname-File name (with path) to output flat field extraction and
;                fitting information
; writeflatmodel-Set this keyword to write flat data image, ivar, and
;                final extraction model image to a file.  Will only
;                work if "flatinfoname" is present also (ASB).
; writearcmodel- Set this keyword to write arc data image, ivar, and
;                final extraction model image to a file.  Will only
;                work if "arcinfoname" is present also (ASB).
;   bbspec         - use bbspec extraction code
;
; OUTPUTS:
;   arcstruct  - Structure array with extracted arc calibration information
;   flatstruct - Structure array with extracted flat calibration information
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   Always pair arcs to the nearest good flat, and flats to the nearest good arc
;   (nearest in time, as defined by the TAI-BEG keyword in the FITS headers).
;
;   Also store SUPERFLATSET from fiberflat, since we need this to remove
;   small scale features present in all spectra

; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;   extract_image
;   fiberflat()
;   fitarcimage
;   fitdispersion
;   fitflatwidth()
;   get_tai
;   reject_arc()
;   reject_flat()
;   sdssproc
;   splog
;   trace320crude()
;   traceset2xy
;   xy2traceset
;
; INTERNAL SUPPORT ROUTINES:
;   create_arcstruct()
;   create_flatstruct()
;
; REVISION HISTORY:
;   24-Jan-2000  Written by D. Schlegel, Princeton
;   27-Nov-2000  Changed to proftype 3, added minflat, maxflat keywords
;    8-Jan-2001  And now back to proftype 1, more robust against bad columns
;   26-Jan-2001  And now let's check both 1&3, and use the better fit
;      Apr-2010  Added "write[flat,arc]model" option (A. Bolton, Utah)
;   25-Jan-2011  Added "twophase" test and switching, A. Bolton, Utah
;   29-Mar-2011  Switched to bundle-wise pure IDL extraction, A. Bolton, Utah
;
;-
;------------------------------------------------------------------------------
function create_arcstruct, narc

   on_error, 0
   compile_opt idl2
   compile_opt idl2, hidden
   
  ftemp = create_struct( name='ARC_STRUCT', $
    'NAME', '', $
    'TAI', 0D, $
    'TSEP', 0D, $
    'QBAD', 0B, $
    'IFLAT', -1, $
    'BESTCORR', 0.0, $
    'NMATCH', 0L, $
    'MEDWIDTH', fltarr(4), $
    'LAMBDA', ptr_new(), $
    'REJLINE', ptr_new(), $
    'XPEAK', ptr_new(), $
    'XDIF_TSET', ptr_new(), $
    'WSET', ptr_new(), $
    'DISPSET', ptr_new(), $
    'FIBERMASK', ptr_new() )

  arcstruct = replicate(ftemp, narc)
  
  return, arcstruct
end
;------------------------------------------------------------------------------
function create_flatstruct, nflat

   on_error, 0
   compile_opt idl2
   compile_opt idl2, hidden
   
  ftemp = create_struct( name='FLAT_STRUCT', $
    'NAME', '', $
    'TAI', 0D, $
    'TSEP', 0D, $
    'QBAD', 0, $
    'IARC', -1, $
    'PROFTYPE', 0, $
    'MEDWIDTH', fltarr(4), $
    'FIBERMASK', ptr_new(), $
    'TSET', ptr_new(), $
    'XSOL', ptr_new(), $
    'WIDTHSET', ptr_new(), $
    'FFLAT', ptr_new(), $
    'SUPERFLATSET', ptr_new() )
    
  flatstruct = replicate(ftemp, nflat)
  
  return, flatstruct
end
;------------------------------------------------------------------------------

pro mlcalib, flatname, arcname, fibermask=fibermask, cartid=cartid, $
    lampfile=lampfile, indir=indir, timesep=timesep, $
    ecalibfile=ecalibfile, plottitle=plottitle, $
    arcinfoname=arcinfoname, flatinfoname=flatinfoname, $
    arcstruct=arcstruct, flatstruct=flatstruct, skipproc=skipproc, $
    minflat=minflat, maxflat=maxflat, bbspec=bbspec, $
    writeflatmodel=writeflatmodel, writearcmodel=writearcmodel, $
    visual=visual, fiberparam=fiberparam, flat=flat, arc=arc

   on_error, 0
   compile_opt idl2
    
  if (~keyword_set(indir)) then indir = '.'
  if (~keyword_set(timesep)) then timesep = 7200
  if (~keyword_set(minflat)) then minflat = 0.8
  if (~keyword_set(maxflat)) then maxflat = 1.2
  
  stime1 = systime(1)
  
  ;---------------------------------------------------------------------------
  ; Determine spectrograph ID and color from first flat file
  ;---------------------------------------------------------------------------
  
  sdssproc, flatname[0], indir=indir, spectrographid=spectrographid, color=color

  nflat = n_elements(flatname)
    
  ;---------------------------------------------------------------------------
  ; LOOP THROUGH FLATS + TRACE
  ;---------------------------------------------------------------------------
  splog, 'LOOP THROUGH FLATS-----------------------------'  
  
  flatstruct = create_flatstruct(nflat)
  
  for iflat=0, nflat-1 do begin
  
    splog, iflat+1, nflat, format='("Tracing flat #",I3," of",I3)'
    
    ;---------------------------------------------------------------------
    ; Read flat-field image
    ;---------------------------------------------------------------------

if ~keyword_set(SKIPPROC) then begin     
    splog, 'Reading flat ', flatname[iflat]
    sdssproc, flatname[iflat], flatimg, flativar, indir=indir, hdr=flathdr, /applybias, /applypixflat, nsatrow=nsatrow, fbadpix=fbadpix,$
      ecalibfile=ecalibfile, minflat=minflat, maxflat=maxflat,/applycrosstalk
      
    configuration=obj_new('configuration', sxpar(flathdr, 'MJD'))

    ; Decide if this flat is bad
    qbadflat = reject_flat(flatimg, flathdr, nsatrow=nsatrow, fbadpix=fbadpix, percent80thresh=configuration->spcalib_reject_calib_percent80thresh())
endif else begin
    flatimg = mrdfits(indir+'/'+flatname[iflat],0,flathdr)

    tmp = getenv('BOSS_SPECTRO_DATA')+strtrim(sxpar(flathdr,'MJD'),2)+'/sdR-'+strtrim(sxpar(flathdr,'CAMERAS'),2)+'-00'+strtrim(sxpar(flathdr,'EXPOSURE'),2)+'.fit.gz'
    sdssproc, tmp, dum, flativar
    qbadflat=0
    configuration=obj_new('configuration', sxpar(flathdr, 'MJD'))
endelse

    flatinfofile = string(format='(a,i8.8,a)',flatinfoname, sxpar(flathdr, 'EXPOSURE'), '.fits')

    ;if set then write out the 2d bias-subtraced, pixel-flatted flat field
    if keyword_set(writefiles) then begin
      rawflatname=str_replace(flatinfofile,'spFlat','spFlat-2d')
      mwrfits, flatimg, rawflatname, /create
      splog, 'Writing 2d flat image to file: ', rawflatname 
    endif
 
    if (~keyword_set(fibermask)) then tmp_fibmask = 0 else tmp_fibmask = fibermask

    if (~qbadflat) then begin
      ;------------------------------------------------------------------
      ; Create spatial tracing from flat-field image
      ;------------------------------------------------------------------
      splog, 'Tracing fibers in ', flatname[iflat]
      xsol = trace320crude(flatimg, flativar, yset=ycen, maxdev=1.0, fibermask=tmp_fibmask, cartid=cartid, xerr=xerr, flathdr=flathdr, $
       padding=configuration->spcalib_trace320crude_padding(), plottitle=plottitle+' Traces '+flatname[iflat],visual=visual,fiberparam=fiberparam)

      ;retrieve bundle and fiber info
      nbundle = fiberparam.nbundle
      bundleid = fiberparam.bundleid
      nfiber=fiberparam.nfiber
      radius = fiberparam.radius
        
      splog, 'Fitting traces in ', flatname[iflat]
      ntrace = (size(xsol, /dimens))[1]
      outmask = 0
      ; Ignore values whose central point falls on a bad pixel
      ; ASB: New recipe for inmask, just masking fully useless rows, since trace320crude has already done clever fill-ins:
      inmask = (total(flativar gt 0., 1) gt 0.) # replicate(1B, ntrace)
      xy2traceset, ycen, xsol, tset, ncoeff=configuration->spcalib_xy2traceset_ncoeff(), maxdev=0.5, outmask=outmask, /double, xerr=xerr, inmask=inmask

      junk = where(outmask EQ 0, totalreject)
      if (totalreject GT configuration->spcalib_rejecttheshold()) then begin
        splog, 'Reject flat ' + flatname[iflat] + ': ' + string(format='(i8)', totalreject) + ' rejected pixels'
        qbadflat = 1
      endif
      
      traceset2xy, tset, ycen, xsol
      flatstruct[iflat].tset = ptr_new(tset)
      flatstruct[iflat].xsol = ptr_new(xsol)
      flatstruct[iflat].fibermask = ptr_new(tmp_fibmask)
    endif else begin
      ;bad flat
      xsol = 0
      flatstruct[iflat].qbad = 1
    endelse
    
    ;----------
    ; Verify that traces are separated by > 3 pixels   ;modify this for 2 pixel radius for sims -> to varying radius
    if (qbadflat EQ 0) then begin
      print, 'CHECKING IF FIBERS SEPARATED BY > 3 pixels - PRODUCES WARNING if TOO CLOSE'
      sep = xsol[*,1:ntrace-1] - xsol[*,0:ntrace-2]
      tooclose = where(sep LT 3.0)
      if (tooclose[0] NE -1) then begin
        splog, 'Reject flat ' + flatname[iflat] + ': Traces not separated by more than '+strtrim(3.0,2)+' pixels'
      endif
    endif
    
    print, 'QBADFLAT: ',qbadflat
    if (~qbadflat) then begin
      ;---------------------------------------------------------------------
      ; Extract the flat-field image to obtain width and flux
      ;---------------------------------------------------------------------

      print, 'EXTRACTING FLAT FIELD----------------------'
      sigma = configuration->spcalib_sigmaguess() ; Initial guess for gaussian width
      highrej = 15
      lowrej = 15
      npoly = 10 ; Fit 1 terms to background
      wfixed = [1,1] ; Fit the first gaussian term + gaussian width

      ;BOSS better with proftype=1
      proftype = configuration->spcalib_extract_image_proftype()  
      splog, 'Extracting flat with proftype=', proftype

      extract_image, flatimg, flativar, xsol, sigma, flux, fluxivar, proftype=proftype, wfixed=wfixed, highrej=highrej, lowrej=lowrej, $
                      npoly=npoly, relative=1, ansimage=ansimage, reject=[0.1, 0.6, 0.6], chisq=chisq3, visual=visual, surve=survey,ymodel=ymodel

      widthset3 = fitflatwidth(flux, fluxivar, ansimage, tmp_fibmask, ncoeff=configuration->spcalib_fitflatwidth_ncoeff(), sigma=sigma, $
                            medwidth=medwidth, mask=configuration->spcalib_fitflatwidth_mask(flux,fluxivar), $
                            inmask=configuration->spcalib_fitflatwidth_inmask(flux,fluxivar,ntrace), /double, nbundle=nbundle,nfiber=nfiber)
      widthset = widthset3
      
      ;relieve some memory
      ansimage = 0
      widthset3 = 0

      junk = where(flux GT 1.0e5, nbright)  ;select out bright pixels     
      splog, 'Using proftype=', proftype
      splog, 'Found ', nbright, ' bright pixels in extracted flat ', flatname[iflat], format='(a,i7,a,a)'
        
      flatstruct[iflat].proftype  = proftype
      flatstruct[iflat].fibermask = ptr_new(tmp_fibmask)
      flatstruct[iflat].widthset = ptr_new(widthset)
      flatstruct[iflat].medwidth  = medwidth
      
    endif
    
    flatstruct[iflat].name = flatname[iflat]
    get_tai, flathdr, tai_beg, tai_mid, tai_end
    flatstruct[iflat].tai = tai_mid
    flatstruct[iflat].qbad = qbadflat
    obj_destroy,configuration
    
    if flatstruct[iflat].qbad eq 1 then splog, 'WARNING: QBAD is 1 for FLAT, FLAT REJECT' else splog, 'GOOD FLAT!, '+flatname[iflat]    
    
    ;write out the flat field extraction after the 1st pass -- extraction and determination of profile widths
    if keyword_set(writefiles) then begin
      flat1extname = str_replace(flatinfofile,'spFlat','spFlat-1flux')
      mwrfits, flux, flat1extname, /create
      flat1modname=str_replace(flatinfofile,'spFlat','spFlat-1model')
      mwrfits, ymodel, flat1modname, /create  
      spawn, ['gzip', '-f', flat1extname], /noshell
      spawn, ['gzip', '-f', flat1modname], /noshell     
    endif

  endfor
  
  ;--------REDUCE THE ARCS
  if ~keyword_set(flat) then reduceArcs, arcname, indir=indir, skipproc=skipproc, arcstruct=arcstruct, flatstruct=flatstruct, arcinfoname=arcinfoname


  ;---------------------------------------------------------------------------
  ; LOOP THROUGH FLATS + CREATE FIBERFLATS
  ;---------------------------------------------------------------------------
  splog, 'CREATING FIBERFLATS----------------------'
  for iflat=0, nflat-1 do begin
  
    splog, iflat+1, nflat, format='("Create fiberflats for flat #",I3," of",I3)'
      
    ;----------
    ; Identify the nearest arc for each flat-field, which must be within TIMESEP seconds and be good.
    iarc = -1
    igood = where(arcstruct.qbad EQ 0)
    if (igood[0] NE -1) then begin
      tsep = min( abs(flatstruct[iflat].tai - arcstruct[igood].tai), ii )
      if (tsep LE timesep AND timesep NE 0) then iarc = igood[ii]
      flatstruct[iflat].tsep = tsep
    endif
    
    if (iarc GE 0) then splog, 'Flat ' + flatname[iflat] + ' paired with arc ' + arcname[iarc] $
    else begin
      splog, 'Flat ' + flatname[iflat] + ' paired with no arc'
      flatstruct[iflat].qbad = 1 ; Flat is bad if no companion arc exists
    endelse
    
    flatstruct[iflat].iarc = iarc
    
    if (~flatstruct[iflat].qbad) then begin
    
      widthset = *(flatstruct[iflat].widthset)
      wset = *(arcstruct[iarc].wset)
      xsol = *(flatstruct[iflat].xsol)
      tmp_fibmask = *(flatstruct[iflat].fibermask)
      proftype = flatstruct[iflat].proftype
      
      ;---------------------------------------------------------------------
      ; Read flat-field image (again)
      ;---------------------------------------------------------------------
      
      ; If there is only 1 flat image, then it's still in memory
      if (nflat GT 1) then begin
        splog, 'Reading flat ', flatname[iflat]
        sdssproc, flatname[iflat], flatimg, flativar, indir=indir, hdr=flathdr, /applybias, /applypixflat, ecalibfile=ecalibfile, minflat=minflat, maxflat=maxflat,/applycrosstalk
      endif
      configuration=obj_new('configuration',sxpar(flathdr, 'MJD'))

      ;---------------------------------------------------------------------
      ; Extract the flat-field image
      ;---------------------------------------------------------------------
      
      traceset2xy, widthset, xx, sigma2   ; sigma2 is real width
      highrej = 15
      lowrej = 15
      npoly = 5 ; Fit 5 terms to background, just get best model
      wfixed = [1,0] ; Do not refit for Gaussian widths, only flux ???

      ;TEMP MANGA SOLUTION - use extract_bundle_image for red, extract_image for blue
      if (color EQ 'red') then begin
            extract_bundle_image, flatimg, flativar, xsol, sigma2, flux, fluxivar, proftype=proftype, wfixed=wfixed, highrej=highrej, lowrej=lowrej, $
              npoly=2L, relative=1, chisq=schisq, ansimage=ansimage2, reject=[0.1, 0.6, 0.6], ymodel=ymodel, nperbun=20L, buffsize=8L, $
              visual=visual, survey=survey,nbundle=nbundle, nfiber=nfiber, radius=radius, bundleid=bundleid
      endif
      
      if (color EQ 'blue') then begin
            extract_image, flatimg, flativar, xsol, sigma2, flux, fluxivar, proftype=proftype, wfixed=wfixed, highrej=highrej, lowrej=lowrej, $
             npoly=npoly, relative=1, ansimage=ansimage, reject=[0.1, 0.6, 0.6], chisq=chisq3, visual=visual, survey=survey,ymodel=ymodel
      endif

      ;write out the flat field extraction after the 2nd pass -- actual flux extraction over blocks or rows 
      if keyword_set(writefiles) then begin
        flat2extname=str_replace(flatinfofile,'spFlat','spFlat-2flux')
        mwrfits, flux, flat2extname, /create
        flat2modname=str_replace(flatinfofile,'spFlat','spFlat-2model')
        mwrfits, ymodel, flat2modname, /create  
        spawn, ['gzip', '-f', flat2extname], /noshell
        spawn, ['gzip', '-f', flat2modname], /noshell 
      endif

      if (keyword_set(bbspec)) then begin
         basisfile = 'spBasisPSF-*-'+strmid(arcstruct[iarc].name,4,11)+'.fits'
         tmproot = 'tmp-'+strmid(flatstruct[iflat].name,4,11)
         bbspec_extract, flatimg, flativar, bbflux, bbfluxivar, basisfile=basisfile, ximg=xsol, ymodel=bb_ymodel, tmproot=tmproot, /batch 

         ; Deal with case of only the first few spectra being re-extracted...
         dims = size(bbflux,/dimens)
         flux[0:dims[0]-1,0:dims[1]-1] = bbflux
         fluxivar[0:dims[0]-1,0:dims[1]-1] = bbfluxivar * (fluxivar[0:dims[0]-1,0:dims[1]-1] GT 0) ; <- Retain old rejection

         outfile = 'ymodel-'+strmid(flatstruct[iflat].name,4,11)+'.fits'
         mwrfits, bb_ymodel, outfile, /create
         mwrfits, ymodel, outfile
      endif
        
      ;---------------------------------------------------------------------
      ; Compute fiber-to-fiber flat-field variations
      ;---------------------------------------------------------------------
      sigma2 = 0
      xsol = 0
      
      fflat = fiberflat(flux, fluxivar, wset, fibermask=tmp_fibmask, /dospline, pixspace=5, plottitle=plottitle+' Superflat '+flatstruct[iflat].name, $
        badflatfracthresh=configuration->spcalib_fiberflat_badflatfracthresh(), minval=configuration->spcalib_fiberflat_minval(flux), visual=visual, $
        slitmap=slitmap, fiberparam=fiberparam, origff=origff, survey=survey, superflatset=superflatset, medval=medval)
        
      if (n_elements(fflat) EQ 1) then begin
        flatstruct[iflat].qbad  = 1
        splog, 'Reject flat ' + flatname[iflat] + ': No good traces'
      endif
      
      flatstruct[iflat].fflat = ptr_new(fflat)
      flatstruct[iflat].superflatset = ptr_new(superflatset)
      flatstruct[iflat].fibermask = ptr_new(tmp_fibmask)
      
      ;------------------------------------------------------------------
      ; Write information on flat field processing
      if (keyword_set(flatinfoname)) then begin
        sxaddpar, flathdr, 'NBRIGHT', nbright, 'Number of bright pixels (>10^5) in extracted flat-field'
        flatinfofile = string(format='(a,i8.8,a)',flatinfoname, sxpar(flathdr, 'EXPOSURE'), '.fits')
        print, 'WRITING FLAT: ',flatinfofile

        mwrfits, *flatstruct[iflat].fflat, flatinfofile, flathdr, /create
        mwrfits, *flatstruct[iflat].tset, flatinfofile
        mwrfits, *flatstruct[iflat].fibermask, flatinfofile
        mwrfits, *flatstruct[iflat].widthset, flatinfofile
        mwrfits, *flatstruct[iflat].superflatset, flatinfofile

        ;Add extensions for where everything is on the same lambda grid.
        getCommonMangaWave, wset, fflat, FILE=flatinfofile
        
        spawn, ['gzip', '-f', flatinfofile], /noshell
        
        ;write the flat image model if requested
        if keyword_set(writeflatmodel) then begin
           flatmodelfile = string(format='(a,i8.8,a)',flatinfoname + 'MODELIMG-', sxpar(flathdr, 'EXPOSURE'), '.fits')
           mwrfits, flatimg, flatmodelfile, /create
           mwrfits, flativar, flatmodelfile
           mwrfits, ymodel, flatmodelfile
           spawn, ['gzip', '-f', flatmodelfile], /noshell
        endif
        ymodel = 0
      endif
      
      obj_destroy,configuration
    endif
  endfor

  ;write out calibration structures for the case where you want to skip the reduction of them next time
  save, flatstruct, arcstruct, filename='calibstruct.sav', description='structures containing the flat and arc parameters'
  
  splog, 'Elapsed time = ', systime(1)-stime1, ' seconds', format='(a,f6.0,a)'
  return
end
;------------------------------------------------------------------------------
