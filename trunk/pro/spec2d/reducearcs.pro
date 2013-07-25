



pro reduceArcs, arcname, indir=indir, SKIPPROC=SKIPPROC, arcstruct=arcstruct, flatstruct=flatstruct, arcinfoname=arcinfoname

   on_error, 0
   compile_opt idl2
   
   narc = n_elements(arcname)
  
    ;---------------------------------------------------------------------------
    ; LOOP THROUGH ARCS + FIND WAVELENGTH SOLUTIONS
    ;---------------------------------------------------------------------------
    splog, 'LOOP THROUGH ARCS-------------------------'
    
    arcstruct = create_arcstruct(narc)
    
    for iarc=0, narc-1 do begin
    
      splog, iarc+1, narc, format='("Extracting arc #",I3," of",I3)'
      
      ;---------------------------------------------------------------------
      ; Read the arc
      ;---------------------------------------------------------------------
  
      if ~keyword_set(SKIPPROC) then begin    
          splog, 'Reading arc ', arcname[iarc]
          
          sdssproc, arcname[iarc], arcimg, arcivar, indir=indir, hdr=archdr, /applybias, /applypixflat, nsatrow=nsatrow, fbadpix=fbadpix, $
            ecalibfile=ecalibfile, minflat=minflat, maxflat=maxflat,/applycrosstalk
          ny = (size(arcimg,/dimens))[1]
            
          configuration=obj_new('configuration', sxpar(archdr, 'MJD'))
          splog, 'Fraction of bad pixels in arc = ', fbadpix
          
          ;----------
          ; Decide if this arc is bad
          qbadarc = reject_arc(arcimg, archdr, nsatrow=nsatrow, fbadpix=fbadpix)
          
      endif else begin
          ;sdssproc, arcname[iarc], arcimg, arcivar, indir=indir, hdr=archdr
          arcimg = mrdfits(indir+'/'+arcname[iarc],0,archdr)
          tmp = getenv('BOSS_SPECTRO_DATA')+strtrim(sxpar(archdr,'MJD'),2)+'/sdR-'+strtrim(sxpar(archdr,'CAMERAS'),2)+'-00'+strtrim(sxpar(archdr,'EXPOSURE'),2)+'.fit.gz'
          sdssproc, tmp, dum, arcivar 
          qbadarc=0  
          ny = (size(arcimg,/dimens))[1]
          fbadpix=0
          configuration=obj_new('configuration', sxpar(archdr, 'MJD')) 
      endelse

     ;----------
     ; Identify the nearest flat-field for this arc, which must be within TIMESEP seconds and be a good flat.      
      get_tai, archdr, tai_beg, tai_mid, tai_end
      tai = tai_mid
      
      iflat = -1
      igood = where(flatstruct.qbad EQ 0)
      if (igood[0] NE -1) then begin
        tsep = min( abs(tai - flatstruct[igood].tai), ii )
        if (tsep LE timesep AND timesep NE 0) then iflat = igood[ii]
      endif
      
      if (iflat GE 0) then begin
        print, 'FOUND FLAT TO PAIR WITH ARC!'
        splog, 'Arc ' + arcname[iarc] + ' paired with flat ' + flatname[iflat]
      endif else begin
        print, 'NO FLAT TO PAIR WITH ARC, SETTING BAD ARC to 1'
        splog, 'Arc ' + arcname[iarc] + ' paired with no flat'
        qbadarc = 1
      endelse
      
      if (~qbadarc) then begin
        xsol = *(flatstruct[iflat].xsol)
        widthset = *(flatstruct[iflat].widthset)
        tmp_fibmask = *(flatstruct[iflat].fibermask)
        proftype = flatstruct[iflat].proftype
        
        ;----------
        ; Calculate possible shift between arc and flat
        xcor = match_trace(arcimg, arcivar, xsol, radius=radius, nfiber=nfiber)
        
        bestlag = median(xcor-xsol)
        if (abs(bestlag) GT 2.0) then begin
          qbadarc = 1
          splog, 'Reject arc: pixel shift is larger than 2 pixel'
          splog, 'Reject arc ' + arcname[iarc] + ': Pixel shift = ', bestlag
        endif
      endif
      
      if (~qbadarc) then begin
        splog, 'Shifting traces with match_trace', bestlag
     
        ;---------------------------------------------------------------------
        ; Extract the arc image
        ;---------------------------------------------------------------------

        traceset2xy, widthset, xx, traced_sigma
        
        highrej = 15
        lowrej = 15
        wfixed = [1,0] ; ASB: Don't fit for width terms.
  
        splog, 'Extracting arc'
        extract_bundle_image, arcimg, arcivar, xcor, traced_sigma, flux, fluxivar, proftype=proftype, wfixed=wfixed, highrej=highrej, lowrej=lowrej, npoly=2L, relative=1, $
          reject=[0.1, 0.6, 0.6], ymodel=ymodel, nperbun=20L, buffsize=8L, visual=visual, survey=survey, nbundle=nbundle, nfiber=nfiber, radius=radius, bundleid=bundleid
          
        ;flag to determine whether or not to do 2-phase arc solution:
        twophase = sxpar(archdr, 'TWOPHASE')
        if keyword_set(twophase) then splog, 'Setting 2-phase readout flag'
  
        if keyword_set(VISUAL) then begin
          getwindow,/open
          h=bytscl(flux,min=0,max=255)
          ml_tvimage, h, /axis, axkeywords={charsize:2,xtitle:'Row',ytitle:'Fiber Number',title:'Extracted Arc Image - '+string(f='(a1,i1)',color, spectrographid)}
        endif
  
        ;---------------------------------------------------------------------
        ; Compute correlation coefficient for this arc image
        ;---------------------------------------------------------------------
  
        splog, 'Searching for wavelength solution'
        aset = 0

        fitarcimage, flux, fluxivar, aset=aset, color=color, lampfile=lampfile, fibermask=tmp_fibmask, bestcorr=bestcorr, $
                    acoeff=configuration->spcalib_arcfitguess_acoeff(color), dcoeff=configuration->spcalib_arcfitguess_dcoeff(color), $
                    wrange=configuration->spcalib_fitarcimage_wrange(color), twophase=twophase, visual=visual, survey=survey, $
                    nbundle=nbundle, nfiber=nfiber, radius=radius
  
        arcstruct[iarc].bestcorr = bestcorr
        
        if ((color EQ 'blue' AND bestcorr LT 0.5) OR (color EQ 'red'  AND bestcorr LT 0.5) ) then begin
          qbadarc = 1
          splog, 'Reject arc ' + arcname[iarc] + ': correlation is only = ' + string(format='(i4)', bestcorr)
        endif
      endif
      
      ;if we don't have a bad arc then continue
      if (~qbadarc) then begin
      
        ;---------------------------------------------------------------------
        ; Compute wavelength calibration
        ;---------------------------------------------------------------------
        arccoeff = configuration->spcalib_arccoeff()
        splog, 'Searching for wavelength solution'
        fitarcimage, flux, fluxivar, xpeak, ypeak, wset, ncoeff=arccoeff, aset=aset, color=color, lampfile=lampfile, fibermask=tmp_fibmask, lambda=lambda, $
                      rejline=rejline, xdif_tset=xdif_tset, acoeff=configuration->spcalib_arcfitguess_acoeff(color), dcoeff=configuration->spcalib_arcfitguess_dcoeff(color), $
                      wrange=configuration->spcalib_fitarcimage_wrange(color), twophase=twophase,visual=visual, survey=survey, nbundle=nbundle, nfiber=nfiber, radius=radius
  
        print, 'Number of Rejected lamp lines: ', n_elements(where(rejline eq 'Reject-offset')), ' out of 45 lamp lines'
  
        ;check for no wset
        if (~keyword_set(wset)) then begin
          print, 'WARNING: NO WSET, SETTING QBADARC TO 1'
          splog, 'Wavelength solution failed'
          qbadarc = 1
        endif else begin
  
          nfitcoeff = configuration->spcalib_ncoeff(color)
          ilamp = where(rejline EQ '')
          dispset = fitdispersion(flux, fluxivar, xpeak[*,ilamp], sigma=configuration->spcalib_sigmaguess(), ncoeff=nfitcoeff, $
                                   xmin=0.0, xmax=ny-1, nfiber=nfiber, medwidth=wsigarr, numbundles=nbundle) 
  
          arcstruct[iarc].dispset = ptr_new(dispset)
          arcstruct[iarc].wset = ptr_new(wset)
          arcstruct[iarc].nmatch = N_elements(lambda)
          arcstruct[iarc].lambda = ptr_new(lambda)
          arcstruct[iarc].rejline = ptr_new(rejline)
          arcstruct[iarc].tsep = tsep
          arcstruct[iarc].xpeak = ptr_new(xpeak)
          arcstruct[iarc].xdif_tset = ptr_new(xdif_tset)
          arcstruct[iarc].fibermask = ptr_new(tmp_fibmask)
          arcstruct[iarc].medwidth = wsigarr
          
          ;------------------------------------------------------------------
          ; Write information on arc lamp processing
          
          if (keyword_set(arcinfoname)) then begin
            sxaddpar, archdr, 'FBADPIX', fbadpix, 'Fraction of bad pixels in raw image'
            sxaddpar, archdr, 'BESTCORR', bestcorr, 'Best Correlation coefficient'
              
            arcinfofile = string(format='(a,i8.8,a)',arcinfoname, sxpar(archdr, 'EXPOSURE'), '.fits')
            
            print, 'WRITING ARC: ',arcinfofile
              
            mwrfits, flux, arcinfofile, archdr, /create
            mwrfits, [transpose(lambda), xpeak], arcinfofile
            mwrfits, *arcstruct[iarc].wset, arcinfofile
            mwrfits, *arcstruct[iarc].fibermask, arcinfofile
            mwrfits, *arcstruct[iarc].dispset, arcinfofile
     
            ;put flux onto common MANGA wavelength grid and add extensions to arc file
            getCommonMangaWave, wset, flux, FILE=arcinfofile
            
            spawn, ['gzip', '-f', arcinfofile], /noshell
            
           ;write arc image model info if requested:
            if keyword_set(writearcmodel) then begin
               arcmodelfile = string(format='(a,i8.8,a)',arcinfoname + 'MODELIMG-', sxpar(archdr, 'EXPOSURE'), '.fits')
               mwrfits, arcimg, arcmodelfile, /create
               mwrfits, arcivar, arcmodelfile
               mwrfits, ymodel, arcmodelfile
               spawn, ['gzip', '-f', arcmodelfile], /noshell
            endif
            ymodel = 0
          endif  ;end of write arc file
        endelse ;end of no wset if statement
      endif  ;end of qbadarc if 
      
      arcstruct[iarc].name = arcname[iarc]
      arcstruct[iarc].tai = tai
      arcstruct[iarc].iflat = iflat
      arcstruct[iarc].qbad = qbadarc
      
      if arcstruct[iarc].qbad eq 1 then splog, 'WARNING: QBAD FOR ARC = 1, ARC REJECT' else splog, 'GOOD ARC!, '+arcname[iarc]
      
      obj_destroy,configuration
    endfor
    
    arcimg = 0
    arcivar = 0
    
  
  
  return

end