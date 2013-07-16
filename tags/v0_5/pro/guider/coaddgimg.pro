;+
;Propose:
;  Stack all guider frames during a science exposure so that we have a time-integrated PSF.
;Usage:
;  pro coaddgimg, mjd, exposure
;Input:
;  mjd  --- mjd
;  exposure --- exposure number, given as an integer
;Output:
;  none
;  The routine will write a coadd image from the guider.
;  the image will be read in by another routine to produce a parametered output
;Revision History
;   Written by Renbin Yan on Jan 10, 2013
;-
pro coaddgimg,mjd,exposure
common subimageblock, stamp,errstamp,nx
common flatinfo, flatfile_inmem, flat,hdr_flat
common darkinfo, darkfile_inmem, dark,hdr_dark

  mjdstr = string(mjd,format='(i0.0)')
  guider_dir = getenv('GCAM_DATA')+'/'+mjdstr+'/'
  outputdir = getenv('GCAM_REDUX')+'/images/'

  expodir = getenv('BOSS_SPECTRO_DATA')+'/'+mjdstr+'/'
  expostr = string(exposure,format='(i8.8)')
  files = 'sdR-??-'+expostr+'.fit*'
  expofiles=file_search(expodir,files,count=ct)
  if ct eq 0 then message,'No files found for exposure:'+expostr+' in path:'+expodir
  
  file = expofiles[0]
  sci_hdr =headfits(file)
  guider1=sxpar(sci_hdr,'GUIDER1')
  guidern=sxpar(sci_hdr,'GUIDERN')
  seqstart = fix(strmid(guider1,10,4))
  seqlast = fix(strmid(guidern,10,4))
  nimg = seqlast-seqstart+1

  filestr = mjdstr+'-'+expostr


  for i=0,nimg-1 do begin
    seqno = string(seqstart+i,format='(i4.4)')

    origin = 'gimg-'+seqno+'.fits'
    proc ='proc-gimg-'+seqno+'.fits'
    exist1 = file_search(guider_dir+origin,count=c1)
    exist2 = file_search(guider_dir+proc,count=c2)
    if c1 ne 1 or c2 ne 1 then begin
       error = 1
       print,'Missing file:,',guider_dir+origin, ' or ',guider_dir+proc
       print,'Coadding stopped and program quit, no file written.'
       return
    endif 
    
    proc_hdr=headfits(guider_dir+proc)
    if i eq 0 then begin 
       dra=sxpar(proc_hdr,'DCNRA')
       ddec=sxpar(proc_hdr,'DCNDEC')
    endif else begin
       if dra ne sxpar(proc_hdr,'DCNRA') or ddec ne sxpar(proc_hdr,'DCNDEC') then message,'DCNRA or DCNDEC changed during this science exposure. Stop.'
    endelse

    img=mrdfits(guider_dir+origin,0,/unsigned,hdr,/silent)
    imgsz = size(img,/dimen)
    if imgsz[0] ne 524 or imgsz[1] ne 512 then begin
      message,'Error: I am expecting the guider frames to be 524x512. But it is not! Does that change?'
    endif

    flatfile0=sxpar(hdr,'FLATFILE')
    if n_elements(flatfile) eq 0 then begin
      flatfile=flatfile0
    endif else begin
      if flatfile0 ne flatfile then message,'Error: Guider images to be combined are using more than one flat frame. Are you sure the guider frames are all taken during one exposure? '
    endelse

    darkfile0 = sxpar(hdr,'DARKFILE')
    if n_elements(darkfile_inmem) eq 0 then begin
       p0=strpos(darkfile0,'gimg')
       darkfile = strmid(darkfile0,p0,strlen(darkfile0)-p0)
       dark=mrdfits(guider_dir+darkfile,0,/unsigned,hdr_dark)
       darkfile_inmem=darkfile0
    endif else if darkfile0 ne darkfile_inmem then begin
       p0 = strpos(darkfile0,'gimg')
       darkfile = strmid(darkfile0,p0,strlen(darkfile0)-p0)
       dark=mrdfits(guider_dir+darkfile,0,/unsigned,hdr_dark)
       darkfile_inmem=darkfile0
    endif

    ;The guide camera has very little dark current. This is shown by comparing a 4 sec image with a 20 sec dark exposure. They are almost identical with a constant offset (the longer exposure image is lower in counts). But the bias has a shape: it gets stronger at small x. Therefore, we should use an offset dark image to bias-subtract the science image. 

    ;bias subtract image
    darksize=size(dark,/dimen)
    darkx = darksize[0] & darky=darksize[1]
    if darkx eq imgsz[0] and darky eq imgsz[1] then begin
      biasdark=median(dark[darky:darkx-6*darky/512,*])
      biasimg =median(img[darky:darkx-6*darky/512,*])
      darksmooth=median(dark,dimen=2)#(fltarr(darky)+1)
      bias_to_sub = float(darksmooth)+(float(biasimg)-float(biasdark))
      bimg = float(img)-bias_to_sub
      bimg = bimg-median(bimg)
    endif else begin
      print,'Warning: Dark image size does not match the guider image size --- cannot use the dark image for dark subtraction.'
      print,'Instead we will use simple median to subtract the image, which is less accurate, especially for fiber 7 and 3' 
      bimg = float(img)-median(img)
    endelse

    if i eq 0 then begin 
       totimg = bimg 
       totexp = sxpar(hdr,'EXPTIME')
    endif else begin 
       totimg=totimg+bimg
       totexp = totexp + sxpar(hdr,'EXPTIME')
    endelse
  endfor

  if n_elements(flatfile_inmem) eq 0 then begin
       p0=strpos(flatfile0,'gimg')
       flatfile = strmid(flatfile0,p0,strlen(flatfile0)-p0)
       flat=mrdfits(guider_dir+flatfile,0,/unsigned,hdr_flat)
       flatfile_inmem=flatfile0
  endif else if flatfile0 ne flatfile_inmem then begin
       p0=strpos(flatfile0,'gimg')
       flatfile = strmid(flatfile0,p0,strlen(flatfile0)-p0)
       flat=mrdfits(guider_dir+flatfile,0,/unsigned,hdr_flat)
       flatfile_inmem=flatfile0
  endif


  ;bias subtract flat using the overscan
  if (size(flat,/dimen))[0] ne 1048 then begin
      message,'Error: I am expecting the size for the flat frame to be 1048 x1024. But it is not!'
  endif
  bflat = flat-median(flat[1024:1038,*])
  ;bin flat by 2x2
  flatsize=size(bflat,/dimen)
  binflat = rebin(bflat,flatsize[0]/2,flatsize[1]/2.)
  thresh = percentiles(binflat,96.5)
  bitmask = byte(binflat*0)
  bitmask[where(binflat le thresh)] = 4B
;stop
  binflat = binflat/median(binflat[where(bitmask eq 0B)])

  ;flat image
  fimg = totimg/binflat*(bitmask eq 0B)+totimg*(bitmask eq 4B)
  k=where(finite(fimg) eq 0,ct_inf)
  if ct_inf gt 0 then fimg[k]=0.0

; readnoise 10.7 ADU, 17.1 e-
  errimg = sqrt(abs(fimg)+nimg*17.^2)*(bitmask eq 0B) + 1.e30*(bitmask eq 4B)
  if ct_inf gt 0 then errimg[k] = 1.e30

  mkhdr,header,fimg
  sxdelpar,header,'DATE'
  sxaddpar,header,'TIMESYS',sxpar(sci_hdr,'TIMESYS'),before='COMMENT'
  sxaddpar,header,'DATE-OBS',sxpar(sci_hdr,'DATE-OBS'),before='COMMENT'
  sxaddpar,header,'GUIDER1',guider1,'The first guider image',before='COMMENT'
  sxaddpar,header,'GUIDERN',guidern,'The last guider image',before='COMMENT'
  sxaddpar,header,'TOTEXPTIME',totexp,'Total guider exposure time',before='COMMENT'
  sxaddpar,header,'OBSTIME',sxpar(sci_hdr,'EXPTIME'),'Total science exposure time',before='COMMENT'
  sxaddpar,header,'DARKFILE',darkfile0,'guider dark exposure',before='COMMENT'
  sxaddpar,header,'FLATFILE',flatfile0,'guider flat exposure',before='COMMENT'
  sxaddpar,header,'CARTID',sxpar(sci_hdr,'CARTID'),'The currently loaded cartridge',before='COMMENT'
  sxaddpar,header,'PLATEID',sxpar(sci_hdr,'PLATEID'),'The currently loaded plate',before='COMMENT'
  sxaddpar,header,'EXPOSURE',sxpar(sci_hdr,'EXPOSURE'),before='COMMENT' 
  sxaddpar,header,'MJD',sxpar(sci_hdr,'MJD'),'APO fMJD day at start of exposure',before='COMMENT'
  sxaddpar,header,'DCNRA',dra,'applied user supplied offset in RA, arcsec',before='COMMENT'
  sxaddpar,header,'DCNDEC',ddec,'applied user supplied offset in Dec, arcsec',before='COMMENT'

  mwrfits,fimg,outputdir+'cogimg-'+filestr+'.fits',header,/create
  mwrfits,errimg,outputdir+'cogimg-'+filestr+'.fits'
  mwrfits,bitmask,outputdir+'cogimg-'+filestr+'.fits'
end

