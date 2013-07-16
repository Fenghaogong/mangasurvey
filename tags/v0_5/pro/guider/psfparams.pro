;+
; Purpose:
;   This routine will fit the guider star images to models and produce a
;   parameterized PSF to be used in spectrophotometry calibration etc.
;-

function psfparams, mjd, exposure
common subimageblock, stamp,errstamp,nx

   gdatadir=getenv('GCAM_DATA')
; Read in the coadded, flat-fielded image corrsponding to the mjd and exposure. 
   gdir=getenv('GCAM_REDUX')
   mjdstr=string(mjd,format='(i0.0)')
   expostr=string(exposure,format='(i8.8)')
   coaddfile = gdir+'/images/'+'cogimg-'+mjdstr+'-'+expostr+'.fits'
   img = mrdfits(coaddfile,0,hdr_g)
   errimg=mrdfits(coaddfile,1)
   bitmask=mrdfits(coaddfile,2)

   guider1=sxpar(hdr_g,'GUIDER1')
   procfile= file_search(gdatadir+'/'+mjdstr+'/',guider1,count=ct)
   if ct eq 1 then info = mrdfits(procfile,6) else message,'No proc file found with name '+procfile

   ; smoothed image for rough fast centroiding.
   smimg=smooth(img,3)

   sz = size(img,/dimen)
   xsize = sz[0] & ysize = sz[1]
   xcoords = findgen(xsize)#(fltarr(ysize)+1)
   ycoords = (fltarr(xsize)+1)#findgen(ysize)
   
;  smooth the bitmask and take out the pixels on the edge.
   smask= smooth(float(bitmask),[3,3]) ne 0
; Find Acquisition fibers, determine centroid of star, find sky region and measure sky background.

   acqu = where(strmatch(info.fiber_type,'ACQUIRE*') and info.exists eq 'T' and finite(info.xcenter),n_acqu)
   for i=0,n_acqu-1 do begin
       p=where(abs(xcoords-info[acqu[i]].xcenter) lt info[acqu[i]].radius $
           and abs(ycoords-info[acqu[i]].ycenter) lt info[acqu[i]].radius)
       tmp = max(smimg[p],maxind)
       xc0 = xcoords[p[maxind]]
       yc0 = ycoords[p[maxind]]
       skyarea = where((xcoords[p]-xc0)^2+(ycoords[p]-yc0)^2 gt 11.^2 and smask[p] eq 0B)
       if i eq 0 then totskypix=p[skyarea] else totskypix=[totskypix,p[skyarea]]
;       stop
   endfor
   sky = median(img[totskypix])

;  sky-subtract the image
   sbimg = img-sky*(smask eq 0B)

; Change the bitmask of the region around the Tritium spot so that we don't zero out that region.
    p = where(abs(xcoords-132) lt 30 and abs(ycoords-406) lt 30)
    smask[p] = 8B

   
; Zero out all unmasked regions.
   bad = where(smask ne 0B and smask ne 8B)
   sbimg[bad] = 0.0

    
    info[acqu].radius=8.5
    good=where(info.exists eq 'T' and strmatch(info.fiber_type,'TRITIUM*') eq 0 and finite(info.xcenter),nstar)
    tritium = where(strmatch(info.fiber_type,'TRITIUM*'))
    neff   = fltarr(nstar)
    sigma  = fltarr(nstar)
    theta  = fltarr(nstar)
    flux   = fltarr(nstar)
    gflux   = fltarr(nstar)
    fitxy  =  fltarr(2,nstar)
    gfitarr = fltarr(6,nstar)
    fwhm   = fltarr(nstar)
    chi2   = fltarr(nstar)
    rmsarr = fltarr(nstar)
  
   for i=0,nstar-1 do begin
       p =where(abs(xcoords-info[good[i]].xcenter) lt info[good[i]].radius $
              and abs(ycoords-info[good[i]].ycenter) lt info[good[i]].radius $
              and ((bitmask and 5B) eq 0))
       neff[i] = djs_neff(sbimg[p])
       sigma[i] = sqrt(neff[i]/4/!pi)
       tmp = max(smimg[p],maxind)

       xbest = xcoords[p[maxind]]
       ybest = ycoords[p[maxind]]

       xran = where(abs(xcoords[*,0]-xbest) lt info[good[i]].radius,nx)
       yran = where(abs(ycoords[0,*]-ybest) lt info[good[i]].radius,ny)

       stamp = sbimg[min(xran):max(xran),min(yran):max(yran)]
       errstamp = errimg[min(xran):max(xran),min(yran):max(yran)]
       sub_xc = xbest-min(xran)
       sub_yc = ybest-min(yran)

;       print,'fiberid:',info[good[i]].fiberid,' sigma:',sigma[i] 
       valid = where(errstamp lt 1.d29,n_valid)
       if keyword_set(single) then begin
          sig_guess = sigma[i] > 1.0
          guess = [sub_xc,sub_yc,sig_guess^2]
          gfitres = mpfitfun('gaussprof',findgen(nx*ny),reform(stamp,nx*ny),reform(errstamp,nx*ny),guess,yfit=fit2,covar=cov,maxiter=100,bestnorm=bestnorm,perror=error,status=status,/quiet)
          yfit = gaussprof(findgen(nx*ny),gfitres,amp=amp)
          fwhm[i] = 2.35*sqrt(gfitres[2])
          gflux[i] = 2*!pi*amp[0]*gfitres[2]
          gfitres=[gfitres,0.0]
          chi2[i] = bestnorm/(n_valid-3)
       endif else begin
          sig_guess = sigma[i] > 1.0
          guess = [sub_xc,sub_yc,(sig_guess*0.8)^2,(sig_guess*1.2)^2]
;          yfit = doublegau(findgen(nx*ny),guess) 

          gfitres = mpfitfun('doublegau',findgen(nx*ny),reform(stamp,nx*ny),reform(errstamp,nx*ny),guess,yfit=fit2,covar=cov,maxiter=100,bestnorm=bestnorm,perror=error,status=status,/quiet)
          yfit = doublegau(findgen(nx*ny),gfitres,amps=amps)
;          print,info[good[i]].fiberid," parameter: ",amps,gfitres[2:3]
          sigs = sqrt(gfitres[2:3])
          sigs = sigs(sort(sigs))
          tx = findgen(100)/100*(sigs[1]-sigs[0])*1.175+sigs[0]*1.175
          ty = amps[0]*exp(-tx*tx/(2*gfitres[2]))+amps[1]*exp(-tx*tx/(2*gfitres[3]))
          fwhm[i] = 2*interpol(tx,ty,total(amps)/2.)
          gflux[i] = 2*!pi*(amps##transpose(gfitres[2:3]))
          chi2[i] = bestnorm/(n_valid-4)
       endelse
       subxcoords= findgen(nx)#(fltarr(ny)+1)
       subycoords= (fltarr(nx)+1)#findgen(ny)
       rsq=(subxcoords-gfitres[0])^2+(subycoords-gfitres[1])^2
       rmsarr[i] = sqrt(total(rsq*stamp)/total(stamp))
       gfitarr[*,i] = [gfitres ,amps]
       fit_xc = gfitres[0]+min(xran)
       fit_yc = gfitres[1]+min(yran)
       fitxy[*,i]=[fit_xc,fit_yc]

       ind = where((xcoords[p]-fit_xc)^2 + (ycoords[p]-fit_yc)^2 lt 6.^2,n_pixel)
       if n_pixel gt 0 then flux[i] = total(sbimg[p[ind]]) else flux[i]=0.0
;       stop
    endfor

; Find Tritium spot and measure its size
    p = where(abs(xcoords-132) lt 30 and abs(ycoords-406) lt 30)
    tmp = max(smimg[p],maxind)
    xbest = xcoords[p[maxind]]
    ybest = ycoords[p[maxind]]
    tri_stamp = sbimg[xbest-5:xbest+5,ybest-5:ybest+5]
    tri_neff = djs_neff(tri_stamp)

    select=addtag(info[good],{mjd:0L,plate:0L,FLATCART:0,neff:0.0,rms:0.0,sigma:0.0,gfit:fltarr(6),mfwhm:0.0,fitxy:fltarr(2),measureflux:0.0,measuregflux:0.0,chi2:0.0,starttime:0.0d,dateobs:' ',exptime:0.0,sexptime:0.0,guider1:' ',guidern:' ',dcnra:0.0,dcndec:0.0})
    select.sigma = sigma*0.428
    select.neff = neff
    select.rms = rmsarr*0.428
    select.gfit = gfitarr
    select.mfwhm = fwhm*0.428
    select.fitxy = fitxy 
    select.measureflux = flux
    select.measuregflux = gflux
    select.chi2  = chi2

    tri_str=addtag(info[tritium],{mjd:0L,plate:0L,FLATCART:0,neff:0.0,rms:0.0,sigma:0.0,gfit:fltarr(6),mfwhm:0.0,fitxy:fltarr(2),measureflux:0.0,measuregflux:0.0,chi2:0.0,starttime:0.0d,dateobs:' ',exptime:0.0,sexptime:0.0,guider1:' ',guidern:' ',dcnra:0.0,dcndec:0.0})
    tri_str.neff = tri_neff
    tri_str.fitxy = [xbest,ybest]
    select = [select,tri_str]

    select.flatcart = sxpar(hdr_g,'CARTID')
    select.plate = sxpar(hdr_g,'PLATEID')
    select.mjd = mjd
    select.exptime = sxpar(hdr_g,'TOTEXPTI')
    select.sexptime = sxpar(hdr_g,'OBSTIME')
    select.guider1 = sxpar(hdr_g,'GUIDER1')
    select.guidern = sxpar(hdr_g,'GUIDERN')

    if mjd le 56288 then unitscale=1/0.0604822 else unitscale=1.0
    select.dcnra = sxpar(hdr_g,'DCNRA')*unitscale
    select.dcndec = sxpar(hdr_g,'DCNDEC')*unitscale

    dateobs=sxpar(hdr_g,'DATE-OBS')
    select.dateobs = dateobs
    timestr = (strsplit(dateobs,'T',/extract))[1]
    parts = strsplit(timestr,':',/extract)
    hr = fix(parts[0]) & min = fix(parts[1])  & sec = double(strmid(parts[2],0,4))
    time = (hr+min/60.+sec/3600.)/24.0
    select.starttime = time
    error = 0
;    print,'Time',systime(/sec)-t0
    return,select
end
