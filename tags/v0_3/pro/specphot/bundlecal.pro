pro bundlecal,exposure

    if n_elements(exposure) gt 1 then begin 
      for i=0,n_elements(exposure)-1 do bundlecal,exposure[i]
      return
    endif

    gcam_redux_dir = getenv('GCAM_REDUX')
    summary = mrdfits(gcam_redux_dir+'gpsfcat.fits',1)

    t = where(summary.exposure eq exposure, have) 
    if have ne 1 then message,'Error: '+string(have, format='(i0.0)')+' observation found. I am expecting just one.'

    observation=summary[t]
    manga_redux_dir = getenv('MANGA_SPECTRO_REDUX')
    mjdstr = string(observation.mjd,format='(i0.0)')
    expstr = string(observation.exposure,format='(i8.8)')
    dir = manga_redux_dir+string(observation.plate,format='(i0.0)')+'/'+mjdstr+'/'
    bfile = file_search(dir,'spFrame-b1-'+expstr+'.*',count=ct_b)
    rfile = file_search(dir,'spFrame-r1-'+expstr+'.*',count=ct_r)

    if ct_b ne 1 or ct_r ne 1 then message, 'Error: Found '+string(ct_b)+' spFrame file for the blue side and '+string(ct_r)+' spFrame file for the red side. I am expecting 1 and 1.'

;    mlframe_read,bfile,objflux=bflux,objivar=bivar,mask=bmask,loglam=bloglam,dispimg=bdispimg,plugmap=plugmap,wset=bwset
;    mlframe_read,rfile,objflux=rflux,objivar=rivar,mask=rmask,loglam=rloglam,dispimg=rdispimg,plugmap=plugmap,wset=rwset

;    correct_dlam,bflux,bivar,bwset,dlam=dloglam
;    correct_dlam,bdispimg,0,bwset,dlam=dloglam,/inverse
;    correct_dlam,rflux,rivar,rwset,dlam=dloglam
;    correct_dlam,rdispimg,0,rwset,dlam=dloglam,/inverse
   
    bundles = ['ma001','ma002','ma003','ma004','ma005','ma006','ma007','ma008']
    if observation.plate eq 6870 then bundles=bundles[1:7]
    for i=0,n_elements(bundles)-1 do begin
       bundlefluxrationew,exposure,bundles[i],newloglam,mratfit,factor=factor,scale=scale,rotation=rotation,xbest=xbest,ybest=ybest
       if i eq 0 then begin 
          npix = n_elements(newloglam)
          single = {loglam:newloglam,bundleid:' ',factor:0.0,scale:0.0,rotation:0.0,xbest:0.0,ybest:0.0,thrupt:fltarr(npix)}
          cat = replicate(single,n_elements(bundles))
          cat.bundleid = bundles
       endif
       cat[i].factor=factor
       cat[i].scale=scale
       cat[i].rotation=rotation
       cat[i].xbest=xbest
       cat[i].ybest=ybest 
       if i eq 0 then totfit = mratfit else totfit = [[totfit],[mratfit]]
    endfor
    cat.thrupt = totfit

    meanfit = total(totfit,2)/8.
    diff = totfit-meanfit#(fltarr(8)+1)
    rms = sqrt(total(diff^2,2)/8.)
    err = rms/sqrt(8.-1)
    plot,10^newloglam#(fltarr(8)+1),totfit,xran=[3500,10500],xst=1,yran=[0,400],xtitle='!6Wavelength (A)',ytitle='observed / model',ps=1,symsize=0.1,charthick=2,charsize=2
    oplot,10^newloglam,meanfit,color=249,thick=2

    plot,10^newloglam,rms/meanfit,yran=[0,1],xran=[3500,10500],xst=1,xtitle='!6Wavelength',ytitle='Fractional scatter or uncertainty of mean',thick=1.5,charthick=1.5,charsize=2
    oplot,10^newloglam,err/meanfit,color=249,thick=1.5
    legend,['RMS scatter/mean','uncertainty of mean /mean'],linest=[0,0],color=[0,249],box=0,charsize=2,charthick=1.5

    npix = (size(totfit,/dimen))[0]
    mwrfits,cat,'thrupt_'+string(exposure,format='(i8.8)')+'.fits',/create
end
