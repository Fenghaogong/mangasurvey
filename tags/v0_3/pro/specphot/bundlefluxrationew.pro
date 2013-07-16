function matchbundleratio2, bundleratio, ratioerr, reffiber, distarr, cenwave=cenwave,dwave=dwave,predictratio=predictratio,usemask=usemask,weight=weight,ivarweight=ivarweight,lambda=lambda
;common com_fluxbundle, cube,xarr,yarr,wavearr
common com_fluxbundle2, psfplane,xarr,wavearr
    if n_elements(cenwave) ne n_elements(dwave) then message,'Error: cenwave and dwave must have the same number of elements.'
    if n_elements(weight) ne n_elements(lambda) then message,'Error: lambda and weight must have the same number of elements.'
    if n_elements(weight) ne (size(ivarweight,/dimen))[0] then message,"Error: ivarweight's first dimension must match weight."
    xpos = interpol(findgen(n_elements(xarr)),xarr,distarr)
    nfiber = (size(distarr,/dimen))[0]
    if nfiber ne (size(ivarweight,/dimen))[1] then message,"Error: ivarweight's second dimension must match nfiber."

    nwave = n_elements(cenwave)
    n_wavearr = n_elements(wavearr)
    covfn=fltarr(nwave,nfiber)
    for i=0,nfiber-1 do begin
      covfn_wavearr = interpolate(psfplane,xpos[i,*],indgen(n_wavearr)) 
      covfn_allwave = interpol(covfn_wavearr,wavearr,lambda)
      weightedflux = covfn_allwave*weight*ivarweight[*,i]
      for j=0,nwave-1 do begin
         lam1=cenwave[j]-dwave[j]/2.
	 lam2=cenwave[j]+dwave[j]/2.
	 ind = where(lambda ge lam1 and lambda lt lam2)
         covfn[j,i] = total(weightedflux[ind])/total(ivarweight[ind,i])
      endfor
    endfor
    predictratio=covfn/(covfn[*,reffiber]#(fltarr(nfiber)+1))
    useid = where(indgen(nfiber) ne reffiber and usemask,nuse)
    diff = (bundleratio[*,useid]-predictratio[*,useid])/ratioerr[*,useid]
    chi2 = total(diff*diff)/(nuse*nwave-1)
    return,chi2
end

function mcmcbundleratio,fluxratio,fluxratioerr,offsetx0=offsetx0,offsety0=offsety0,fiberoffx=fiberoffx,fiberoffy=fiberoffy,nfiber=nfiber,n_wavearr=n_wavearr,reffiber=reffiber,centerwave=centerwave,dwave=dwave,innermask=innermask,weight=weight,ivarweight=ivarweight,lambda=lambda,scale=scale,rotation=rotation,xbest=xbest,ybest=ybest,nstep=nstep
      stepsize = [0.06,1,0.04] ; [positional move radius, rotation in degree, scale of ADR]
      if NOT keyword_set(nstep) then nstep =1000L
      random=randomu(seed,nstep*4)
      rr = random[0:nstep-1]
      theta = random[nstep:2*nstep-1]*!pi*2
      xstep = rr*cos(theta)
      ystep = rr*sin(theta)
      anglestep = random[2*nstep:3*nstep-1]-0.5;*stepsize[1]*!dtor
      scalestep = random[3*nstep:4*nstep-1]-0.5;*stepsize[2])
      offsetx = offsetx0
      offsety = offsety0
      dxarr = fiberoffx#(fltarr(n_wavearr)+1)-(fltarr(nfiber)+1)#offsetx
      dyarr = fiberoffy#(fltarr(n_wavearr)+1)-(fltarr(nfiber)+1)#offsety
      distarr = sqrt(dxarr*dxarr+dyarr*dyarr)
      oldchi2 = matchbundleratio2(fluxratio,fluxratioerr,reffiber,distarr,cenwave=centerwave,dwave=dwave,usemask=innermask,weight=weight,ivarweight=ivarweight,lambda=lambda)
      xposarr = fltarr(nstep+1)
      yposarr = fltarr(nstep+1)
      scalearr = fltarr(nstep+1)
      rotationarr = fltarr(nstep+1)
      chi2arr = fltarr(nstep+1)
      chi2arr[0] = oldchi2
      newscale = 1.0
      newrotation = 0.0
      newx = fiberoffx
      newy = fiberoffy 
      npoint=0
      badstep_ct=0
      for i=0,nstep-1 do begin  
         if i mod 1000 eq 0 then print,i
         tmpscale = newscale*10^(scalestep[i]*stepsize[2])
         tmprotation = newrotation + anglestep[i]*stepsize[1]*!dtor
         rotmatrix = [[cos(tmprotation),-sin(tmprotation)],[sin(tmprotation),cos(tmprotation)]]
         offsetx=offsetx0*tmpscale
         offsety=offsety0*tmpscale
         xy2= [[offsetx],[offsety]]#rotmatrix
         offsetx = xy2[*,0]
         offsety = xy2[*,1]
         dx = xstep[i]*stepsize[0]
         dy = ystep[i]*stepsize[0]
         dxarr = (newx+dx)#(fltarr(n_wavearr)+1)-(fltarr(nfiber)+1)#offsetx
         dyarr = (newy+dy)#(fltarr(n_wavearr)+1)-(fltarr(nfiber)+1)#offsety
         distarr = sqrt(dxarr*dxarr+dyarr*dyarr)
         
         newchi2 = matchbundleratio2(fluxratio,fluxratioerr,reffiber,distarr,cenwave=centerwave,dwave=dwave,usemask=innermask,weight=weight,ivarweight=ivarweight,lambda=lambda)
         if newchi2 lt oldchi2 or randomu(seed,1) lt exp(oldchi2-newchi2) then begin
             npoint += 1
       	     xposarr[npoint] = xposarr[npoint-1]+dx
  	     yposarr[npoint] = yposarr[npoint-1]+dy
             scalearr[npoint] = tmpscale
             rotationarr[npoint] = tmprotation
	     chi2arr[npoint] = newchi2
	     oldchi2 = newchi2
             newx = newx+dx
	     newy = newy+dy 
             newscale = tmpscale
             newrotation = tmprotation
         endif else begin
             badstep_ct +=1
             if badstep_ct gt 300 then begin
                stepsize = stepsize*0.6
                if stepsize[0] lt 0.003 and stepsize[1] lt 0.25 and stepsize[2] lt 0.004 then break
                badstep_ct=0
             endif 
         endelse
      endfor
      print,'Iteration stopped on step ',i
      chi2arr = chi2arr[0:npoint]
      xposarr = xposarr[0:npoint]
      yposarr = yposarr[0:npoint]
      scalearr = scalearr[0:npoint]
      rotationarr = rotationarr[0:npoint]
      tmp = min(chi2arr,ind)
      xbest = xposarr[ind]
      ybest = yposarr[ind]
      scale = scalearr[ind]
      rotation = rotationarr[ind]
      fchi2 = chi2arr[ind]
      return,fchi2
end


function modelthrupt, distarr, wavelength
common com_fluxbundle2, psfplane,xarr,wavearr

;   if n_elements(distarr) ne n_elements(wavelength) then message,'Error: distarr and wavelength must have the same number of elements.'
;   xpos = interpol(findgen(n_elements(xarr)),xarr,xoff)
;   ypos = interpol(findgen(n_elements(yarr)),yarr,yoff)
   nfiber = (size(distarr,/dimen))[0]
   if (size(distarr,/dimen))[1] ne n_elements(wavelength) then message,"Error:distarr's first dimension must agree with the size of wavelength."
   xpos = interpol(findgen(n_elements(xarr)),xarr,distarr)
   wavepos = interpol(findgen(n_elements(wavearr)),wavearr,wavelength)
   for i=0,nfiber-1 do begin
     covfn = interpolate(psfplane,xpos[i,*],wavepos)
     if i eq 0 then totcovfn = covfn else totcovfn+=covfn
   endfor
   return,totcovfn
end

pro bundlefluxrationew, plate, mjd, expseq,totresult;, newloglam,mratfit,factor=factor,scalebest=scalebest,rotationbest=rotationbest,xbest=xbest,ybest=ybest
;common com_fluxbundle, cube,xarr,yarr,wavearr
common com_fluxbundle2, psfplane,xarr,wavearr

    gcam_redux_dir = getenv('GCAM_REDUX')
    manga_redux_dir = getenv('MANGA_SPECTRO_REDUX')
    summary = mrdfits(gcam_redux_dir+'gpsfcat.fits',1)
      
    mjdstr = string(mjd,format='(i0.0)')
    dir = manga_redux_dir+string(plate,format='(i0.0)')+'/'+mjdstr+'/'

    nexp = n_elements(expseq)
    allexp = where(summary.plate eq plate and summary.mjd eq mjd,nallexp)
    if max(expseq) gt nallexp-1 or min(expseq) lt 0 then message,'Error: wrong exposure sequence index specified.'
    exposures = summary[allexp[expseq]].exposure  
    observations = summary[allexp[expseq]]

    for i=0,nexp-1 do begin
;       t = where(summary.exposure eq exposures[i], have) 
;       if have ne 1 then message,'Error: '+string(have, format='(i0.0)')+' observation found. I am expecting just one.'

    ;   observation=summary[t]
       expstr = string(exposures[i],format='(i8.8)')
       bfile = file_search(dir,'spFrame-b1-'+expstr+'.*',count=ct_b)
       rfile = file_search(dir,'spFrame-r1-'+expstr+'.*',count=ct_r)

       if ct_b ne 1 or ct_r ne 1 then message, 'Error: Found '+string(ct_b)+' spFrame file for the blue side and '+string(ct_r)+' spFrame file for the red side. I am expecting 1 and 1.'

       mlframe_read,bfile,objflux=bflux,objivar=bivar,mask=bmask,loglam=bloglam,dispimg=bdispimg,plugmap=plugmap,wset=bwset
       mlframe_read,rfile,objflux=rflux,objivar=rivar,mask=rmask,loglam=rloglam,dispimg=rdispimg,plugmap=plugmap,wset=rwset

       correct_dlam,bflux,bivar,bwset,dlam=dloglam
       correct_dlam,bdispimg,0,bwset,dlam=dloglam,/inverse
       correct_dlam,rflux,rivar,rwset,dlam=dloglam
       correct_dlam,rdispimg,0,rwset,dlam=dloglam,/inverse
       if i eq 0 then begin
          bsize=size(bflux,/dimen)
	  rsize=size(rflux,/dimen)
	  totbflux=fltarr([bsize,nexp]) & totbivar=totbflux 
	  totbloglam=dblarr([bsize,nexp]) & totbdispimg=totbloglam
	  totbmask=lonarr([bsize,nexp])
	  totrflux=fltarr([rsize,nexp]) & totrivar=totrflux 
	  totrloglam=dblarr([rsize,nexp]) & totrdispimg=totrloglam
	  totrmask=lonarr([rsize,nexp])
       endif
       totbflux[*,*,i] = bflux
       totbivar[*,*,i] = bivar
       totbdispimg[*,*,i] = bdispimg
       totbloglam[*,*,i] = bloglam
       totbmask[*,*,i] = bmask
       totrflux[*,*,i] = rflux
       totrivar[*,*,i] = rivar
       totrdispimg[*,*,i] = rdispimg
       totrloglam[*,*,i] = rloglam
       totrmask[*,*,i] = rmask
    endfor ; for 3 frames, this will use 200MB of memory.

    slmread=mlreadslm(plate,mjd,slitmap)
;    bundleid = 'ma001'
  all=where(strmatch(slitmap.ifuname,'ma*'))
  uu =uniq(slitmap[all].ifuname,sort(slitmap[all].ifuname))
  bundleids=slitmap[all[uu]].ifuname
  if plate eq 6870 then bundleids=bundleids[where(bundleids ne 'ma001')]
  nbundle=n_elements(bundleids)
  for ib=0,nbundle-1 do begin
    bind = where(slitmap.ifuname eq bundleids[ib],nfiber)
    ss = sort(slitmap[bind].fnum)
    bind = bind[ss]

    status=mlreadbm(bundleids[ib],mjd,bmap)
; check bmap match with slitmap[bind] in fiberid
    if min(slitmap[bind].fnum-bmap.ise) ne 0 or max(slitmap[bind].fnum-bmap.ise) ne 0 then message,'Error: bmap and slitmap[bind] do not have matching fnum! some sorting is needed.'

    dloglam=1.d-4
    maxlam = 10600.
    minlam = 3400.
    nlam = ceil((alog10(maxlam)-alog10(minlam))/dloglam)

    newloglam = alog10(minlam)+dloglam*dindgen(nlam)

;    idsky = where(strmatch(plugmap.objtype,'SKY*'))
;    p=where(slitmap[plugmap[idsky].fiberid-1].fsize eq 2.0)
;    idsky2a = idsky[p]
;    combine1fiber,bloglam[*,idsky2a],bflux[*,idsky2a],bivar[*,idsky2a],finalmask=bmask[*,idsky2a],newloglam=newloglam,newflux=bnewflux,newivar=bnewivar,andmask=bnewmask
;    combine1fiber,rloglam[*,idsky2a],rflux[*,idsky2a],rivar[*,idsky2a],finalmask=rmask[*,idsky2a],newloglam=newloglam,newflux=rnewflux,newivar=rnewivar,andmask=rnewmask

;    skyivar = bnewivar+rnewivar
;    skyflux = (bnewflux*bnewivar+rnewflux*rnewivar)/skyivar
;    bad = where(skyivar eq 0,nbad)
;    if nbad gt 0 then skyflux[bad] = (bnewflux[bad]+rnewflux[bad])/2.

    lam0 = 3500.
;    dwave = 1000.
;    nwave = fix((10000-lam0)/dwave)
;    centerwave =lam0+(findgen(nwave)+0.5)*dwave
    wave1 = [3500.,4000.,4500.,5000.,5500.,6500.,7500.,9000.]
    wave2 = [4000.,4500.,5000.,5500.,6500.,7500.,9000.,10500.]
    nwave = n_elements(wave1)
    centerwave = (wave1+wave2)/2.
    dwave = wave2-wave1
    skylevel = fltarr(nwave)
    skyerr = fltarr(nwave)

;    for j=0,nwave-1 do begin
;       lam1 = alog10(lam0+j*dwave)
;       lam2 = alog10(lam0+(j+1)*dwave)
;       indwave=where(newloglam gt lam1 and newloglam lt lam2 and skyivar ne 0,ngood)
;       if ngood eq 0 then message,'No good flux obtained.'
;       skylevel[j] = median(skyflux[indwave]) 
;       skyerr[j] = medianerr(skyflux[indwave])
;    endfor
    totfluxarr = fltarr(nlam,nfiber,nexp)
    totivararr = fltarr(nlam,nfiber,nexp)
    totdisparr = fltarr(n_elements(newloglam),nfiber,nexp)
    for ie=0,nexp-1 do begin
      for i=0,nfiber-1 do begin
        combine1fiber,totbloglam[*,bind[i],ie],totbflux[*,bind[i],ie],totbivar[*,bind[i],ie],finalmask=totbmask[*,bind[i],ie],indisp=totbdispimg[*,bind[i],ie],newloglam=newloglam,newflux=bnewflux,newivar=bnewivar,andmask=bnewmask,newdisp=bnewdisp
        combine1fiber,totrloglam[*,bind[i],ie],totrflux[*,bind[i],ie],totrivar[*,bind[i],ie],finalmask=totrmask[*,bind[i],ie],indisp=totrdispimg[*,bind[i],ie],newloglam=newloglam,newflux=rnewflux,newivar=rnewivar,andmask=rnewmask,newdisp=rnewdisp
        totivararr[*,i,ie] = bnewivar+rnewivar
        totfluxarr[*,i,ie] = (bnewflux*bnewivar+rnewflux*rnewivar)/totivararr[*,i,ie]
        totdisparr[*,i,ie] = (bnewdisp*bnewivar+rnewdisp*rnewivar)/totivararr[*,i,ie]
        bad = where(totivararr[*,i,ie] eq 0,nbad)
        if nbad gt 0 then begin 
           totfluxarr[bad,i,ie] = (bnewflux[bad]+rnewflux[bad])/2.
 	   totdisparr[bad,i,ie] = (bnewdisp[bad]+rnewdisp[bad])/2.
        endif
      endfor
    endfor
    fiberflux = total(totfluxarr*totivararr,1)/total(totivararr,1)
;    fiberflux = total(fluxarr,1)
;    ss = sort(total(fiberflux,2))
    tmp = max(total(fiberflux,2),reffiber)
;    reffiber = ss[nfiber-1]
    if slitmap[bind[reffiber]].fnum ne (nfiber+1)/2 then begin
        print,'Warning: the brightest fiber is not the central fiber. Correct?'
        stop
    endif
;    skyresi = median(fluxarr[*,ss[0:6]],dimen=2)

    platescale = 3.62730/60. ; (mm/arcsec) as given in Gunn et al. 2006
    fiberx = bmap.x0mm/platescale
    fibery = bmap.y0mm/platescale
;    innermask = (sqrt(fiberx*fiberx+fibery*fibery) lt 3) and bmap.gbu ge 0

    middlefiber = fix(nfiber/2)
    fiberx0 = fiberx[middlefiber]
    fibery0 = fibery[middlefiber]
    innermask = (sqrt((fiberx-fiberx0)^2+(fibery-fibery0)^2) lt 3) and bmap.gbu ge 0
    innerfiber = where(innermask,ninner)
    if innermask[reffiber] eq 0 then message,'Error: reffiber is not one of the inner 7 fibers.'

    fnum = slitmap[bind].fnum
;innermask = (fnum ge 84 and fnum le 85) or (fnum ge 101 and fnum le 105)
;stop

    if nfiber gt 19 then begin 
      inner7 = where(innermask)
      uu = uniq(slitmap[bind[inner7]].blockid,sort(slitmap[bind[inner7]].blockid))
      uublockid = slitmap[bind[inner7[uu]]].blockid
      for i=0,n_elements(uu)-1 do begin
        id = where(slitmap[bind].blockid eq uublockid[i])
	if i eq 0 then totid = id else totid=[totid,id]
      endfor
    endif else totid=indgen(19)
    ntotid = n_elements(totid)

    fluxarr = fltarr(nwave,nfiber,nexp)
    fluxerr = fltarr(nwave,nfiber,nexp)
    sbfluxarr = totfluxarr*0.0
    sbivararr = totivararr*0.0
    for ie=0,nexp-1 do begin
      ss = sort(fiberflux[totid,ie])
      skyin = bind[totid[ss[0:ntotid-1-12]]]
      combine1fiber,totbloglam[*,skyin,ie],totbflux[*,skyin,ie],totbivar[*,skyin,ie],finalmask=totbmask[*,skyin,ie],newloglam=newloglam,newflux=bnewflux,newivar=bnewivar,andmask=bnewmask
      combine1fiber,totrloglam[*,skyin,ie],totrflux[*,skyin,ie],totrivar[*,skyin,ie],finalmask=totrmask[*,skyin,ie],newloglam=newloglam,newflux=rnewflux,newivar=rnewivar,andmask=rnewmask

      skyivar2 = bnewivar+rnewivar
      skyflux2 = (bnewflux*bnewivar+rnewflux*rnewivar)/skyivar2
      bad = where(skyivar2 eq 0,nbad)
      if nbad gt 0 then skyflux2[bad] = (bnewflux[bad]+rnewflux[bad])/2.

      sbfluxarr[*,*,ie] = totfluxarr[*,*,ie]-skyflux2#(fltarr(nfiber)+1)
      skyivar_2d = skyivar2#(fltarr(nfiber)+1)
      good = where(totivararr[*,*,ie] ne 0.0 and skyivar_2d ne 0.0)
      tmpivar = sbivararr[*,*,ie]
      tmpivar[good] = 1/(1/(totivararr[*,*,ie])[good]+1/skyivar_2d[good])
      sbivararr[*,*,ie]=tmpivar
      for j=0,nwave-1 do begin
         lam1 = alog10(wave1[j])
         lam2 = alog10(wave2[j])
         indwave=where(newloglam gt lam1 and newloglam lt lam2,ngood)
         if ngood eq 0 then message,'No good flux obtained.'
         fluxarr[j,*,ie] = total(sbfluxarr[indwave,*,ie]*sbivararr[indwave,*,ie],1)/total(sbivararr[indwave,*,ie],1)
         fluxerr[j,*,ie] = sqrt(1/total(sbivararr[indwave,*,ie],1))
      endfor
;    fluxarr = fluxarr - skyresi#(fltarr(nfiber+1))
    endfor

;    tmp = max(total(fluxarr,1)*innermask,reffiber)

;    ind1 = where(10^newloglam gt 3700 and 10^newloglam lt 3750)
;    snr=median(sbfluxarr[ind1]*sqrt(sbivararr[ind1]))/sqrt(3725*(10^(1.e-4)-1)) ; S/N per Angstrom
;stop
    comivar=total(sbivararr[*,reffiber,*],3)
;    comflux=total(sbfluxarr[*,reffiber,*]*sbivararr[*,reffiber,*],3)/comivar
    comflux=total(sbfluxarr[*,reffiber,*],3)/nexp
    ivarproduct=sbivararr[*,reffiber,0]
    for i=1,nexp-1 do ivarproduct=ivarproduct*sbivararr[*,reffiber,i]
    t=where(ivarproduct eq 0,nbad)
    if nbad gt 0 then comivar[t] = 0.0

    tmp = typingmodule(comflux,newloglam,comivar,totdisparr[*,reffiber,0],bind[reffiber],observations[0],plugmap,kindx=kindx,modflux=modflux)

    thru=mrdfits('meanthrupt.fits',1)
    init_thru=interpol(thru.thrupt,thru.loglam,newloglam)
    

    for ie=0,nexp-1 do begin
      fluxratio = fluxarr[*,*,ie]/(fluxarr[*,reffiber,ie]#(fltarr(nfiber)+1))
      fluxratioerr = abs(fluxratio)*sqrt((fluxerr[*,*,ie]/fluxarr[*,*,ie])^2+(fluxerr[*,reffiber,ie]/fluxarr[*,reffiber,ie])^2#(fltarr(nfiber)+1))

;    stop
      tmp=min(abs(centerwave-5300),refwave)

; interpolate to get the expected flux per fiber at this dither position.    

;    staroffx = (191.10996-191.11030)*3600.
;    staroffy = (17.510236-17.507112)*3600.
      staroffx = 0.0
      staroffy = 0.0
      fiberoffx = bmap.xpmm/platescale+observations[ie].dcnra-staroffx
      fiberoffy = bmap.ypmm/platescale+observations[ie].dcndec-staroffy


      nsize = 1
      factorarr=1+(findgen(nsize)-nsize/2)*0.05
;    totchi2arr = fltarr(nsize)
;    scalearr = findgen(15)/10.+0.6
;    scalearr = [0.9,1.0,1.1]
;    nscl = n_elements(scalearr) 
;    chi2sizearr=fltarr(nscl,nsize); chi2sizearr = fltarr(nwave,nsize)
      xcen = fltarr(nsize) ; xcen_arr = fltarr(nwave,nsize)
      ycen = fltarr(nsize) ; ycen_arr = fltarr(nwave,nsize)
      scale = fltarr(nsize)
      rotation = fltarr(nsize)
      fchi2 = fltarr(nsize)
      for is=0,nsize-1 do begin
        factor = factorarr[is]
        params=observations[ie].params
        params[2:3] = params[2:3]*factor
        params[0:1] = params[0:1]/factor^2
        fcpsf1d,params,observations[ie].ha-(slitmap[bind[0]].ra-slitmap[bind[0]].cenra),plugmap[bind[0]].dec,plugmap[bind[0]].xfocal,plugmap[bind[0]].yfocal,psfplane,wavearr=wavearr,xarr=xarr,magnify=2,/resetwave

        if is eq 0 then status=mldar((observations[ie].ha-(slitmap[bind[0]].ra-slitmap[bind[0]].cenra))/15.,plugmap[bind[0]].dec,wavearr,parangle,offsetx0,offsety0,waveref=centerwave[refwave],raobj=slitmap[bind[0]].ra,racen=slitmap[0].cenra,deccen=slitmap[0].cendec,/distort)
        n_wavearr = n_elements(wavearr)
;      fchi2[is]=mcmcbundleratio(fluxratio,$
;                  fluxratioerr,offsetx0=offsetx0,$
;                  offsety0=offsety0,fiberoffx=fiberoffx,$
;                  fiberoffy=fiberoffy,nfiber=nfiber,$
;                  n_wavearr=n_wavearr,reffiber=reffiber,$
;                  centerwave=centerwave,dwave=dwave,innermask=innermask,$
;                  weight=modflux*init_thru,ivarweight=sbivararr,$
;                  lambda=10^newloglam,xbest=xbest,ybest=ybest,scale=scalebest,$
;                  rotation=rotationbest,nstep=2000)
        fchi2[is]=mcmcbundleratio(fluxratio[*,innerfiber],$
                  fluxratioerr[*,innerfiber],offsetx0=offsetx0,$
                  offsety0=offsety0,fiberoffx=fiberoffx[innerfiber],$
                  fiberoffy=fiberoffy[innerfiber],nfiber=ninner,$
                  n_wavearr=n_wavearr,reffiber=(where(innerfiber eq reffiber))[0],$
                  centerwave=centerwave,dwave=dwave,innermask=intarr(ninner)+1,$
                  weight=modflux*init_thru,ivarweight=sbivararr[*,innerfiber],$
                  lambda=10^newloglam,xbest=xbest,ybest=ybest,scale=scalebest,$
                  rotation=rotationbest,nstep=2000)

        xcen[is] = xbest
        ycen[is] = ybest
        scale[is] = scalebest
        rotation[is] = rotationbest
;      offsetx=offsetx0*scale[is]
;      offsety=offsety0*scale[is]
;      rotmatrix = [[cos(rotation[is]),-sin(rotation[is])],[sin(rotation[is]),cos(rotation[is])]]
;      xy2= [[offsetx],[offsety]]#rotmatrix
;      offsetx = xy2[*,0]
;      offsety = xy2[*,1]
;      dxarr = (fiberoffx+xcen[is])#(fltarr(n_wavearr)+1)-(fltarr(nfiber)+1)#offsetx
;      dyarr = (fiberoffy+ycen[is])#(fltarr(n_wavearr)+1)-(fltarr(nfiber)+1)#offsety
;      distarr = sqrt(dxarr*dxarr+dyarr*dyarr)
;      chi2 = matchbundleratio2(fluxratio,fluxratioerr,reffiber,distarr,cenwave=centerwave,dwave=dwave,usemask=innermask,weight=modflux*init_thru,ivarweight=sbivararr,lambda=10^newloglam,predict=predict)
;stop
;      endfor
;      chi2sizearr[*,is]=fchi2
;      xcen_arr[*,is] = xcen
;      ycen_arr[*,is] = ycen
;      totchi2arr[is] = fchi2 ; total(fchi2)
;      chi2 = matchbundleratio(fluxratio[3,*],fluxratioerr[3,*],reffiber,fiberoffx+xcen[3],fiberoffy+ycen[3],cenwave=centerwave[3],dwave=dwave[3],usemask=innermask,weight=modflux*init_thru,ivarweight=sbivararr,lambda=10^newloglam,predict=predict)
      endfor

      if nsize ge 5 then begin
        res = svdfit(factorarr,fchi2,3)
        factor = -res[1]/(2*res[2])

;   tmp = min(fchi2,ind_s)
;   xcen_final = xcen_arr[ind_s]
;   ycen_final = ycen_arr[ind_s]
       
        params=observations[ie].params
        params[2:3] = params[2:3]*factor
        params[0:1] = params[0:1]/factor^2
        fcpsf1d,params,observations[ie].ha-(slitmap[bind[0]].ra-slitmap[bind[0]].cenra),plugmap[bind[0]].dec,plugmap[bind[0]].xfocal,plugmap[bind[0]].yfocal,psfplane,wavearr=wavearr,xarr=xarr,magnify=2,/resetwave
        bestchi2=mcmcbundleratio(fluxratio[*,innerfiber],fluxratioerr[*,innerfiber], $
           offsetx0=offsetx0,offsety0=offsety0,fiberoffx=fiberoffx[innerfiber],$
           fiberoffy=fiberoffy[innerfiber],nfiber=ninner,n_wavearr=n_wavearr,$
           reffiber=(where(innerfiber eq reffiber))[0],centerwave=centerwave,$
           dwave=dwave,innermask=intarr(ninner)+1,weight=modflux*init_thru,$
           ivarweight=sbivararr[*,innerfiber],lambda=10^newloglam,xbest=xbest,$
           ybest=ybest,scale=scalebest,rotation=rotationbest,nstep=3000)
        if bestchi2 gt min(fchi2) then begin
          bestchi2 = min(fchi2,ind_s)
          factor = factorarr[ind_s]
          params=observations[ie].params
          params[2:3] = params[2:3]*factor
          params[0:1] = params[0:1]/factor^2
          fcpsf1d,params,observations[ie].ha-(slitmap[bind[0]].ra-slitmap[bind[0]].cenra),plugmap[bind[0]].dec,plugmap[bind[0]].xfocal,plugmap[bind[0]].yfocal,psfplane,wavearr=wavearr,xarr=xarr,magnify=2,/resetwave
          scalebest = scale[ind_s]
          rotationbest=rotation[ind_s]
          xbest= xcen[ind_s]
          ybest= ycen[ind_s]
        endif 
      endif else begin
        bestchi2 = min(fchi2,ind_s)
        factor = factorarr[ind_s]
        scalebest = scale[ind_s]
        rotationbest=rotation[ind_s]
        xbest= xcen[ind_s]
        ybest= ycen[ind_s]
        params=observations[ie].params
        params[2:3] = params[2:3]*factor
        params[0:1] = params[0:1]/factor^2
        fcpsf1d,params,observations[ie].ha-(slitmap[bind[0]].ra-slitmap[bind[0]].cenra),plugmap[bind[0]].dec,plugmap[bind[0]].xfocal,plugmap[bind[0]].yfocal,psfplane,wavearr=wavearr,xarr=xarr,magnify=2,/resetwave
      endelse  

      offsetx = offsetx0*scalebest
      offsety = offsety0*scalebest
      rotmatrix = [[cos(rotationbest),-sin(rotationbest)],[sin(rotationbest),cos(rotationbest)]]
      xy2 = [[offsetx],[offsety]]#rotmatrix
      offsetx = xy2[*,0]
      offsety = xy2[*,1]
      dxarr = (fiberoffx+xbest)#(fltarr(n_wavearr)+1)-(fltarr(nfiber)+1)#offsetx
      dyarr = (fiberoffy+ybest)#(fltarr(n_wavearr)+1)-(fltarr(nfiber)+1)#offsety
      distarr = sqrt(dxarr*dxarr+dyarr*dyarr)

      snratio = total(sbfluxarr[*,innerfiber,ie]*sqrt(sbivararr[*,innerfiber,ie]),1)
      tmp = max(snratio,ind)
      bright= where(snratio gt 0.5*max(snratio),nbright)
      covfn = modelthrupt(distarr[innerfiber[bright],*],wavearr)
;stop
;      covfn= modelthrupt(distarr[reffiber,*],wavearr)
;   chi2 = matchbundleratio2(fluxratio,fluxratioerr,reffiber,distarr,cenwave=centerwave,dwave=dwave,usemask=innermask,weight=modflux*init_thru,ivarweight=sbivararr,lambda=10^newloglam,predict=predict)
      covfnall = interpol(covfn,wavearr,10^newloglam,/spline)

      k = where(innermask)
;   dist = sqrt((fiberoffx+xcen_final)^2+(fiberoffy+ycen_final)^2)
;   waveind = interpol(findgen(n_elements(wavearr)),wavearr,centerwave[3])
;   xx = findgen(400)/100.
;   xind = interpol(findgen(n_elements(xarr)),xarr,xx)
;   yind = interpol(findgen(n_elements(yarr)),yarr,0.0)
;   psf = interpolate(cube,xind,yind,waveind,/grid)
;   plot,xx,psf

;   norm = interpol(psf,xx,dist[reffiber])
;   oplot,dist[k],fluxratio[3,k]*norm,ps=1
;stop
;         id1 = where(chi2arr gt tmp+0.9 and chi2arr lt tmp+1.1,nid1)
;         xid1rec[0:(nid1<100)-1,idw]=xposarr[id1]
;         yid1rec[0:(nid1<100)-1,idw]=yposarr[id1]

;   print,observation.exposure,centerwave[refwave],xcen[refwave],ycen[refwave]
;   status=mldar((observation.ha-(slitmap[bind[0]].ra-slitmap[bind[0]].cenra))/15.,plugmap[bind[0]].dec,centerwave,parangle,offsetx,offsety,waveref=centerwave[refwave],raobj=slitmap[bind[0]].ra,racen=slitmap[0].cenra,deccen=slitmap[0].cendec,/distort)
;   xrange=minmax([-xcen+xcen[refwave],offsetx])
;   yrange=minmax([-ycen+ycen[refwave],offsety])
;   psplot,'dar_off.ps',/square,/color
;   !x.margin=[10,3]
;   plot,[-xcen+xcen[refwave],offsetx],[-ycen+ycen[refwave],offsety],ps=1,/nodata,xtitle='d_RA (arcsec)',ytitle='d_Dec (arcsec)',xcharsize=1.3,ycharsize=1.3,title=bundleid
;   oplot,-xcen+xcen[refwave],-ycen+ycen[refwave],ps=-1,symsize=2
;   oplot,offsetx,offsety,ps=-4,linest=2,symsize=2,color=249
;   legend,['Measured DAR','Theoretical DAR'],linest=[0,2],psym=[1,4],symsize=[2,2],/bottom,box=0
;   legend,[string(observation.exposure,format='(i8.8)'),bundleid],box=0,/right,/top
;stop
;   for i=0,nwave-1 do begin &$
;     t=where(xid1rec[*,i] ne 0.0) &$
;     oplot,-xid1rec[t,i]+xcen[refwave],-yid1rec[t,i]+ycen[refwave],ps=1,symsize=0.1 &$
;   endfor

;   x1 = interpol(fiberoffx[reffiber]+xcen,centerwave,10^newloglam,/spline)
;   y1 = interpol(fiberoffy[reffiber]+ycen,centerwave,10^newloglam,/spline)
   ; fraction of light going down the central fiber at each wavelength
;   fraction = modelthrupt(x1,y1,10^newloglam)

      data = total(reform(sbfluxarr[*,innerfiber[bright],ie], nlam, nbright),2)
      ivarproduct = product(reform(sbivararr[*,innerfiber[bright],ie],nlam,nbright),2)
      ivar = total(reform(sbivararr[*,innerfiber[bright],ie],nlam,nbright),2)*(ivarproduct gt 0)
      model = modflux*covfnall
      mratio = data/model
      mrativar = ivar*model^2
      mrativar = mrativar *(1- spflux_masklines(newloglam,/stellar))
      bb=where(newloglam lt alog10(6800.))
      b_sset=spflux_bspline(newloglam[bb],mratio[bb],mrativar[bb],everyn=10)
      rr=where(newloglam gt alog10(6800.))
      r_sset=spflux_bspline(newloglam[rr],mratio[rr],mrativar[rr],everyn=1.5)
      mratfit = newloglam*0.0
      mratfit[bb]=bspline_valu(newloglam[bb],b_sset)
      mratfit[rr]=bspline_valu(newloglam[rr],r_sset)

      result={exposure:0L,ifuname:' ',factor:0.0,xbest:0.0,ybest:0.0,$
              scalebest:0.0,rotationbest:0.0,chi2:0.0,model:fltarr(nlam),modflux:fltarr(nlam),$
              data:fltarr(nlam),mratio:fltarr(nlam),mrativar:fltarr(nlam),$
	      newloglam:dblarr(nlam)}
      if n_elements(totresult) eq 0 then begin
	 totresult = replicate(result,nexp*nbundle)
      endif
      result.exposure=exposures[ie]
      result.ifuname=bundleids[ib]
      result.factor=factor
      result.xbest=xbest
      result.ybest=ybest
      result.scalebest=scalebest
      result.rotationbest=rotationbest
      result.chi2 =bestchi2
      result.model=model
      result.modflux=modflux
      result.data =sbfluxarr[*,reffiber,ie]
      result.mratio=mratio
      result.mrativar=mrativar
      result.newloglam=newloglam
      totresult[ie*nbundle+ib]=result
    endfor
  endfor
  mwrfits,totresult,'fluxcal'+string(plate,format='(i0.0)')+'_'+mjdstr+'.fits',/create
end
