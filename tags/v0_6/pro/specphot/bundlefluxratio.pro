function matchbundleratio, bundleratio, ratioerr, reffiber, dxarr,dyarr, cenwave=cenwave,dwave=dwave,predictratio=predictratio,usemask=usemask
common com_fluxbundle, cube,xarr,yarr,wavearr
    if n_elements(dxarr) ne n_elements(dyarr) then message, 'Error: dxarr and dyarr must have the same number of elements.'
    if n_elements(cenwave) ne n_elements(dwave) then message,'Error: cenwave and dwave must have the same number of elements.'
    xpos = interpol(findgen(n_elements(xarr)),xarr,dxarr)
    ypos = interpol(findgen(n_elements(yarr)),yarr,dyarr)
    nfiber = n_elements(dxarr)
    nwave = n_elements(cenwave)
    thrupt=fltarr(nwave,nfiber)
    for i=0,nfiber-1 do begin
      thruptallwave = reform(interpolate(cube,xpos[i],ypos[i],indgen(n_elements(wavearr)),/grid))
      for j=0,nwave-1 do begin
         lam1=cenwave[j]-dwave[j]/2.
	 lam2=cenwave[j]+dwave[j]/2.
	 ind = where(wavearr ge lam1 and wavearr le lam2)
         thrupt[j,i] = mean(thruptallwave[ind])
;	 thrupt[j,i] = median(thruptallwave[ind])
;        thrupt[j,i] = int_tabulated(wavearr[ind],thruptallwave[ind])
      endfor
    endfor
    predictratio=thrupt/(thrupt[*,reffiber]#(fltarr(nfiber)+1))
    useid = where(indgen(nfiber) ne reffiber and usemask and bundleratio gt 0.01,nuse)
    chi2 = total(((bundleratio[*,useid]-predictratio[*,useid])^2)/(ratioerr[*,useid]^2))/(nuse*nwave-1.)
    return,chi2
end

function modelthrupt, xoff,yoff, wavelength
common com_fluxbundle, cube, xarr,yarr,wavearr

   if n_elements(xoff) ne n_elements(yoff) or n_elements(xoff) ne n_elements(wavelength) then message,'Error: xoff, yoff, and wavelength must have the same number of elements.'
   xpos = interpol(findgen(n_elements(xarr)),xarr,xoff)
   ypos = interpol(findgen(n_elements(yarr)),yarr,yoff)
   thrupt = interpolate(cube,xpos,ypos,wavelength)
   return,thrupt
end

pro bundlefluxratio, exposure, bundleid,newloglam,mratfit
common com_fluxbundle, cube,xarr,yarr,wavearr

    gcam_redux_dir = getenv('GCAM_REDUX')
    summary = mrdfits(gcam_redux_dir+'gpsfcat.fits',1)
      
;    readdither,'dithertable.txt',ditherposarr,dcnraarr,dcndecarr

;    p = where(ditherposarr eq ditherpos, havedither)
;    if havedither ne 1 then message,'There is no such dither positions: D'+string(ditherpos,format='(i0.0)')

;    dx = abs(summary.dcnra-dcnraarr[p])
;    dy = abs(summary.dcndec-dcndecarr[p])

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

    mlframe_read,bfile,objflux=bflux,objivar=bivar,mask=bmask,loglam=bloglam,dispimg=bdispimg,plugmap=plugmap,wset=bwset
    mlframe_read,rfile,objflux=rflux,objivar=rivar,mask=rmask,loglam=rloglam,dispimg=rdispimg,plugmap=plugmap,wset=rwset

    correct_dlam,bflux,bivar,bwset,dlam=dloglam
    correct_dlam,bdispimg,0,bwset,dlam=dloglam,/inverse
    correct_dlam,rflux,rivar,rwset,dlam=dloglam
    correct_dlam,rdispimg,0,rwset,dlam=dloglam,/inverse

    slmread=mlreadslm(observation.plate,observation.mjd,slitmap)
;    bundleid = 'ma001'
    bind = where(slitmap.ifuname eq bundleid,nfiber)
    ss = sort(slitmap[bind].fnum)
    bind = bind[ss]

    status=mlreadbm(bundleid,observation.mjd,bmap)
; check bmap match with slitmap[bind] in fiberid
    if min(slitmap[bind].fnum-bmap.ise) ne 0 or max(slitmap[bind].fnum-bmap.ise) ne 0 then message,'Error: bmap and slitmap[bind] do not have matching fnum! some sorting is needed.'

;    params=observation.params
;    factor = 1.0
;    params[2:3] = params[2:3]*factor
;    params[0:1] = params[0:1]/factor^2
;    fcpsf,params,observation.ha-(slitmap[bind[0]].ra-slitmap[bind[0]].cenra),plugmap[bind[0]].dec,plugmap[bind[0]].xfocal,plugmap[bind[0]].yfocal,cube,wavearr=wavearr,xarr=xarr,yarr=yarr,magnify=2,/nodar

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
    fluxarr = fltarr(nwave,nfiber)
    fluxerr = fltarr(nwave,nfiber)
    totfluxarr = fltarr(nlam,nfiber)
    totivararr = fltarr(nlam,nfiber)
    totdisparr = fltarr(n_elements(newloglam),nfiber)
    for i=0,nfiber-1 do begin
       combine1fiber,bloglam[*,bind[i]],bflux[*,bind[i]],bivar[*,bind[i]],finalmask=bmask[*,bind[i]],indisp=bdispimg[*,bind[i]],newloglam=newloglam,newflux=bnewflux,newivar=bnewivar,andmask=bnewmask,newdisp=bnewdisp
       combine1fiber,rloglam[*,bind[i]],rflux[*,bind[i]],rivar[*,bind[i]],finalmask=rmask[*,bind[i]],indisp=rdispimg[*,bind[i]],newloglam=newloglam,newflux=rnewflux,newivar=rnewivar,andmask=rnewmask,newdisp=rnewdisp
       totivararr[*,i] = bnewivar+rnewivar
       totfluxarr[*,i] = (bnewflux*bnewivar+rnewflux*rnewivar)/totivararr[*,i]
       totdisparr[*,i] = (bnewdisp*bnewivar+rnewdisp*rnewivar)/totivararr[*,i]
       bad = where(totivararr[*,i] eq 0,nbad)
       if nbad gt 0 then begin 
          totfluxarr[bad,i] = (bnewflux[bad]+rnewflux[bad])/2.
	  totdisparr[bad,i] = (bnewdisp[bad]+rnewdisp[bad])/2.
       endif
    endfor
    fiberflux = total(totfluxarr*totivararr,1)/total(totivararr,1)
;    fiberflux = total(fluxarr,1)
    ss = sort(fiberflux)
;    skyresi = median(fluxarr[*,ss[0:6]],dimen=2)

    platescale = 3.62730/60. ; (mm/arcsec) as given in Gunn et al. 2006
    fiberx = bmap.x0mm/platescale
    fibery = bmap.y0mm/platescale
    innermask = (sqrt(fiberx*fiberx+fibery*fibery) lt 3) and bmap.gbu ge 0


    if nfiber gt 19 then begin 
      inner7 = where(innermask)
      uu = uniq(slitmap[bind[inner7]].blockid,sort(slitmap[bind[inner7]].blockid))
      uublockid = slitmap[bind[inner7[uu]]].blockid
      for i=0,n_elements(uu)-1 do begin
        id = where(slitmap[bind].blockid eq uublockid[i])
	if i eq 0 then totid = id else totid=[totid,id]
      endfor
      ntotid = n_elements(totid)
      ss = sort(fiberflux[totid])
      skyin = bind[totid[ss[0:ntotid-1-12]]]
    endif else skyin = bind[ss[0:6]]

    combine1fiber,bloglam[*,skyin],bflux[*,skyin],bivar[*,skyin],finalmask=bmask[*,skyin],newloglam=newloglam,newflux=bnewflux,newivar=bnewivar,andmask=bnewmask
    combine1fiber,rloglam[*,skyin],rflux[*,skyin],rivar[*,skyin],finalmask=rmask[*,skyin],newloglam=newloglam,newflux=rnewflux,newivar=rnewivar,andmask=rnewmask

    skyivar2 = bnewivar+rnewivar
    skyflux2 = (bnewflux*bnewivar+rnewflux*rnewivar)/skyivar2
    bad = where(skyivar2 eq 0,nbad)
    if nbad gt 0 then skyflux2[bad] = (bnewflux[bad]+rnewflux[bad])/2.

    sbfluxarr = totfluxarr-skyflux2#(fltarr(nfiber)+1)
    skyivar_2d = skyivar2#(fltarr(nfiber)+1)
    good = where(totivararr ne 0.0 and skyivar_2d ne 0.0)
    sbivararr = totivararr*0.0
    sbivararr[good] = 1/(1/totivararr[good]+1/skyivar_2d[good])
    for j=0,nwave-1 do begin
       lam1 = alog10(wave1[j])
       lam2 = alog10(wave2[j])
       indwave=where(newloglam gt lam1 and newloglam lt lam2,ngood)
       if ngood eq 0 then message,'No good flux obtained.'
       fluxarr[j,*] = total(sbfluxarr[indwave,*]*sbivararr[indwave,*],1)/total(sbivararr[indwave,*],1)
       fluxerr[j,*] = sqrt(1/total(sbivararr[indwave,*],1))
    endfor
;    fluxarr = fluxarr - skyresi#(fltarr(nfiber+1))

    tmp = max(total(fluxarr,1)*innermask,reffiber)
;stop

    tmp = typingmodule(sbfluxarr[*,reffiber],newloglam,sbivararr[*,reffiber],totdisparr[*,reffiber],bind[reffiber],observation,plugmap,kindx=kindx,modflux=modflux)

    fluxratio = fluxarr/(fluxarr[*,reffiber]#(fltarr(nfiber)+1))
    fluxratioerr = abs(fluxratio)*sqrt((fluxerr/fluxarr)^2+(fluxerr[*,reffiber]/fluxarr[*,reffiber])^2#(fltarr(nfiber)+1))

;    stop

; interpolate to get the expected flux per fiber at this dither position.    

    fiberoffx = bmap.xpmm/platescale+observation.dcnra
    fiberoffy = bmap.ypmm/platescale+observation.dcndec


;    mldar((observation.ha-(slitmap[bind[0]].ra-slitmap[bind[0]].cenra))/15.,slitmap[bind[0]].dec,centerwave,parangle,dar_dx,dar_dy,waveref=centerwave[refwave])
    nsize = 7
    factorarr=1+(findgen(nsize)-nsize/2)*0.05
    totchi2arr = fltarr(nsize)
    chi2sizearr = fltarr(nwave,nsize)
    xcen_arr = fltarr(nwave,nsize)
    ycen_arr = fltarr(nwave,nsize)
    for is=0,nsize-1 do begin
      factor = factorarr[is]
      params=observation.params
      params[2:3] = params[2:3]*factor
      params[0:1] = params[0:1]/factor^2
      fcpsf,params,observation.ha-(slitmap[bind[0]].ra-slitmap[bind[0]].cenra),plugmap[bind[0]].dec,plugmap[bind[0]].xfocal,plugmap[bind[0]].yfocal,cube,wavearr=wavearr,xarr=xarr,yarr=yarr,magnify=2,/nodar,/resetwave

      xcen = fltarr(nwave)
      ycen = fltarr(nwave)
      fchi2 = fltarr(nwave)
      xid1rec = fltarr(100,nwave)
      yid1rec = fltarr(100,nwave)
      for idw=0,nwave-1 do begin
         chi2 = matchbundleratio(fluxratio[idw,*],fluxratioerr[idw,*],reffiber,fiberoffx,fiberoffy,cenwave=centerwave[idw],dwave=dwave[idw],usemask=innermask)
         nstep =1000
         steps=randomu(seed,nstep*2)
         rr = steps[0:nstep-1]
         theta = steps[nstep:2*nstep-1]*!pi*2
         xstep = rr*cos(theta)
         ystep = rr*sin(theta)
         stepsize = 0.1
         oldchi2=chi2
         newx = fiberoffx
         newy = fiberoffy
         xposarr = fltarr(nstep+1)
         yposarr = fltarr(nstep+1)
         chi2arr = fltarr(nstep+1)
         chi2arr[0] = oldchi2
         npoint=0
         for i=0,nstep-1 do begin
;           print,'MCMC step ',i
           dx = xstep[i]*stepsize
           dy = ystep[i]*stepsize
           newchi2 = matchbundleratio(fluxratio[idw,*],fluxratioerr[idw,*],reffiber,newx+dx,newy+dy,cenwave=centerwave[idw],dwave=dwave[idw],usemask=innermask)
           if randomu(seed,1) lt exp(oldchi2-newchi2) then begin
             npoint += 1
       	     xposarr[npoint] = xposarr[npoint-1]+dx
  	     yposarr[npoint] = yposarr[npoint-1]+dy
	     chi2arr[npoint] = newchi2
	     oldchi2 = newchi2
             newx = newx+dx
	     newy = newy+dy 
           endif
         endfor
         chi2arr = chi2arr[0:npoint]
         xposarr = xposarr[0:npoint]
         yposarr = yposarr[0:npoint]
;         print,'X,Y:',fiberoffx[9]+xcen[idw], fiberoffy[9]+ycen[idw],format='(a,f6.2,f6.2)'
;         chi2 = matchbundleratio(fluxratio[idw,*],fluxratioerr[idw,*],reffiber,fiberoffx+xposarr[ind],fiberoffy+yposarr[ind],cenwave=centerwave[idw],dwave=dwave[idw],usemask=innermask,predict=predict)
         tmp = min(chi2arr,ind)
         xcen[idw] = xposarr[ind]
         ycen[idw] = yposarr[ind]
         fchi2[idw] = chi2arr[ind]
      endfor
      chi2sizearr[*,is]=fchi2
      xcen_arr[*,is] = xcen
      ycen_arr[*,is] = ycen
      totchi2arr[is] = total(fchi2)
   endfor

   tmp = min(totchi2arr,ind_s)
   xcen_final = xcen_arr[*,ind_s]
   ycen_final = ycen_arr[*,ind_s]
stop
;         id1 = where(chi2arr gt tmp+0.9 and chi2arr lt tmp+1.1,nid1)
;         xid1rec[0:(nid1<100)-1,idw]=xposarr[id1]
;         yid1rec[0:(nid1<100)-1,idw]=yposarr[id1]


   tmp=min(abs(centerwave-5300),refwave)
   print,observation.exposure,centerwave[refwave],xcen[refwave],ycen[refwave]
   status=mldar((observation.ha-(slitmap[bind[0]].ra-slitmap[bind[0]].cenra))/15.,plugmap[bind[0]].dec,centerwave,parangle,offsetx,offsety,waveref=centerwave[refwave])
   xrange=minmax([-xcen+xcen[refwave],offsetx])
   yrange=minmax([-ycen+ycen[refwave],offsety])
;   psplot,'dar_off.ps',/square,/color
   !x.margin=[10,3]
   plot,[-xcen+xcen[refwave],offsetx],[-ycen+ycen[refwave],offsety],ps=1,/nodata,xtitle='d_RA (arcsec)',ytitle='d_Dec (arcsec)',xcharsize=1.3,ycharsize=1.3,title=bundleid
   oplot,-xcen+xcen[refwave],-ycen+ycen[refwave],ps=-1,symsize=2
   oplot,offsetx,offsety,ps=-4,linest=2,symsize=2,color=249
   legend,['Measured DAR','Theoretical DAR'],linest=[0,2],psym=[1,4],symsize=[2,2],/bottom,box=0
   legend,[string(observation.exposure,format='(i8.8)'),bundleid],box=0,/right,/top
;stop
   for i=0,nwave-1 do begin &$
     t=where(xid1rec[*,i] ne 0.0) &$
     oplot,-xid1rec[t,i]+xcen[refwave],-yid1rec[t,i]+ycen[refwave],ps=1,symsize=0.1 &$
   endfor

   x1 = interpol(fiberoffx[reffiber]+xcen,centerwave,10^newloglam,/spline)
   y1 = interpol(fiberoffy[reffiber]+ycen,centerwave,10^newloglam,/spline)
   ; fraction of light going down the central fiber at each wavelength
   fraction = modelthrupt(x1,y1,10^newloglam)

   model = modflux*fraction
   mratio = sbfluxarr[*,reffiber]/model
   mrativar = sbivararr[*,reffiber]*model^2
   mrativar = mrativar *(1- spflux_masklines(newloglam,/stellar))
   bb=where(newloglam lt alog10(6000.))
   b_sset=spflux_bspline(newloglam[bb],mratio[bb],mrativar[bb],everyn=10)
   rr=where(newloglam gt alog10(6000.))
   r_sset=spflux_bspline(newloglam[rr],mratio[rr],mrativar[rr],everyn=1.5)
   mratfit = newloglam*0.0
   mratfit[bb]=bspline_valu(newloglam[bb],b_sset)
   mratfit[rr]=bspline_valu(newloglam[rr],r_sset)

end
