pro snrstar, exposure, bundleid,newloglam,mratfit
;common com_fluxbundle, cube,xarr,yarr,wavearr
common com_fluxbundle2, psfplane,xarr,wavearr

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

    dloglam=1.d-4
    maxlam = 10600.
    minlam = 3400.
    nlam = ceil((alog10(maxlam)-alog10(minlam))/dloglam)

    newloglam = alog10(minlam)+dloglam*dindgen(nlam)

    bundleid = ['ma001','ma002','ma003','ma004','ma005','ma006','ma007','ma008']

    snr = fltarr(n_elements(bundleid))
    mags=fltarr(5,n_elements(bundleid))
    for ib=0,n_elements(bundleid)-1 do begin
      bind = where(slitmap.ifuname eq bundleid[ib],nfiber)
      ss = sort(slitmap[bind].fnum)
      bind = bind[ss]

      status=mlreadbm(bundleid[ib],observation.mjd,bmap)
; check bmap match with slitmap[bind] in fiberid
    if min(slitmap[bind].fnum-bmap.ise) ne 0 or max(slitmap[bind].fnum-bmap.ise) ne 0 then message,'Error: bmap and slitmap[bind] do not have matching fnum! some sorting is needed.'
      mags[*,ib] = plugmap[bind[0]].mag



    wave1 = [3500.,4000.,4500.,5000.,5500.,6500.,7500.,9000.]
    wave2 = [4000.,4500.,5000.,5500.,6500.,7500.,9000.,10500.]
    nwave = n_elements(wave1)
    centerwave = (wave1+wave2)/2.
    dwave = wave2-wave1
    skylevel = fltarr(nwave)
    skyerr = fltarr(nwave)

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
;    innermask = (sqrt(fiberx*fiberx+fibery*fibery) lt 3) and bmap.gbu ge 0

    fiberx0 = fiberx[ss[nfiber-1]]
    fibery0 = fibery[ss[nfiber-1]]
    innermask = (sqrt((fiberx-fiberx0)^2+(fibery-fibery0)^2) lt 3) and bmap.gbu ge 0
;    innermask = intarr(nfiber)
;    innermask[ss[nfiber-5:nfiber-1]]=1

    fnum = slitmap[bind].fnum
;innermask = (fnum ge 84 and fnum le 85) or (fnum ge 101 and fnum le 105)

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

    ind1 = where(10^newloglam gt 3700 and 10^newloglam lt 3750)
    snr[ib]=median(sbfluxarr[ind1,reffiber]*sqrt(sbivararr[ind1,reffiber]))/sqrt(3725*(10^(1.e-4)-1)) ; S/N per Angstrom
   endfor
   stop
end
