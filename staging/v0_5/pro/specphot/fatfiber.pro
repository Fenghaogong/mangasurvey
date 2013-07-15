pro fatfiber, exposure

    if n_elements(exposure) gt 1 then begin
       for i=0,n_elements(exposure)-1 do fatfiber,exposure[i]
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


    ditherpos = [[-0.417,-0.721],[0.833,0.0],[-0.417,0.721]]

    obspos = [observation.dcnra,observation.dcndec]
    distsq = (ditherpos[0,*]-obspos[0])^2+(ditherpos[1,*]-obspos[1])^2
    mindist = min(distsq,ditherid)
    dithername = string(ditherid+1,format='(i0.0)')

    status=mlreadslm(observation.plate,observation.mjd,slitmap)

    fatid = where(strmatch(slitmap.ifuname, 'STD5D'+dithername+'*'),nfat)
    if nfat eq 0 then message,'No FAT fiber found centered for this position.'

    mlframe_read,bfile,fatid,objflux=bflux,objivar=bivar,mask=bmask,loglam=bloglam,dispimg=bdisp,wset=bwset,plugmap=plugmap
    mlframe_read,rfile,fatid,objflux=rflux,objivar=rivar,mask=rmask,loglam=rloglam,dispimg=rdisp,wset=rwset

    correct_dlam,bflux,bivar,bwset,dlam=dloglam
    correct_dlam,bdisp,0,bwset,dlam=dloglam,/inverse
    correct_dlam,rflux,rivar,rwset,dlam=dloglam
    correct_dlam,rdisp,0,rwset,dlam=dloglam,/inverse

;    dloglam=1.d-4
    maxlam = 10500.
    minlam = 3500.
    nlam = ceil((alog10(maxlam)-alog10(minlam))/dloglam)

    newloglam = alog10(minlam)+dloglam*dindgen(nlam)

    totfluxarr = fltarr(nlam,nfat)
    totivararr = totfluxarr*0.0
    totdisparr = totfluxarr*0.0

    for i=0,nfat-1 do begin
       combine1fiber,bloglam[*,i],bflux[*,i],bivar[*,i],finalmask=bmask[*,i],indisp=bdisp[*,i],newloglam=newloglam,newflux=bnewflux,newivar=bnewivar,andmask=bnewmask,newdisp=bnewdisp
       combine1fiber,rloglam[*,i],rflux[*,i],rivar[*,i],finalmask=rmask[*,i],indisp=rdisp[*,i],newloglam=newloglam,newflux=rnewflux,newivar=rnewivar,andmask=rnewmask,newdisp=rnewdisp
       totivararr[*,i] = bnewivar+rnewivar
       totfluxarr[*,i] = (bnewflux*bnewivar+rnewflux*rnewivar)/totivararr[*,i]
       totdisparr[*,i] = (bnewdisp*bnewivar+rnewdisp*rnewivar)/totivararr[*,i]
       bad = where(totivararr[*,i] eq 0, nbad)
       if nbad gt 0 then begin
          totfluxarr[bad,i] = (bnewflux[bad]+rnewflux[bad])/2.
	  totdisparr[bad,i] = (bnewdisp[bad]+rnewdisp[bad])/2.
       endif
    endfor

    skyid = where(strmatch(slitmap.ifuname,'SKY5*'),nfatsky) 
    if nfatsky eq 0 then message,'No FAT sky fibers.'
    mlframe_read,bfile,skyid,objflux=bskyflux,objivar=bskyivar,mask=bskymask,loglam=bskyloglam,dispimg=bskydisp,wset=bskywset
    mlframe_read,rfile,skyid,objflux=rskyflux,objivar=rskyivar,mask=rskymask,loglam=rskyloglam,dispimg=rskydisp,wset=rskywset

    correct_dlam,bskyflux,bskyivar,bskywset,dlam=dloglam
    correct_dlam,bskydisp,0,bskywset,dlam=dloglam,/inverse
    correct_dlam,rskyflux,rskyivar,rskywset,dlam=dloglam
    correct_dlam,rskydisp,0,rskywset,dlam=dloglam,/inverse

    sbfluxarr = totfluxarr*0.0
    sbivararr = sbfluxarr
    
    uu = uniq(slitmap[fatid].blockid,sort(slitmap[fatid].blockid))
    uublockid = slitmap[fatid[uu]].blockid
    for i=0,n_elements(uu)-1 do begin
       s_in_block = where(slitmap[skyid].blockid eq uublockid[i],nskyinblock)
       if nskyinblock eq 0 then message,'No fat sky fibers in this block.'
       combine1fiber,bskyloglam[*,s_in_block],bskyflux[*,s_in_block],bskyivar[*,s_in_block],finalmask=bskymask[*,s_in_block],newloglam=newloglam,newflux=bskynewflux,newivar=bskynewivar,andmask=bnewskymask
       combine1fiber,rskyloglam[*,s_in_block],rskyflux[*,s_in_block],rskyivar[*,s_in_block],finalmask=rskymask[*,s_in_block],newloglam=newloglam,newflux=rskynewflux,newivar=rskynewivar,andmask=rnewskymask
       skyivar = bskynewivar+rskynewivar
       skyflux = (bskynewflux*bskynewivar+rskynewflux*rskynewivar)/skyivar
       bad = where(skyivar eq 0,nbad)
       if nbad gt 0 then skyflux[bad] = (bskynewflux[bad]+rskynewflux[bad])/2.

       p = where(slitmap[fatid].blockid eq uublockid[i],nfatin)
       sbfluxarr[*,p] = totfluxarr[*,p]-skyflux#(fltarr(nfatin)+1)
       sbivararr[*,p] = 1/(1/totivararr[*,p]+1/(skyivar#(fltarr(nfatin)+1)))
       bad = where(totivararr[*,p] eq 0 or skyivar#(fltarr(nfatin)+1) eq 0,nbad)
       if nbad gt 0 then sbivararr[bad] = 0.0
    endfor
      
    totfit = fltarr(nlam,nfat)
    for i=0,nfat-1 do begin
       tmp = typingmodule(sbfluxarr[*,i],newloglam,sbivararr[*,i],totdisparr[*,i],i,observation,plugmap,kindx=kindx,modflux=modflux)
;       stop
       model = modflux
       mratio = sbfluxarr[*,i]/model
       mrativar = sbivararr[*,i]*model^2
       mrativar = mrativar*(1-spflux_masklines(newloglam,/stellar))
       bb = where(newloglam lt alog10(6000.))
       b_sset = spflux_bspline(newloglam[bb],mratio[bb],mrativar[bb],everyn=10)
       rr=where(newloglam gt alog10(6000.))
       r_sset=spflux_bspline(newloglam[rr],mratio[rr],mrativar[rr],everyn=1.5)
       totfit[bb,i]=bspline_valu(newloglam[bb],b_sset)
       totfit[rr,i]=bspline_valu(newloglam[rr],r_sset)
    endfor
    
    single = {loglam:newloglam,fiberid:-1,thrupt:fltarr(nlam)}
    cat = replicate(single,nfat)
    cat.fiberid = fatid+1
    cat.thrupt = totfit
    mwrfits,cat,'thrupt5'+string(exposure,format='(i8.8)')+'.fits',/create
end
