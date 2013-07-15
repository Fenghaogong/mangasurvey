pro specphotstat, obj

;   obj=mrdfits('fluxcal.fits',1)

   uu = uniq(obj.exposure,sort(obj.exposure))
   nexp = n_elements(uu)
   expos = obj[uu].exposure

   wave1=[3800,4950,6700]
   wave2=[3850,5000,6750]
   centerwave=(wave1+wave2)/2.

   lambda=10^obj[0].newloglam
   newloglam = obj[0].newloglam

   bb=where(lambda lt 6200.)
   rr=where(lambda gt 6200.)
   mratfit = obj.mratio*0.0
   for i=0,n_elements(obj)-1 do begin
      b_sset=spflux_bspline(newloglam[bb],obj[i].mratio[bb],obj[i].mrativar[bb],everyn=10)
      r_sset=spflux_bspline(newloglam[rr],obj[i].mratio[rr],obj[i].mrativar[rr],everyn=1.5)
;      mratfit = newloglam*0.0
      mratfit[bb,i]=bspline_valu(newloglam[bb],b_sset)
      mratfit[rr,i]=bspline_valu(newloglam[rr],r_sset)
   endfor

;   p1=where(10^obj[0].newloglam gt wave1[0] and 10^obj[0].newloglam lt wave2[0])
;   p2=where(10^obj[0].newloglam gt wave1[1] and 10^obj[0].newloglam lt wave2[1])
;   p3=where(10^obj[0].newloglam gt wave1[2] and 10^obj[0].newloglam lt wave2[2])
   snarray = obj.mratio*sqrt(obj.mrativar)

   for i=0,nexp-1 do begin
     d1 = where(obj.exposure eq expos[i],nbundle)

     meanthrupt = total(obj[d1].mratio,2)/nbundle
     diff = obj[d1].mratio-meanthrupt#(fltarr(nbundle)+1)
     rms = sqrt(total(diff*diff,2)/nbundle)
     err = rms/sqrt(nbundle-1)
     for j=0,2 do begin
        p = where(lambda gt wave1[j] and lambda lt wave2[j])
        sn = median((snarray[p,*])[*,d1],dimen=1)
;        print,sn
        good = where(sn gt 2,ngood)
        if ngood lt 4 then begin 
           good = indgen(nbundle)
           ngood = nbundle
           print,'force using all bundle.'
        endif
        rms = sqrt(total(diff[*,good]*diff[*,good],2)/(ngood-1))
        err = rms/sqrt(float(ngood))
        print,'Dither '+string(i,format='(i0.0)'),centerwave[j],median(err[p]/meanthrupt[p]), ' from ',ngood,' bundles'
     endfor
   endfor

   for i=0,nexp-1 do begin
     d1 = where(obj.exposure eq expos[i],nbundle)

     meanthrupt = total(mratfit[*,d1],2)/nbundle
     diff = mratfit[*,d1]-meanthrupt#(fltarr(nbundle)+1)
     rms = sqrt(total(diff*diff,2)/nbundle)
     err = rms/sqrt(nbundle-1)
     for j=0,2 do begin
        p = where(lambda gt wave1[j] and lambda lt wave2[j])
        sn = median((snarray[p,*])[*,d1],dimen=1)
;        print,sn
        good = where(sn gt 2,ngood)
        if ngood lt 4 then begin 
           good = indgen(nbundle)
           ngood = nbundle
           print,'force using all bundle.'
        endif
        rms = sqrt(total(diff[*,good]*diff[*,good],2)/ngood)
        err = rms/sqrt(ngood-1)
        print,'Dither '+string(i,format='(i0.0)'),centerwave[j],median(err[p]/meanthrupt[p]), ' from ',ngood,' bundles'
     endfor
   endfor

   for i=0,nexp-1 do begin
     flux = fltarr(3,nbundle)
     d1=where(obj.exposure eq expos[i],nbundle)
     for j=0,2 do begin
        p = where(lambda gt wave1[j] and lambda lt wave1[j]+20)
        flux[j,*] =  median(obj[d1].mratio[p], dimen=1)
     endfor
     ratio1 = flux[0,*]/flux[2,*]
     ratio2 = flux[1,*]/flux[2,*]
     print,stddev(ratio1)/mean(ratio1)/sqrt(nbundle)
     print,stddev(ratio2)/mean(ratio2)/sqrt(nbundle)
   endfor
end
   

