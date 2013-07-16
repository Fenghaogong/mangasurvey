;+
; Name: guider_psfcat
; Purpose: top level routine to collect psf, transparency, and dither info 
;           obtained from guider for each science exposure
;-
pro guider_psfcat

   dir=getenv('MANGA_SPECTRO_DATA')

   list = file_search(dir,'sdR-b1-*.fit*',count=ct)

   zerocat={mjd:0L,exposure:0L,plate:0,flatcart:0,dcnra:0.0,dcndec:0.0,params:dblarr(4),fwhm:0.0,transpar:0.0,guider1:' ',guidern:' ',SEXPTIME:0.0,GEXPTIME:0.0,alt:0.0,az:0.0}
   totcat = replicate(zerocat,ct)

   lat = 32+46/60.+49/3600.d
   lon = -(105.+49/60.+13/3600.d)
   j=0
   for i=0,ct-1 do begin
     hdr=headfits(list[i])
     flavor = sxpar(hdr,'FLAVOR')
     cart = sxpar(hdr,'CARTID')
     if strmatch(flavor,'science*',/fold) and cart eq 1 then begin
        single=zerocat
        mjd = sxpar(hdr,'MJD')
        expno=sxpar(hdr,'EXPOSURE')
        print,expno,' ',flavor
;        coaddgimg,mjd,expno
        cat=psfparams(mjd,expno)
        p = where(cat.mfwhm gt 0.8 and cat.fiberid le 16 and cat.focusoffset eq 0.0)
        ss = sort(cat[p].mfwhm)
        pick = p[ss[0]]
        single.fwhm =cat[pick].mfwhm
        struct_assign,cat[pick],single,/nozero
        single.exposure = expno
        single.gexptime = cat[pick].exptime
        sigmas = sqrt(cat[pick].gfit[2:3])*0.428 ; arcsec
        totflux = 2*!pi*(cat[pick].gfit[4]*sigmas[0,*]^2+cat[pick].gfit[5]*sigmas[1,*]^2)
        single.params[0:1]=cat[pick].gfit[4:5]/totflux[0]
        single.params[2:3]=sqrt(cat[pick].gfit[2:3])*0.428
        appendguidemag,cat
        transpar=10^(-0.4*(25.6666 - (cat.ugriz[1]+cat.ugriz[2])/2.))*(cat.measureflux/cat.exptime)
        single.transpar=median(transpar[p])
        ra = sxpar(hdr,'RADEG')
        dec = sxpar(hdr,'DECDEG')
        jd = (sxpar(hdr,'TAI-BEG')+sxpar(hdr,'EXPTIME')/2.)/86400.+2400000.5
        eq2hor,ra,dec,jd,alt,az,ha,lat=lat,lon=lon,altitude=2788.
        single.alt = alt
        single.az = az
        totcat[j++]=single
     endif
   endfor
   totcat=totcat[0:j-1]
   outputdir=getenv('GCAM_REDUX')
   mwrfits,totcat,outputdir+'gpsfcat.fits',/create
end
