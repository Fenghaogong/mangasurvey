pro posoffset,plate,mjd

   platenum = string(plate,format='(i0.0)')
   path=strcompress(getenv('MANGACORE_DIR')+'/'+getenv('MANGACORE_VER')+'/drill/',/remove_all)
   drillfile=strcompress(path+'plDrillPos-'+platenum+'.par')
   drillcmmfile=strcompress(path+'plate'+platenum+'_CMM.par')

   slmread=mlreadslm(plate,mjd,slitmap)

   ifus=where(strmatch(slitmap.ifuname,'ma*'))
   uu = uniq(slitmap[ifus].ifuname,sort(slitmap[ifus].ifuname))
   bundleids=slitmap[ifus[uu]].ifuname
   

   focalscale = 16.5338; arcsec/mm
   
   nbundle=n_elements(bundleids)
   xoff = fltarr(nbundle)
   yoff = fltarr(nbundle)
   for i=0,nbundle-1 do begin
     bind=where(slitmap.ifuname eq bundleids[i])
     objra = slitmap[bind[0]].ra
     objdec=slitmap[bind[0]].dec
     status=mlreaddrillfiles(drillfile,drillcmmfile,objra,objdec,xdrill,ydrill)
     xoff[i] = xdrill/1000.*focalscale
     yoff[i] = ydrill/1000.*focalscale
   endfor
   for i=0,nbundle-1 do print,bundleids[i],xoff[i],yoff[i]
end
