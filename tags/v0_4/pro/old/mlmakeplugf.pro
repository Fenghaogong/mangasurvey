; function mlmakeplugf.pro
;
; Uses plPlugMapP.par together with plateInput file to generate
; plPlugMapF.par, which lists the positions of individual fibers
; within bundles for the automapper.  Positions are approximate
; and assume perfect bundle construction.
;
; Author:
;   David R. Law (Dunlap Institute; drlaw@di.utoronto.ca)
;   First written: 12/22/2012
;   Last modified: 01/05/2013
;
; REVISION HISTORY:
;   v1: 22-Dec-2012  D. Law
;       First written.

pro mlmakeplugf,plugmapP,plugmapF=plugmapF,plateinput=plateinput

; Original plugmap
plugmap=yanny_readone(plugmapP,'PLUGMAPOBJ',hdr=hdr)

; Figure out design id
designid=yanny_par(hdr,'designid')

; Figure out corresponding plateInput file based on design id
; (assume in same folder) if not specified
temp=stregex(plugmapP,'plPlugMapP')
path=strmid(plugmapP,0,temp)
if (keyword_set(plateinput)) then plateinput=plateinput $
else plateinput=strcompress(path+'plateInput-'+string(designid)+'.par')

; plateinput defines which IFU is in which hole
plateinput=yanny_readone(plateinput,'STRUCT1')

; Define plate center and drill hour angle
racen=yanny_par(hdr,'raCen')
deccen=yanny_par(hdr,'decCen')
ha=(yanny_par(hdr,'ha'))[0]

; Figure out which entries in plugP are IFUs
galaxy=plugmap[where(plugmap.objType eq 'GALAXY')]
ngal=(size(galaxy))[1]

; Loop through IFUs
for i=0,ngal-1 do begin
  ; For each IFU, determine hole ra, dec and the bundle_id
  ; in plateinput format of what IFU is in it.
  inputra=plateinput[i].ra
  inputdec=plateinput[i].dec
  bundleid=plateinput[i].bundle_id
  inputtype=strcompress(string(bundleid),/remove_all)
  len=strlen(inputtype)
  ; Translate bundle_id into the ifusize (e.g., 127) and
  ; ifu id (e.g., 1, for 127_1)
  ifusize=long(strmid(inputtype,0,len-1))
  ifuid=long(strmid(inputtype,len-1,1))

  basera=plateinput[i].ra
  basedec=plateinput[i].dec
  ; Figure out what the matching entry in plugP is and copy
  ; its plugmap entry
  origentry=galaxy[where(abs(galaxy.ra - basera) lt 0.0001)]
  ; Add the 'bundle_id' tag
  origentry=create_struct(origentry,'bundle_id',bundleid)
  ; Replicate to the number of fibers in bundle
  newentry=replicate(origentry,ifusize)

  ; Locations are approximate, so just use cart01 ifus for positions
  bmapfile=strcompress(getenv('MANGAROOT')+'/manga3d/'+getenv('MANGAVER')+'/database/bundlemaps/cart01_'+string(ifusize)+'_'+string(ifuid)+'.bm', /remove_all)
  bmap=yanny_readone(bmapfile,'BUNDLEMAP',hdr=bhdr)

  ; For each fiber within IFU, tweak ra, dec, xfocal, yfocal
  ; according to offset of fiber from bundle center
  for j=0,ifusize-1 do begin
    newentry[j].ra+= -bmap[j].xrel/60./3600./cos(newentry[j].dec*!PI/180.)
    newentry[j].dec+=bmap[j].yrel/60./3600.
    ad2xyfocal,newentry[j].ra,newentry[j].dec,xf,yf,racen=racen,deccen=deccen,lst=(ha+newentry[j].ra),lambda=5500.
    newentry[j].xFocal=xf
    newentry[j].yFocal=yf
  endfor

  ; Collate new information into plugfentries
  if (i eq 0) then plugfentries=newentry $
  else plugfentries=struct_append(plugfentries,newentry)
endfor

; Determine output file name
if (keyword_set(plugmapF)) then plugmapF=plugmapF $
else plugmapF=mlstrreplace(plugmapP,'plPlugMapP','plPlugMapF')

; Write plugF plugmap
yanny_write,plugmapF,ptr_new(plugfentries)

end
