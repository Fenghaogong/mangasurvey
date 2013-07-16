;+
; function mlmakeexp
;
; Generates a single .exp exposure information file for all science 
; frames in a given folder using header info.  Ordinarily would 
; do the whole thing automatically, but given plugging and mapping 
; issues do the bundle assignment by hand, so replicate lines from 
; sample exp map.  Stored with science frames.
;
; Author:
;   David R. Law (Dunlap Institute; drlaw@di.utoronto.ca)
;   First written: 12/30/2012
;   Last modified: 01/25/2013
;
; REVISION HISTORY:
;   v1: 30-DEC-2012  D. Law
;       First written.
;   v1.1: 23-JAN-2013  D. Law
;       Name format tweak to include leading '00' in exposure number.
;       Modification to run automatically
;   v1.2: 25-JAN-2013 D. Law
;       Modification to include taistart in output
;-

function mlmakeexp,inputpath,outputname,flavor

status=0L

; Dialog window to select the science frames
; and define directory.  Filter to include only b1
; sci and spFrame files.
if flavor eq 'sos' then frames=file_search(inputpath+'sci*b1*.fits')
if flavor eq 'full' then frames=file_search(inputpath+'spFrame*b1*.fits.gz')

; Check that some files were found.  If not, error -30
if (size(frames))[0] eq 0 then begin
  splog,'ERROR: No science data in directory',inputpath
  return,-30L
endif

nexp=(size(frames))[1]

; Read sample exp file
samplefile=strcompress(getenv('MANGACORE_DIR')+'/'+getenv('MANGACORE_VER')+'/templates/exp-56280.par',/remove_all)
; Read observing parameter block
sampleparam=yanny_readone(samplefile,'OPARAM',hdr=ehdr)

; Replicate sample to number of exposures
obsparam=replicate(sampleparam,nexp)

; Get basic info
head=headfits(frames[0])
plugname=fxpar(head,'NAME')

; Find the corresponding plugmap file from mangadb
plugfile=strcompress(getenv('MANGACORE_DIR')+'/'+getenv('MANGACORE_VER')+'/plugmaps/plPlugMapM-'+plugname+'.par',/remove_all)
temp=findfile(plugfile,count=ct)
if ct eq 0 then begin
  splog,'ERROR: Cannot find plugmap',plugfile
  return,-36L
endif

; Read plugM file
plugm=yanny_readone(plugfile,'PLUGMAPOBJ',hdr=phdr)

; Define central plate coordinates from plugmap
cenra=yanny_par(phdr,'raCen')
cendec=yanny_par(phdr,'decCen')
obsparam[*].cenra=cenra
obsparam[*].cendec=cendec

; Define drill HA from plugmap
hadrill=(yanny_par(phdr,'ha'))[0]
; Convert to hours:min:sec
mlhasex,hadrill/15.,hadrillsex
obsparam[*].drillha=hadrillsex

; Loop over exposures populating information
for i=0,nexp-1 do begin
  head=headfits(frames[i])

  ; Figure out exposure number
  expnum=fxpar(head,'EXPOSURE')
  obsparam[i].expnum=strcompress('00'+string(expnum),/remove_all)
  ; Exposure time in seconds
  exptime=fxpar(head,'EXPTIME')
  obsparam[i].exptime=exptime
  ; LST
  tai=fxpar(head,'TAI-BEG'); MDJ(TAI) in seconds at exposure start
  ; Note that taistart is recorded as a string in obsparam to ensure
  ; all significant digits are retained
  obsparam[i].taistart=strn(tai,format='(d15.2)')
  jd = 2400000.5D + tai / (24.D*3600.D)
  longitude = 360. - 105.820417d0
  ct2lst, lstdecimal, longitude, junk, jd
  mlhasex,lstdecimal,lst
  obsparam[i].LSTstart=lst
  ; Hour angle
  cenha=lstdecimal*15.-cenra
  mlhasex,cenha/15.,cenhasex
  obsparam[i].cenha=cenhasex
  ; Date
  dobs=fxpar(head,'DATE-OBS')
  date=(strsplit(dobs,'T',/extract))[0]
  obsparam[i].date=date
  ; UT start
  ut=(strsplit(dobs,'T',/extract))[1]
  obsparam[i].UTstart=ut
  ; MJD
  mjd=fxpar(head,'MJD')
  obsparam[i].mjd=strcompress(mjd,/remove_all)
  ; Plate
  plate=fxpar(head,'PLATEID')
  obsparam[i].plate=strcompress(plate,/remove_all)
  ; Cart
  cart=fxpar(head,'CARTID')
  if (cart lt 10) then cart=strcompress('0'+string(cart),/remove_all) $
  else cart=strcompress(string(cart),/remove_all)
  obsparam[i].cart=cart
  ; Airmass
  am=tai2airmass(cenra,cendec,tai=tai)
  obsparam[i].airmass=am
  ; Dither position.  Currently set to -1, not used.
  dither=-1
  obsparam[i].dposn=dither
  ; seeing
  seeing=1.6
  obsparam[i].seeing=seeing
  ; cloud cover
  ccpct=0.0
  obsparam[i].ccpct=ccpct
  ; temperature
  tcel=fxpar(head,'AIRTEMP')
  obsparam[i].tcel=tcel
  ; humidity
  rh = fxpar(head,'HUMIDITY')
  obsparam[i].rhumid=rh
  ; Pressure, convert from in Hg to mb
  pres=fxpar(head,'PRESSURE')*1013.25/29.9213
  obsparam[i].pres=pres
endfor

; Determine output filename
outname=strcompress(inputpath+outputname,/remove_all)
yanny_write,outname,[ptr_new(obsparam)],hdr=ehdr

return,status
end

