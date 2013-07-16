; function mlmakeslitmap
;
; Generates a slitmap for a given plugM file and sticks it in the database
;
; Author:
;   David R. Law (Dunlap Institute; drlaw@di.utoronto.ca)
;   First written: 01/05/2013
;   Last modified: 01/16/2013
;
; REVISION HISTORY:
;   v1: 05-JAN-2013  D. Law
;       First written.
;   v1.1: 16-Jan-2013 D. Law
;       Revised to automate and set up for new directory system

pro mlmakeslitmap,plate,mjd

; Find plugM file from svn archive
plugMtarget=strcompress(getenv('MANGACORE_DIR')+'/'+getenv('MANGACORE_VER')+'/plugmaps/plPlugMapM-'+string(plate)+'-',/remove_all)
plugMfile=file_search(plugMtarget+'*')
plugMfile=mlmatchmjd(plugMfile,plugMtarget,mjd)
; Read plugmap
plugM=yanny_readone(plugMfile,'PLUGMAPOBJ',hdr=phdr)
plugM=plugM[where(plugM.holeType eq 'MANGA')]


; Find plateHoles file from svn archive
; Need to do something more clever here when plateid >9999
platefile=strcompress(getenv('MANGACORE_DIR')+'/'+getenv('MANGACORE_VER')+'/platedesign/plateholes/plateHolesSorted-00'+string(plate)+'.par',/remove_all)
platefile=file_search(platefile)
; If platefile not found, exit
if (size(platefile))[1] ne 1 then begin
    print,'ERROR: No platefile found in archive directory as'
    print,platefile
    exit
endif
splog,'Using ',platefile
; Read plateHoles
holes=yanny_readone(platefile,'STRUCT1')

; Define cart #, format to string with leading 0 for #<10
cart=yanny_par(phdr,'cartridgeId')
if cart le 9 then cart=strcompress('0'+string(cart),/remove_all) $
else cart=strcompress(cart,/remove_all)


; Define the slitmap file to produce in svn archive
slitfile=strcompress(getenv('MANGACORE_DIR')+'/'+getenv('MANGACORE_VER')+'/slitmaps/slitmap-'+string(plate)+'-'+string(mjd)+'.par',/remove_all)


; Figure out which slitmap template is appropriate
slittemplate=strcompress(getenv('MANGACORE_DIR')+'/'+getenv('MANGACORE_VER')+'/templates/slitmap-cart'+cart+'-',/remove_all)
slittemplatefile=file_search(slittemplate+'*')
slittemplatefile=mlmatchmjd(slittemplatefile,slittemplate,mjd)

; Read BLOCKS entries
blocks=yanny_readone(slittemplatefile,'BLOCKS',hdr=shdr)
; Read in SMAP entries
smap=yanny_readone(slittemplatefile,'SMAP')
; Ordinarily we might take one line and copy it as needed,
; but for test run number doesnt change, so just modify
; lines

; Define central plate coordinates from plugmap
cenra=yanny_par(phdr,'raCen')
cendec=yanny_par(phdr,'decCen')
smap[*].cenra=cenra
smap[*].cendec=cendec

; Do all blocks in one loop
for i=1,560 do begin
  indexS=where(smap.fiberid eq i)
  indexM=where(plugM.fiberid eq i)
  if (indexM eq -1) then begin; Dead fiber
    smap[indexS].plugstatus=0
    smap[indexS].ra=0.
    smap[indexS].dec=0.
  endif else begin ; Live fiber
    smap[indexS].ra=plugM[indexM].ra
    smap[indexS].dec=plugM[indexM].dec
    smap[indexS].mag=plugM[indexM].mag
    ftype=plugM[indexM].objType
    ; Determine details for calibration fibers from
    ; plateHolesSorted file
    if ((ftype eq 'SPECTROPHOTO_STD')or(ftype eq 'SKY')) then begin
      ; Figure out corresponding line in plateHolesSorted
      ; based on ra and dec
      holeindex=where((holes.target_ra eq plugM[indexM].ra)and(holes.target_dec eq plugM[indexM].dec))
      type=holes[holeindex].targettype
      if (type eq 'sky') then smap[indexS].ifuname = 'SKY2'
      if (type eq 'sky3') then smap[indexS].ifuname = 'SKY3'
      if (type eq 'sky5') then smap[indexS].ifuname = 'SKY5'
      if (type eq 'standard_d1') then smap[indexS].ifuname = 'STD2D1'
      if (type eq 'standard_d2') then smap[indexS].ifuname = 'STD2D2'
      if (type eq 'standard_d3') then smap[indexS].ifuname = 'STD2D3'
      if (type eq 'standard3_d1') then smap[indexS].ifuname = 'STD3D1'
      if (type eq 'standard3_d2') then smap[indexS].ifuname = 'STD3D2'
      if (type eq 'standard3_d3') then smap[indexS].ifuname = 'STD3D3'
      if (type eq 'standard5_d1') then smap[indexS].ifuname = 'STD5D1'
      if (type eq 'standard5_d2') then smap[indexS].ifuname = 'STD5D2'
      if (type eq 'standard5_d3') then smap[indexS].ifuname = 'STD5D3'
    endif
  endelse
endfor

; Add some header keywords
nlines=(size(shdr))[1]
newhdr=strarr(nlines+2)
newhdr[0]=shdr
newhdr[nlines]=strcompress('plate '+string(plate))
newhdr[nlines+1]=strcompress('mjd '+string(mjd))

yanny_write,slitfile,[ptr_new(blocks),ptr_new(smap)],hdr=newhdr

return
end
