; function mlcheckbm
;
; Purpose is to check a bundlemap (.bm) file to make sure that it's ok.
; I.e., it has all of the fibers it should, none are missed,
; double-identified, or incorrectly identified
;
; Author:
;   David R. Law (Dunlap Institute; drlaw@di.utoronto.ca)
;   First written: 11/15/2012
;   Last modified: 01/21/2013
;
; REVISION HISTORY:
;   v1: 15-Nov-2012  D. Law
;       First written.
;   v1.1: 21-Jan-2013 D. Law
;       Tweaked logging, modified for new metrology format
function mlcheckbm,ifuname,bmap,bhdr

; Status starts off OK
status=0

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Check that there are the correct number of fibers in bundle,
; Check that the cart, ifusize, and ifuid flags in file match filename
; If fails, error code -12
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
nactual=(size(bmap))[1]
ntarget=mlgetbundlesize(ifuname)

if ifuname ne yanny_par(bhdr,'ifuname') then begin
  splog,strcompress('ERROR: Bad IFU name in '+ifuname+' found '+yanny_par(bhdr,'ifuname'))
  status=-12L
endif

if ntarget ne long(yanny_par(bhdr,'nfiber')) then begin
  splog,strcompress('ERROR: Bad IFU size in '+ifuname+' expected '+string(ntarget)+' found '+yanny_par(bhdr,'nfiber'))
  status=-12L
endif

if nactual ne ntarget then begin
  splog,strcompress('ERROR: Expected '+string(ntarget)+' fibers in bundle '+ifuname+', found '+string(nactual))
  status=-12L
endif

if status ne 0 then return,status

nfiber=nactual

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Check that fnum values are unique and go from 1-N
; If fails, error code -13
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
for i=1,nfiber do begin
  temp=where(bmap.isp eq (i-1)); Check spiral numbering (0 indexed)
  if (size(temp))[1] ne 1 then begin
    splog,strcompress('Bad isp fiber numbering on fiber'+string(i)+', duplicated or missing')
    splog,''
    status=-13L
  endif
  temp=where(bmap.ise eq i); Check serpentine numbering (1 indexed)
  if (size(temp))[1] ne 1 then begin
    splog,strcompress('Bad ise fiber numbering on fiber'+string(i)+', duplicated or missing')
    splog,''
    status=-13L
  endif
endfor

if status ne 0 then return,status

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Check that gbu flags are in the range -1 --> 1
; If fails, error code -15
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
if max(bmap.gbu) gt 1. then begin
  splog,'GBU flag greater than 1 specified.'
  status=-15L
endif
if min(bmap.gbu) lt -1 then begin
  splog,'GBU flag less than -1 specified.'
  status=-15L
endif

if status ne 0 then return,status

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Check that xpmm, ypmm positions are at least 130 um away
; from each other (no superimposed fiber cores if gbu=1).  Check none further than
; 120% of bundle diameter.
; If fails, error code -16
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
maxrad=150.*sqrt(nfiber/3.14)/1000.
for i=0,nfiber-2 do begin
  for j=i+1,nfiber-1 do begin
    offset=sqrt((bmap[i].xpmm-bmap[j].xpmm)^2+(bmap[i].ypmm-bmap[j].ypmm)^2)
    ; If offset too small or too big, and neither fiber is flagged as 'bad',
    ; there's a problem
    if (((min(offset) lt 0.130)or(max(offset) gt 2*maxrad))and(bmap[i].gbu ne -1)and(bmap[j].gbu ne -1)) then begin
      splog,strcompress('Possible location error in fibers ise: '+string(bmap[i].ise)+' and '+string(bmap[j].ise))
      status=-16L
    endif
  endfor
endfor

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Plot it visually
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
mlplotbm,bmap,bhdr

return,status
end
