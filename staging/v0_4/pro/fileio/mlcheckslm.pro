; function mlcheckslm
;
; Purpose is to check a slitmap file to make sure that it's ok.
; I.e., it has all of the fibers it should, none are missed,
; double-identified, or incorrectly identified
;
; Author:
;   David R. Law (Dunlap Institute; drlaw@di.utoronto.ca)
;   First written: 11/14/2012
;   Last modified: 01/21/2013
;
; REVISION HISTORY:
;   v1: 14-Nov-2012  D. Law
;       First written.
;   v1.1: 15-Jan-2013 D. Law
;       Tweaked logging, changed format
function mlcheckslm,slitmap,blocks,shdr

; Status starts off OK
status=0

; Number of blocks
nblocks=(size(blocks))[1]
; Number of fibers
nfiber=fix(total(blocks.blocks_blocksize))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Check that there are the correct number of fibers in each block
; If fails, error code -23
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
for i=0,nblocks-1 do begin
  block=blocks[i].blocks_blockid
  ninblock=blocks[i].blocks_blocksize

  temp=slitmap[where(slitmap.blockid eq block)]
  if (size(temp))[1] ne ninblock then begin
    splog,strcompress('Bad block '+ string(block) +' length')
    splog,strcompress('Should be '+ string(ninblock) +', is '+string((size(temp))[1]))
    splog,''
    status=-23L
  endif
endfor

if status ne 0 then $
  return,status

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Check that each fiber number along slit is unique and reasonable
; (i.e., starts at 1, goes to nfiber)
; If fails, error code -24
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
for i=0,nfiber-1 do begin
  temp=slitmap[i].fiberid
  if (temp ne i+1) then begin
    splog,strcompress('Fiber number '+string(i+1)+' out of order- should be '+string(i+1)+' is '+string(temp))
    splog,''
    status=-24L
  endif

  temp=where(slitmap.fiberid eq i+1)
  if (size(temp))[1] ne 1 then begin
    splog,strcompress('Bad fiber numbering on fiber'+string(i)+', duplicated or missing')
    splog,''
    status=-24L
  endif
endfor

if status ne 0 then $
  return,status

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; IFUs can be identified by the 'ma' in ifuname.
; Figure out the number of unique identifiers with 'ma' in their
; ifuname; this is the number of unique IFUs.
; If fails, error code -25
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

ifunames=strarr(nfiber); Starts off nfiber in size, will resize later
nifu=0; Will contain number of IFUs
for i=0,nfiber-1 do begin
  if (strpos(slitmap[i].ifuname,'ma') ne -1) then temp=slitmap[i].ifuname; It's an IFU
  if where(ifunames eq temp) eq -1 then begin; Name wasn't already in the ifuname array
    ifunames[nifu]=temp
    nifu+=1
  endif
endfor
ifunames=ifunames[0:nifu-1]; Resize array

; Figure out how many fibers should be in each IFU
ifusiz=intarr(nifu)
for i=0,nifu-1 do ifusiz[i]=mlgetbundlesize(ifunames[i])

; Check that each fiber number for each bundle is unique and
; reasonable (i.e., in the range 1-ifusiz, no duplicates, no missing fibers)
for i=0,nifu-1 do begin
  temp=slitmap[where(slitmap.ifuname eq ifunames[i])].fnum; Vector of fnum for the bundle

  if (max(temp) ne ifusiz[i]) then begin ; Check maximum fiber number
    splog,strcompress('Maximum fiber number in bundle '+string(ifunames[i])+' is '+string(max(temp))+', should be '+string(ifusiz[i]))
    splog,''
    status=-25L
  endif

  if (min(temp) ne 1) then begin ; Check minimum fiber number
    splog,strcompress('Minimum fiber number in bundle '+string(ifunames[i])+' is '+string(min(temp))+', should be 1')
    splog,''
    status=-25L
  endif

  if (size(temp))[1] ne ifusiz[i] then begin ; Check total number of fibers
    splog,strcompress('Total number of fibers in bundle '+string(ifunames[i])+' is '+string((size(temp))[1])+', should be '+string(ifusiz[i]))
    splog,''
    status=-25L
  endif

  for j=1,ifusiz[i] do begin ; Check for duplicate fiber numbers
    if (size(where(temp eq j)))[1] ne 1 then begin
      splog,strcompress('Bad fnum in bundle '+string(ifunames[i])+', fiber '+string(j)+' is duplicated or missing.')
      splog,''
      status=-25L
    endif
  endfor
endfor

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Check plugstatus flags valid.  If one invalid, error code -26
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

temp=where((slitmap.plugstatus ne 0)and(slitmap.plugstatus ne 1))
if (size(temp))[0] ne 0 then begin
  splog,'Bad plugstatus value in slitmap'
  status=-26L
endif

return,status
end
