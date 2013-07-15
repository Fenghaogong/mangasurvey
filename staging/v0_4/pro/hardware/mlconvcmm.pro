; Program mlconvcmm
;
; Reads information from UWash CMM files and convert to a Yanny-style
; format for where the holes in a plate were actually drilled.
;
; Author:
;   David R. Law (Dunlap Institute; drlaw@di.utoronto.ca)
;   First written: 12/07/2012
;   Last modified: 12/10/2012
;
; REVISION HISTORY:
;   v1: 07-Dec-2012  D. Law
;       First written.
;
pro mlconvcmm, cmmfile

; Determine output file from input filename
temp=strsplit(cmmfile,'_',/extract)
plateno=temp[1]
outfile='plate'+plateno+'_CMM.par'

; Read header parameters first
openr,drillun,cmmfile,/get_lun

openw,lun,outfile,/get_lun

temp=''
readf,drillun,temp
printf,lun,temp
readf,drillun,temp
printf,lun,temp
readf,drillun,temp

readf,drillun,temp
printf,lun,'FitOffsetX '+(strsplit(temp,' ',/extract))[1]
printf,lun,'FitOffsetY '+(strsplit(temp,' ',/extract))[2]

for i=0,7 do begin
  readf,drillun,temp
  printf,lun,temp
endfor
readf,drillun,temp
readf,drillun,temp

; Read data vectors
readcol,cmmfile,nomx,nomy,measx,measy,residx,residy,residr,nomdia,diaerr,qpresidx,qpresidy,qpresidr,skipline=14

; Setup Yanny file
printf,lun,''
printf,lun,'typedef struct {'
printf,lun,'double nomx;'
printf,lun,'double nomy;'
printf,lun,'double measx;'
printf,lun,'double measy;'
printf,lun,'double residx;'
printf,lun,'double residy;'
printf,lun,'double residr;'
printf,lun,'double nomdia;'
printf,lun,'double diaerr;'
printf,lun,'double qpresidx;'
printf,lun,'double qpresidy;'
printf,lun,'double qpresidr;'
printf,lun,'} CMMPOS;'
printf,lun,''

print,size(nomx)
for i=0,(size(nomx))[1]-1 do begin
  printf,lun,strcompress('CMMPOS '+string(nomx[i])+string(nomy[i])+string(measx[i])+string(measy[i])+string(residx[i])+string(residy[i])+string(residr[i])+string(nomdia[i])+string(diaerr[i])+string(qpresidx[i])+string(qpresidy[i])+string(qpresidr[i]))
endfor

close,lun
close,drillun
free_lun,lun
free_lun,drillun

return
end
