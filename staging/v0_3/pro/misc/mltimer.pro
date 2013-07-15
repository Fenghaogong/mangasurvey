; program mltimer
;
; Timer for use in manga analysis code, very basic.
; Fails when the two input time differ in month or year.
;
; Author:
;   David R. Law (Dunlap Institute; drlaw@di.utoronto.ca)
;   First written: 10/11/2012
;   Last modified: 10/17/2012
;
; REVISION HISTORY:
;   v1: 17-Oct-2012  D. Law
;       First written
;
function mltimer,time1,time2

temp1=strsplit(time1,' ',/extract)
temp2=strsplit(time2,' ',/extract)

date1=temp1[2]+0
date2=temp2[2]+0

temp1a=temp1[3]
temp2a=temp2[3]

temp1b=strsplit(temp1a,':',/extract)
temp2b=strsplit(temp2a,':',/extract)

time1sec=temp1b[0]*(3600.D)+temp1b[1]*(60.D)+temp1b[2]
time2sec=temp2b[0]*(3600.D)+temp2b[1]*(60.D)+temp2b[2]

if (date1 ne date2) then time2sec=time2sec+(date2-date1)*(3600.D)*24.

dt=time2sec-time1sec
return,dt
end

