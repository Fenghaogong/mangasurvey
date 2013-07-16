; function mlfluxcor_placeholder
;
; Placeholder function for flux correction
;
; Returns 0 if everything ok, returns an error code if there was a problem.
;
; Author:
;   David R. Law (Dunlap Institute; drlaw@di.utoronto.ca)
;   First written: 01/23/2013
;   Last modified: 02/01/2013
;
; REVISION HISTORY:
;   v1: 23-Jan-2013  D. Law
;       Placeholder function written
;   v1.1: 01-Feb-2013  D. Law
;       Force replacement
;   v1.2: 05-Feb-2013  D. Law
;       Require providing output name

function mlfluxcor_placeholder, Cfile1, Ffile1


; Write out sky subtracted files
; In placeholder function, simply copy them to
; new names
;Ffile1=str_replace(Cfile1,'spCFrame','spFFrame')
spawn,strcompress('rm -f '+Ffile1)
spawn,strcompress('cp '+Cfile1+' '+Ffile1)

return,0
end
