; function mlgetbundlesize
;
; Reference code stating how many fibers are in a bundle
; of the given 'ma00x' name.
;
; Author:
;   David R. Law (Dunlap Institute; drlaw@di.utoronto.ca)
;   First written: 01/17/2013
;   Last modified: 01/17/2013
;
; REVISION HISTORY:
;   v1: 17-Jan-2013  D. Law
;       First written.
function mlgetbundlesize,ifuname

; Default is negative, undefined bundle name
nfiber=-1

; Test-run bundles
if ifuname eq 'ma001' then nfiber=19
if ifuname eq 'ma002' then nfiber=61
if ifuname eq 'ma003' then nfiber=127
if ifuname eq 'ma004' then nfiber=19
if ifuname eq 'ma005' then nfiber=19
if ifuname eq 'ma006' then nfiber=19
if ifuname eq 'ma007' then nfiber=19
if ifuname eq 'ma008' then nfiber=127

return,nfiber
end
