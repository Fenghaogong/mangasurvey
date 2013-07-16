;+
; program mlrandplateloc
;
; Come up with nrand random locations on a plate
; This is used as input to optimization code that
; figures out the scale and rotation applied by the guider
; to best track objects on the plate.
;
; Required input:
;   racen: RA of plate center in decimal degrees
;   deccen: DEC of plate center in decimal degrees
;   nrand: Number of points to generate
;
; Optional input:
;   rngseed: Starting seed for random number generator
;
; Returned parameters:
;   rarand: RA values of random locations in decimal degrees
;   decrand: DEC values of random locations in decimal degrees
;
; Author:
;   David R. Law (Dunlap Institute; drlaw@di.utoronto.ca)
;   First written: 05/29/2013
;   Last modified: 06/10/2013
;
; REVISION HISTORY:
;   v1: 29-May-2013  D. Law
;       First written.
;   v1.1: 10-June-2013 D. Law
;       Added option to use fixed RNG seed.
;-

pro mlrandplateloc,racen,deccen,nrand,rarand,decrand,rngseed=rngseed

platediam=3.0; degrees

; If keyword rngseed was set, use that as seed to RNG
if (keyword_set(rngseed)) then begin
  ; Invert uniform probability distribution in area to get
  ; proper distribution in radius
  rrand=sqrt(randomn(rngseed,nrand,/uniform))*platediam/2.
  ; Distribution is uniform in angle
  thetarand=randomn(rngseed+1,nrand,/uniform)*2*!PI

; Otherwise use systime_seed
endif else begin
  ; Invert uniform probability distribution in area to get
  ; proper distribution in radius
  rrand=sqrt(randomn(systime_seed,nrand,/uniform))*platediam/2.
  ; Distribution is uniform in angle
  thetarand=randomn(systime_seed+1,nrand,/uniform)*2*!PI
endelse

rarand=rrand*cos(thetarand)/cos(deccen*!PI/180.)+racen
decrand=rrand*sin(thetarand)+deccen

return
end
