function mlparang, HA, decl

lat = TEN(32, 46 ,49)/!RADEG  ; latitude of Apache Point Observatory, in radians

declr=decl/!RADEG; convert declination to radians
HAr=HA/24.*2.*!pi; convert hour angle to radians

; Zenith distance
cos_z = SIN(lat)*SIN(declr) + COS(lat)*COS(declr)*COS(HAr)

; Calculate parallactic angles: Fillipenko 1982 formula only used inverse
; sin, therefore considerable ambiguity about whether one should use 
; value, 180-value, or -(180+value).  D.R. Law has updated this computation to
; use more complete inverse tangent calculation, remedying ambiguity.

par_ycomp=sin(HAr)*cos(lat)*cos(declr)
par_xcomp=sin(lat)-sin(declr)*cos_z
parangle=ATAN(par_ycomp,par_xcomp)*!RADEG

return,parangle
end
