; DRLCosmo_Funcs.pro
; Implementation of some useful cosmological calculations
; for a variety of universe models.
; Many formulae based on Hogg, astro-ph/9905116
; NB: Seems to be working, but has not been thoroughly checked
; esp. at z>30
;
; David R. Law
; drlaw@astro.ucla.edu
; February 2010
; Last modified: February 15, 2010
;
; Standard cosmology:  h=0.732; OmM=0.238; OmL=0.762; OmK=0.;

; Given two redshifts, returns the absolute value of the velocity difference
; between them in the moving frame in km/s.
; For instance, given redshift of a source and redshift of an absorption line,
; calculates outflow velocity in the source frame.
function DRLCosmo_VOutflow,z1,z2
  c=2.99792458D5; // km/s
  v1=((1.+z1)*(1.+z1)-1.)/((1.+z1)*(1.+z1)+1.); // Currently in units of c
  v2=((1.+z2)*(1.+z2)-1.)/((1.+z2)*(1.+z2)+1.); // Currently in units of c
  if (v1 gt v2) then begin ;// Use convention v2>v1 to return positive output
    Temp=v1;
    v1=v2;
    v2=Temp;
  endif
  return, c*(v2-v1)/(1.-v1*v2);// Relativistic velocity addition formula
end

; Useful function
function DRLCosmo_EofZ,z,OmM,OmL,OmK
  return, sqrt(OmM*((1.+z)^3)+OmK*((1.+z)^2)+OmL)
end

; Hubble distance in Mpc
function DRLCosmo_DH,h
  c=2.9979D10; // cm/s
  H0=100.*h*1D5; // cm/s/Mpc
  return, c/H0;
end

; Line-of-sight comoving distance in Mpc
function DRLCosmo_DC,h,z,OmM,OmL,OmK
  zstep=1D-6;
  zprime=zstep;
  zintegral=0.;  

  while (zprime lt z) do begin
    zintegral+=zstep/DRLCosmo_EofZ(zprime,OmM,OmL,OmK);
    zprime+=zstep;
    if (zprime gt 1D-2) then zstep=1D-5;
    if (zprime gt 1.) then zstep=1D-3;
    if (zprime gt 30.) then zstep=zprime/10.;
  endwhile

  return, DRLCosmo_DH(h)*zintegral;
end

; Transverse comoving distance in Mpc
function DRLCosmo_DM,h,z,OmM,OmL,OmK
  epsilon=1D-6;
  dm=0.;

  ; If curvature is zero
  if (abs(OmK) lt epsilon) then dm=DRLCosmo_DC(h,z,OmM,OmL,OmK)

  ; If curvature is positive
  if (OmK ge epsilon) then dm=DRLCosmo_DH(h)*sinh(sqrt(OmK)*DRLCosmo_DC(h,z,OmM,OmL,OmK)/DRLCosmo_DH(h))/sqrt(OmK);

  ; If curvature is negative
  if (OmK lt -epsilon) then dm=DRLCosmo_DH(h)*sin(sqrt(-OmK)*DRLCosmo_DC(h,z,OmM,OmL,OmK)/DRLCosmo_DH(h))/sqrt(-OmK);

  return, dm;
end

; Angular diameter distance in Mpc
function DRLCosmo_DA,h,z,OmM,OmL,OmK
  return, DRLCosmo_DM(h,z,OmM,OmL,OmK)/(1.+z);
end

; Luminosity distance in Mpc
function DRLCosmo_DL,h,z,OmM,OmL,OmK
  return, (1.+z)*DRLCosmo_DM(h,z,OmM,OmL,OmK);
end

