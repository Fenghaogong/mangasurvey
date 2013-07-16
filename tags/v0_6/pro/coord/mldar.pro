;+
; function mldar
;
; DAR calculations for MANGA simulations
; Simple function call returning offset and
; parallactic angle for specified HA, decl,
; and wavelength.  Recall that BLUE wavelengths
; look HIGHER in the sky.
;
; Simple calculations are based on the original code by
; Enrico Marchetti, ESO, January 2001
;
; Now additionally includes a more complex option with
; the /distort flag that incorporates the effects of
; chromatic distortions in the Sloan telescope optics.
; This option requires specifying the RACEN and DECCEN
; coordinates of the plate center, in addition to having
; the BOSS plate design code in your IDL include path.
;
; Required input numbers:
; HA: Hour angle of object in decimal hours
; decl: Target declination of object in decimal degrees
; wave: Wavelength of observation in Angstroms
;
; Optional inputs:
; waveREF: Reference wavelength in Angstroms (5500 Angstroms default)
; TC: Temperature (in degrees C)
; RH: Relative humidity (in percent)
; P: Atmospheric pressure (in millibars)
; RAOBJ: RA of the target in decimal degrees (only needed if /distort set)
; RACEN: Central RA of the plate in decimal degrees (only needed if /distort set)
; DECCEN: Central DEC of the plate in decimal degrees (only needed if /distort set)
; /distort: If set, include effects of chromatic distortion in Sloan optics
;   (NB: Requires RA, RACEN, DECCEN to be set also)
; /debug: If set, prints airmass of observation
;
; Returns magnitude of DAR in arcseconds
;
; Optional parameter returns:
; parangle: parallactic angle in decimal degrees
; offsetX: net DAR X offset (xfocal) in arcsec
; offsetY: net DAR Y offset (yfocal) in arcsec
; altitude: altitude of observation in decimal degrees
;  (apparent altitude, i.e. does not correct for refraction)
;
; Author:
;   David R. Law (Dunlap Institute; drlaw@di.utoronto.ca)
;   First written: 07/01/2012
;   Last modified: 03/25/2013
;
; Modification history:
;   v1: Implemented with commentary.  Hacked from central routines
;       originally coded by Enrico Marchetti, ESO, January 2001
;   v2: Made temperature, RH, pressure optional inputs.  Default
;       parameters are set to characteristic values for APO based on
;       2011 meteorological data available at http://weather.apo.nmsu.edu/logs/zip.html
;   v3: 02-Oct-2012 D. Law
;       Changed function name to mldar from mangalib_darcalc to
;       conform with new filename conventions.
;   v3.1: 06-Nov-2012 D. Law
;       Added optional return of altitude.  Inverted sign on offsetY.
;   v4.0: 25-Feb-2013 D. Law
;       Major addition of DAR calculated incorporating optical distortions
;       in the Sloan optics.  Use by optional parameter flag /distort,
;       and providing values for RACEN and DECCEN.
;-
function mldar,HA,decl,wave,parangle,offsetX,offsetY,altitude,waveREF=waveREF,TC=TC,RH=RH,P=P,RAOBJ=RAOBJ,RACEN=RACEN,DECCEN=DECCEN,DISTORT=DISTORT,debug=debug
lat = TEN(32, 46 ,49)/!RADEG  ; latitude of Apache Point Observatory, in radians

; Set default values
if (keyword_set(waveREF)) then waveREFm=waveREF/10000. $
else waveREFm=5500./10000.; Reference guiding wavelength (microns)
if (keyword_set(TC)) then TC=TC $
else TC=10.5D; Temperature (Celsius)
if (keyword_set(RH)) then RH=RH $
else RH=24.5D; Relative humidity (%)
if (keyword_set(P)) then P=P $
else P=730.0D; Pressure (mbar)
if (keyword_set(RAOBJ)) then RAOBJ=RAOBJ $
else RAOBJ=0.; Target RA (only matters for /distort)
if (keyword_set(RACEN)) then RACEN=RACEN $
else RACEN=0.; Plate central RA (only matters for /distort)
if (keyword_set(DECCEN)) then DECCEN=DECCEN $
else DECCEN=decl; Plate central declination (only matters for /distort)

; Define some output variables
DR=0.

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; If /distort specified, do calculations using BOSS plate design code
if (keyword_set(DISTORT)) then begin
  lst=HA*15.+RAOBJ; Recover LST of observation in decimal degrees
  if (lst gt 360.) then lst=lst-360.

  focalscale=16.5338 ; Plate scale in arcsec/mm

  ; Run ad2xyfocal once at the guide wavelength to get base plate position
  ad2xyfocal,RAOBJ,decl,xfocal0,yfocal0,racen=RACEN,deccen=DECCEN,airtemp=TC,pressure=P,lst=lst,lambda=waveREFm*10000.
  ; Run ad2xyfocal again at the wavelength of interest to get offset position
  ad2xyfocal,RAOBJ,decl,xfocal,yfocal,racen=RACEN,deccen=DECCEN,airtemp=TC,pressure=P,lst=lst,lambda=wave
  ; Difference results to get the offsetx in xfocal and yfocal in arcsec
  offsetx=(xfocal-xfocal0)*focalscale
  offsety=(yfocal-yfocal0)*focalscale
  ; Total differential refraction
  DR=sqrt(offsetx*offsetx+offsety*offsety)

  ; Calculate *effective* parallactic angle after optical distortions
  parangle=ATAN(offsety,offsetx)*!RADEG
  ; Reference it to 0=North, 90=East
  parangle=90.-parangle
  if (parangle gt 360.) then parangle=parangle-360.
  if (parangle lt 0.) then parangle=parangle+360.

  ; Calculate altitude for reference
  cos_z = SIN(lat)*SIN(decl/!RADEG) + COS(lat)*COS(decl/!RADEG)*COS(HA/24.*2.*!pi)
  z = ACOS(cos_z)*!RADEG ;convert to degrees
  altitude=90.-z
endif

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; If /distort not set, do calculations using standard ESO routines
if (not(keyword_set(DISTORT))) then begin
  declr=decl/!RADEG; convert declination to radians
  HAr=HA/24.*2.*!pi; convert hour angle to radians
  wavem=wave/10000.; convert wavelength to microns

  ; Zenith distance
  cos_z = SIN(lat)*SIN(declr) + COS(lat)*COS(declr)*COS(HAr)
  z = ACOS(cos_z)*!RADEG ;convert to degrees
  altitude=90.-z

  if (keyword_set(debug)) then print,'AM = ',1./cos_z

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; Calculate magnitude of DAR
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ZD=z*!DtoR ; zenith distance in radians
  T=TC+273.16D; temperature in Kelvin

  PS=-10474.0+116.43*T-0.43284*T^2+0.00053840*T^3
  P2=RH/100.0*PS
  P1=P-P2
  D1=P1/T*(1.0+P1*(57.90D*1.0E-8-(9.3250D*1.0E-4/T)+(0.25844D/T^2)))
  D2=P2/T*(1.0+P2*(1.0+3.7E-4*P2)*(-2.37321E-3+(2.23366/T)-(710.792/T^2)+(7.75141E4/T^3)))
  S0=1.0/waveREFm & S=1.0/wavem
  N0_1=1.0E-8*((2371.34+683939.7/(130-S0^2)+4547.3/(38.9-S0^2))*D1+$
     (6487.31+58.058*S0^2-0.71150*S0^4+0.08851*S0^6)*D2)
  N_1=1.0E-8*((2371.34+683939.7/(130-S^2)+4547.3/(38.9-S^2))*D1+$
     (6487.31+58.058*S^2-0.71150*S^4+0.08851*S^6)*D2)
  DR=Tan(ZD)*(N0_1-N_1)*206264.8

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; Calculate parallactic angle
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ; Calculate parallactic angles: Fillipenko 1982 formula only used inverse
  ; sin, therefore considerable ambiguity about whether one should use 
  ; value, 180-value, or -(180+value).  D.R. Law has updated this computation to
  ; use more complete inverse tangent calculation, remedying ambiguity.

  par_ycomp=sin(HAr)*cos(lat)*cos(declr)
  par_xcomp=sin(lat)-sin(declr)*cos_z
  parangle=ATAN(par_ycomp,par_xcomp)*!RADEG

  offsetX=-DR*sin(parangle/!RADEG)
  offsetY=-DR*cos(parangle/!RADEG)
endif

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

return,DR
end
