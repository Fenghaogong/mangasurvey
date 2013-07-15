;; ----------------------------------------------------------------------------------------------------------------------------------
;; NAME: MGEFLUX
;;
;; PURPOSE: This function computes 2D flux integrals in specified (circular) apertures for Multi-Gaussian Expansions.  It will
;; also integrate 2D SDSS (smoothed and truncated) versions of Exponential and de Vaucouleurs profiles, as well as standard Sersic
;; profiles with any value of Sersic-n.
;;
;; CALLING SEQUENCE: flux = mgeflux(radius, amp_mge [or Ie], sigma_mge [or Re], sigma_y_mge, XC=, YC=, /LUM_IS_AMP, /PRECISE, $
;;                                  /EXP_SDSS, /DEV_SDSS, sersic_n=, QPROF=1.0)
;;
;; ARGUMENTS:
;;   Radius:   Vector of radii specifying the apertures over which the flux will be integrated.  A value of -1 for an MGE
;;             integrates to infinity.
;;   amp_mge:  Either a vector of MGE amplitudes or, for a defined profile, the surface brightness at Re (Ie).  MGE amplitudes
;;             are defined such that the total 2D luminosity is given by: L2D = 2*!pi*sigma^2*amp_mge
;;   sigma_mge:  Either a vector of corresponding MGE sigma values or, for a defined profile, the effective radius, Re.
;;   sigma_y_mge:  Optional vector of MGE sigma values along the minor-axis.
;;   
;; KEYWORDS:
;;   /LUM_IS_AMP:  Set this keyword to specify that amp_mge values are actually 1D luminosities: L1D = amp_mge * sqrt(2*!pi*sigma^2)
;;   /PRECISE:     Set this keyword to used increased precision in the MGE numerical integration (but obviously, slower)
;;   /EXP_SDSS:    Integrate the SDSS truncated Exponential profile
;;   /DEV_SDSS:    Integrate the SDSS truncated de Vaucouleurs profile
;;   SERSIC_N:     Integrate a Sersic profile with this value of Sersic-n 
;;   XC:           Offset in x-direction between each (or all) gaussians and the center of the aperture.  If a scalar, then the
;;                 offset is applied to all terms.  If a vector, should have same number of elements as MGE terms.
;;   YC:           Same as XC for the y-direction.
;;   QPROF:        The axis ratio of the EXP_SDSS, DEV_SDSS, or Sersic profiles that will be integrated
;;
;; OUTPUT:
;;   Function returns the integrated flux at every radius specified.
;;
;; HISTORY:
;;   Created by K. Bundy (6/25/12) - Version 1.0
;;   Version 2.0, K. Bundy (8/7/12) - Generalized to off-centered apertures for MGE flux integral
;;
;;
;;  x' = x cos(theta) - y sin(theta)
;;  y' = x sin(theta) + y cos(theta)

;; For 2D Gaussian integral in Cartesian coordinates: Define 1D integrand
function gauss2D_1Dint, x
  common gauss2D_common, A0, sigma, sigma_y, qgauss, Rlim, xcen, ycen

  n_sum = n_elements(A0)

  if n_elements(xcen) + n_elements(ycen) EQ 0 then begin  ;; The aperture is centered at (0,0)

     if n_sum EQ 1 then begin  ;; Single component, and maybe we'll want to use the axis ratio (qgauss)
        if qgauss GT 0 then integrand = 2* A0 * sigma * qgauss * sqrt(!pi/2.0) * exp( -x^2/(2.0*sigma^2)) * erf(sqrt( (Rlim^2-x^2)/(2.0*sigma^2*qgauss^2) )) $
        else integrand = 2*A0 * sigma_y * sqrt(!pi/2.0) * exp( -x^2/(2.0*sigma^2)) * erf(sqrt( (Rlim^2-x^2)/(2.0*sigma_y^2) ))
     endif
     
     if n_sum GT 1 then begin  ;; MGE     
        integrand = total(2*A0 * sigma_y * sqrt(!pi/2.0) * exp( -x^2/(2.0*sigma^2)) * erf(sqrt( (Rlim^2-x^2)/(2.0*sigma_y^2) )))
     endif

  endif else begin  ;; The aperture is offset from center of the MGE

     ;; In case only one of xcen or ycen is defined, but not the other
     if n_elements(xcen) EQ 0 then xcen = fltarr(n_sum)
     if n_elements(ycen) EQ 0 then ycen = fltarr(n_sum)

     if n_sum EQ 1 then begin  ;; Single component, and maybe we'll want to use the axis ratio (qgauss)
        if qgauss GT 0 then integrand = A0 * sigma * qgauss * sqrt(!pi/2.0) * exp( -(x + xcen)^2/(2.0*sigma^2)) * $
                                        (erf( (ycen + sqrt(Rlim^2-x^2))/sqrt(2.0*sigma^2*qgauss^2)) - erf( (ycen - sqrt(Rlim^2-x^2))/sqrt(2.0*sigma^2*qgauss^2))) $
        else integrand = A0 * sigma_y * sqrt(!pi/2.0) * exp( -(x + xcen)^2/(2.0*sigma^2)) * $
                         (erf( (ycen + sqrt(Rlim^2-x^2))/sqrt(2.0*sigma_y^2)) - erf((ycen - sqrt(Rlim^2-x^2))/sqrt(2.0*sigma_y^2)))
     endif
     
     if n_sum GT 1 then begin  ;; MGE     
        integrand = total(A0 * sigma_y * sqrt(!pi/2.0) * exp( -(x + xcen)^2/(2.0*sigma^2)) * $
                          (erf((ycen + sqrt(Rlim^2-x^2))/sqrt(2.0*sigma_y^2)) - erf((ycen - sqrt(Rlim^2-x^2))/sqrt(2.0*sigma_y^2))))
     endif

  endelse

  ;stop
  return, integrand
end


;; For 2D Sersic integral: function returns angle integration limits
function thetalims, r
  return, [0, 2*!pi]
end


;; For 2D Gauss integral in polar coordinates: Define integrand
function gauss2D, r, theta

  common gauss2D_common, A0, sigma, sigma_y, qgauss, Rlim, xcen, ycen

integrand = A0 * r * exp(-(r^2/(2*sigma^2))*(1 + (1.0/qgauss^2 - 1)*(sin(theta))^2))

return, integrand
end


;; For 2D Sersic flux integral: Define integrand
function sersic2D, r, theta
  common sersic2D_common, I0, Re, n, bn, q

  ;integrand = I0 * r*exp(-bn * (r/Re)^(1.0/n)*(1 + (1.0/q^2 - 1)*(sin(theta))^2)^(1./(2.0*n)))
  integrand = I0 * r*exp(-bn * ((r/Re)^(1.0/n)*(1 + (1.0/q^2 - 1)*(sin(theta))^2)^(1./(2.0*n)) - 1))  ;; Makes I0 in Ie (surface brightness at Re)

  return, integrand
end

;; Implicit equation for bn
function sersic_bn, bn_arr
common sersic2D_common, I0, Re, n, bn, q
  return, gamma(2.0*n) - 2*igamma(2.0*n, bn_arr)*gamma(2.0*n)
end


;; SDSS: Smoothly truncated exponential profile, 2D polar coordinates 
;;   R: Radius, physical units (same units as those defining Re)
;;   theta: Polar angle, radians
;;   
function fexpcut, r, theta
  common sdss_profs, Ie_Exp, q_Exp, Re_Exp, Ie_Dev, q_Dev, Re_Dev  ;; Ie is surface brightness at Re
  if n_elements(q_exp) EQ 0 then q_Exp = 1.0

  q_Exp = float(q_Exp)
  rre = r / float(Re_exp)
  rvec = rre * sqrt((1 + (1.0/q_Exp^2 - 1)*(sin(theta)^2)))  ;; Positional vector in polar coordinates, assuming y-direction inclination axis ratio, q

  EXPFAC = -1.67835D  ;; Defines Ie_Exp and Re

  fexp = Ie_Exp * exp(ExpFac * (rvec-1.0))

  ExpOut = 4.0               ;; Profile goes to zero at this radius (units of Re)
  ExpCut = 3.0*ExpOut/4.0    ;; Profile truncation begins at this radius (units of Re)
  
  wr_cutoff = where(rre GT ExpCut, nr_cutoff)
  wr_zero = where(rre GT ExpOut, nr_zero)

  fexpcut = fexp
  if nr_cutoff GT 0 then fexpcut[wr_cutoff] = fexpcut[wr_cutoff] *  ( 1 - ((rvec[wr_cutoff] - ExpCut)/(ExpOut - ExpCut))^2)^2.0
  ;;if nr_cutoff GT 0 then fexpcut[wr_cutoff] = fexpcut[wr_cutoff] *  ( 1 - ((rvec[wr_cutoff] - ExpCut)/(ExpOut - ExpCut))^2)^1.0  ;; NOTE: To reproduce phFitobj.h, last exponent must be 1.0 (not 2.0)
  if nr_zero GT 0 then fexpcut[wr_zero] = 0

  return, fexpcut
end

;; Define integrand for the flux integral of the SDSS smoothly truncated exponential profile
function fexpcut_radintgrnd, r, theta
  return, r * fexpcut(r, theta)
end



;; SDSS deV: Smoothly truncated deVaucouleurs profile, 2D polar coordinates 
;;   R: Radius, physical units (same units as those defining Re)
;;   theta: Polar angle, radians
;;   
function fDeVcut, r, theta
  common sdss_profs, Ie_Exp, q_Exp, Re_Exp, Ie_Dev, q_Dev, Re_Dev  ;; Ie is surface brightness at Re
  if n_elements(q_Dev) EQ 0 then q_Dev = 1.0

  q_Dev = float(q_Dev)
  rre = r / float(Re_Dev)
  rvec = rre * sqrt((1 + (1.0/q_Dev^2 - 1)*(sin(theta)^2)))  ;; Positional vector in polar coordinates, assuming y-direction inclination axis ratio, q

  DEFAC = -7.66925D  ;; Defines Ie_Dev and Re

  ;fdeV = Ie_DeV * exp(DEFAC *(rvec^0.25 - 1.0))
  fdeV = Ie_DeV * exp(DEFAC *((rvec^2 + 0.0004)^0.125 - 1.0))

  deVOut = 8.0D
  deVCut = 7*DEVOUT/8.0

  wr_cutoff = where(rre GT deVCut, nr_cutoff)
  wr_zero = where(rre GT deVOut, nr_zero)

  fdeVcut = fdeV
  if nr_cutoff GT 0 then fdeVcut[wr_cutoff] = fdeVcut[wr_cutoff] *  ( 1 - ((rvec[wr_cutoff] - deVCut)/(deVOut - deVCut))^2)^2.0
  ;;if nr_cutoff GT 0 then fdeVcut[wr_cutoff] = fdeVcut[wr_cutoff] *  ( 1 - ((rvec[wr_cutoff] - deVCut)/(deVOut - deVCut))^2)^1.0 ;; NOTE: To reproduce phFitobj.h, last exponent must be 1.0 (not 2.0)
  if nr_zero GT 0 then fdeVcut[wr_zero] = 0

  return, fdeVcut
end

;; Define integrand for the flux integral of the SDSS smoothly truncated exponential profile
function fdeVcut_radintgrnd, r, theta
  return, r * fdeVcut(r, theta)
end


;; ----------------------------------------------------------------------------------------------------------------------------------
;; MAIN PROGRAM
;; ----------------------------------------------------------------------------------------------------------------------------------
function mgeflux, radius, amp, sigma_mge, sigma_y_mge, XC = XC, YC = YC, LUM_IS_AMP = LUM_IS_AMP, PRECISE = PRECISE, $
                  exp_sdss = exp_sdss, dev_sdss=dev_sdss, sersic_n=sersic_n, QPROF = QPROF

;; Error checking
if n_params() LT 3 then begin
   print, 'Usage: flux = mgeflux(radius, amp_mge [or Ie], sigma_mge [or Re], sigma_y_mge, XC=, YC=, /LUM_IS_AMP, /PRECISE'
   print, '                      /EXP_SDSS, /DEV_SDSS, sersic_n=, QPROF=1.0'
   return, -1
endif

;; Check for too many profiles selected
do_profile_integral = (n_elements(EXP_SDSS) + n_elements(DEV_SDSS) + n_elements(SERSIC_N))
if do_profile_integral GT 1 then begin
   print, 'ERROR: Select one profile shape at a time.  Choices are: /EXP_SDSS, /DEV_SDSS, SERSIC_N='
   return, -1
endif

common gauss2D_common, A0, sigma, sigma_y, qgauss, Rlim, xcen, ycen
common sdss_profs, Ie_Exp, q_Exp, Re_Exp, Ie_Dev, q_Dev, Re_Dev  ;; Ie is surface brightness at Re
common sersic2D_common, I0, Re, n, bn, q

;; Clear common block values, if previously defined
if n_elements(xcen) GT 0 then aa = temporary(xcen)
if n_elements(ycen) GT 0 then aa = temporary(ycen)

n_terms = n_elements(amp)
if keyword_set(XC) then XCEN = fltarr(n_terms) + XC
if keyword_set(YC) then YCEN = fltarr(n_terms) + YC

napers = n_elements(radius)
flux = fltarr(napers)

if do_profile_integral EQ 0 then begin  ;; If true, we're integrating an MGE
   if keyword_set(QPROF) then qgauss = QPROF else $
      if n_elements(sigma_y_mge) GT 0 then qgauss = -1 else qgauss = 1.0

   w_defined = where(amp GT 0, n_defined)  ;; Ignore terms with zero or negative amplitude
   if n_defined EQ 0 then return, -1

   amp = amp[w_defined]
   sigma_mge = sigma_mge[w_defined]

   if keyword_set(LUM_IS_AMP) then amp_to_lum = sqrt(2*!pi*sigma_mge^2) else amp_to_lum = 1.0

   A0 = amp[*] / amp_to_lum
   sigma = sigma_mge[*]       
   if n_elements(sigma_y_mge) GT 0 then sigma_y = sigma_y_mge[w_defined] else sigma_y = sigma
 
   !EXCEPT = 0  ;; Suppress math errors

   if n_elements(xcen) + n_elements(ycen) EQ 0 then begin  ;; No x,y offset
      if radius[0] EQ -1 then begin                             ;; Integrate gaussians to infinity
         if keyword_set(QPROF) then $
            flux = total(2*!pi*A0*sigma^2*qgauss) else flux = total(2*!pi*A0*sigma^2*(sigma_y/sigma))
      endif else begin
         
         for m=0, napers-1 do begin
            Rlim = radius[m]
            if keyword_set(PRECISE) then flux[m] = 2*QROMO('gauss2D_1dint', 0, Rlim) else $
               flux[m] = 2*QROMO('gauss2D_1dint', 0, Rlim, K=3, eps=1e-3)          ;; Fastest is K=3 EPS=1e-3        
         endfor
      endelse

   endif else begin  ;; The aperture is offset in x,y
      for m=0, napers-1 do begin
         Rlim = radius[m]
         if keyword_set(PRECISE) then flux[m] = QROMO('gauss2D_1dint', -Rlim, Rlim) else $
            flux[m] = QROMO('gauss2D_1dint', -Rlim, Rlim, K=3, eps=1e-3)          ;; Fastest is K=3 EPS=1e-3        
      endfor
   endelse
   !EXCEPT = 1
endif else begin  ;; Instead, we're integrating a predefined profile

   ;; Define common block (sdss_profs) variables
   Ie_Exp = amp
   Re_Exp = sigma_mge
   Ie_deV = amp
   Re_deV = sigma_mge

   if keyword_set(QPROF) then begin
      q_Exp = QPROF
      q_Dev = QPROF
      q = QPROF
   endif else begin
      q_Exp = 1.0
      q_Dev = 1.0
      q = 1.0
   endelse

   if keyword_set(EXP_SDSS) then for m=0, napers-1 do flux[m] = INT_2D('fexpcut_radintgrnd', [0,radius[m]], 'thetalims',96, /double)
   if keyword_set(DEV_SDSS) then for m=0, napers-1 do flux[m] = INT_2D('fdevcut_radintgrnd', [0,radius[m]], 'thetalims',96, /double)

   if keyword_set(SERSIC_N) then begin
      n = SERSIC_N      
      bn = newton((2.0*n - 0.324)>0.1, 'sersic_bn')  ;; Exact solution, interpolated
      I0 = amp
      Re = sigma_mge

      for m=0, napers-1 do flux[m] = INT_2D('sersic2D', [0,radius[m]], 'thetalims',96, /double)
   endif

endelse


return, flux

end
