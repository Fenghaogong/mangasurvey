;; ----------------------------------------------------------------------------------------------------
;; NAME: MGSERSIC.PRO
;;
;; PURPOSE: Returns an array of size [2, nterms] multi-Gaussian expansion (MGE) using the MGE_FIT_SECTORS package described in
;; Cappellari (2002, MNRAS, 333, 400) for Sersic profiles with arbitrary n-value and effective radius.  For use with SynMags
;; (synthetic aperture photometry).  Please download MGE_FIT_SECTORS prior to using MGSERSIC.
;;
;; CALLING SEQUENCE: 
;;  mge = mgsersic(n_sersic, Re_sersic, NSTEP=NSTEP, RANGE=RANGE, RARR=RARR, SDSS=SDSS)
;;
;; ARGUMENTS:
;;  n_sersic:   An arbitrary Sersic-n value describing the Sersic shape.  Scalar.
;;  Re_sersic:  Effective radius, arbitrary units since length scales will be scaled later.  For a
;;              radius vector ranging from 1 to 1000, a value of Re_sersic=100 might be appropriate.
;;
;; KEYWORDS:
;;  NSTEP:      Number of samples in radius during the fit.  Integer, default = 300.
;;  RANGE:      Minimum and maximum radii to fit.  2-element vector, should avoid 0.  Default is
;;              [1,1000].
;;  SDSS:       Toggle to use the SDSS variants of the Exponential and de Vaucouleurs profiles.
;;
;; NOTES:
;;  Fits are performed in 1-dimension under assumptions described in Cappellari et al. 2002.  A 2D
;;  implementation as used in Bundy et al. 2012 and descibred by Hogg et al. (in prep) provides
;;  better fits to surface brightness.
;;
;; HISTORY:
;;  Created by K. Bundy - 6/25/12
;; ----------------------------------------------------------------------------------------------------


;; Implicit equation for bn (deV profile)
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

;; SDSS Exp with smoother truncation: Provides better MGE fits to SDSS-like profiles
;;   R: Radius, physical units (same units as those defining Re)
;;   theta: Polar angle, radians
;;   
function fexpcut2, r, theta
  common sdss_profs, Ie_Exp, q_Exp, Re_Exp, Ie_Dev, q_Dev, Re_Dev  ;; Ie is surface brightness at Re
  if n_elements(q_exp) EQ 0 then q_Exp = 1.0

  q_Exp = float(q_Exp)
  rre = r / float(Re_exp)
  rvec = rre * sqrt((1 + (1.0/q_Exp^2 - 1)*(sin(theta)^2)))  ;; Positional vector in polar coordinates, assuming y-direction inclination axis ratio, q

  EXPFAC = -1.67835D  ;; Defines Ie_Exp and Re

  fexp = Ie_Exp * exp(ExpFac * (rvec-1.0))

  ExpOut = 4.0               ;; Profile goes to zero at this radius (units of Re)
  ExpCut = 3.0*ExpOut/4.0    ;; Profile truncation begins at this radius (units of Re)
  ExpCut2 = (ExpCut+ExpOut)/2.0
  
  wr_cutoff = where(rre GT ExpCut, nr_cutoff)
  wr_cutoff2 = where(rre GT ExpCut2, nr_cutoff2)
  wr_zero = where(rre GT ExpOut, nr_zero)

  fexpcut = fexp
  if nr_cutoff GT 0 then fexpcut[wr_cutoff] = fexpcut[wr_cutoff] *  exp(-(rvec[wr_cutoff]-expcut)^2)
  if nr_cutoff2 GT 0 then fexpcut[wr_cutoff2] = fexpcut[wr_cutoff2] *  exp(-(rvec[wr_cutoff2]-expcut2)^2)
  ;; if nr_cutoff GT 0 then fexpcut[wr_cutoff] = fexpcut[wr_cutoff] *  ( 1 - ((rvec[wr_cutoff] - ExpCut)/(ExpOut - ExpCut))^2)^2.0
  ;;if nr_cutoff GT 0 then fexpcut[wr_cutoff] = fexpcut[wr_cutoff] *  ( 1 - ((rvec[wr_cutoff] - ExpCut)/(ExpOut - ExpCut))^2)^1.0  ;; NOTE: To reproduce phFitobj.h, last exponent must be 1.0 (not 2.0)
  ;;if nr_zero GT 0 then fexpcut[wr_zero] = 0

  return, fexpcut
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

;; SDSS deV with smoother truncation: Provides better MGE fits to SDSS-like profiles
;;   R: Radius, physical units (same units as those defining Re)
;;   theta: Polar angle, radians
;;   
function fDeVcut2, r, theta
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
  deVCut2 = (deVCut+deVOut)/2.0

  wr_cutoff = where(rre GT deVCut, nr_cutoff)
  wr_cutoff2 = where(rre GT deVCut2, nr_cutoff2)
  wr_zero = where(rre GT deVOut, nr_zero)

  fdeVcut = fdeV
  if nr_cutoff GT 0 then fdeVcut[wr_cutoff] = fdeVcut[wr_cutoff] *  exp(-(rvec[wr_cutoff]-devcut)^2)
  if nr_cutoff2 GT 0 then fdeVcut[wr_cutoff2] = fdeVcut[wr_cutoff2] *  exp(-(rvec[wr_cutoff2]-devcut2)^2)
  ;;if nr_cutoff GT 0 then fdeVcut[wr_cutoff] = fdeVcut[wr_cutoff] *  ( 1 - ((rvec[wr_cutoff] - deVCut)/(deVOut - deVCut))^2)^1.0 ;; NOTE: To reproduce phFitobj.h, last exponent must be 1.0 (not 2.0)
  ;;if nr_zero GT 0 then fdeVcut[wr_zero] = 0

  return, fdeVcut
end


;; ----------------------------------------------------------------------------------------------------------------------------------
;;  Main Function
;; ----------------------------------------------------------------------------------------------------------------------------------
;; Function decomposes Sersic profiles into Multi-Gaussian expansions (MGE)
function mgsersic, n_sersic, Re_sersic, NSTEP=NSTEP, RANGE=RANGE, RARR=RARR, SDSS=SDSS, REFIT=REFIT

if n_params() LT 2 then begin
   print
   print, ' Usage:  mgsersic, n_sersic, Re_sersic, NSTEP=300, RANGE=[1,1000], RARR=, /SDSS, /REFIT'
   print
   return, -1
endif

if not(keyword_set(NSTEP)) then nstep = 300
if not(keyword_set(RANGE)) then range = [1,1000]

;; MGEs computed by David Hogg and received on June 26th, 2012
;;   Hogg amplitudes are defined as B, where each Gaussian term has the prefactor A = B / (2*!Pi*sigma^2), where A (the prefactor)
;;   is the amplitude definition used here and sigma is the original value before scaling by Re_sersic.  Hogg models are fit to
;;   Exp and deV profiles with Re=1.0 and Ie=1.0.

;;   Standard Profiles, sigmas:
sigma_exp_hogg = Re_sersic*sqrt([2.57275869e-03, 1.89240445e-02, 8.28373784e-02, 2.82980458e-01, 8.31927834e-01, 2.25475415e+00])
sigma_dev_hogg = Re_sersic*sqrt([1.26864828e-06, 2.25834145e-05, 2.13623362e-04, 1.54482026e-03, 9.85339708e-03, 6.10055095e-02, 4.08100570e-01, 3.70795781e+00])

;;   Modified SDSS ("Lupton") Profiles, sigmas: 
sigma_exp_hoggsdss = Re_sersic*sqrt([1.20078965e-03, 8.84526493e-03, 3.91463084e-02, 1.39976817e-01, 4.60962500e-01, 1.50159566e+00])
sigma_dev_hoggsdss = Re_sersic*sqrt([2.23759216e-04, 1.00220099e-03, 4.18731126e-03, 1.69432589e-02, 6.84850479e-02, 2.87207080e-01, 1.33320254e+00, 8.40215071e+00])

;;   Standard Profiles, amplitudes: 
amp_exp_hogg = [7.35002089e-03, 9.48127712e-02, 6.35754740e-01, 2.60086402e+00, 5.42848758e+00, 3.16431709e+00] / (2*!pi*(sigma_exp_hogg/Re_sersic)^2) ;; comma at last entry bug?
amp_dev_hogg = [2.62203068e-03, 2.50014619e-02, 1.34130447e-01, 5.13260990e-01, 1.52005095e+00, 3.56204929e+00, 6.44845002e+00, 8.10104509e+00] / (2*!pi*(sigma_dev_hogg/Re_sersic)^2)  ;; comma at last entry bug?

;;   Modified SDSS ("Lupton") Profiles, amplitudes: 
amp_exp_hoggsdss = [2.34853813e-03, 3.07995260e-02, 2.23364214e-01, 1.17949102e+00, 4.33873750e+00, 5.99820770e+00] / (2*!pi*(sigma_exp_hoggsdss/Re_sersic)^2)
amp_dev_hoggsdss = [4.26347652e-02, 2.40127183e-01, 6.85907632e-01, 1.51937350e+00, 2.83627243e+00, 4.46467501e+00, 5.72440830e+00, 5.60989349e+00] / (2*!pi*(sigma_dev_hoggsdss/Re_sersic)^2)

if keyword_set(SDSS) and not(n_sersic EQ 1.0 OR n_sersic EQ 4.0) then begin
   print, ' MGSERSIC ERROR: Modified SDSS profiles must be either Exp of deV (n=1 or 4)'
   return, -1
endif

case n_sersic OF 
   1.0: begin
      if keyword_set(SDSS) then sol = transpose([[amp_exp_hoggsdss], [sigma_exp_hoggsdss]]) else $
         sol = transpose([[amp_exp_hogg], [sigma_exp_hogg]])
   end
   4.0: begin
      if keyword_set(SDSS) then sol = transpose([[amp_dev_hoggsdss], [sigma_dev_hoggsdss]]) else $
         sol = transpose([[amp_dev_hogg], [sigma_dev_hogg]])
   end

   else: begin
      if not(keyword_set(RANGE)) then RANGE = [1,1000]
      R = getrange(RANGE[0],RANGE[1],nstep,/LOG)

      ;;model_sersic =  (pixscale)^2*I0_c * exp( -bn_c*(im_circle*pixscale / Re_c)^(1.0/n_c)*(1 + (1.0/(q_c)^2-1)*sin2_theta)^(1.0/(2*n_c)) )
      
      ;;  Determine exact "bn" numerically 
      common sersic2D_common, I0, Re, n, bn, q
      n = n_sersic
      Re = Re_sersic
      bn = newton((2*n - 0.324)>0.1, 'sersic_bn')  ;; Exact solution, interpolated
      
      profile = exp( -bn*((R/Re_sersic)^(1.0/n_sersic)-1) )  ;; Makes I0 into Ie (surface brightness at Re)
      
      mge_fit_1d, R, profile, NGAUSS=16, SOL=sol
      
      ;; Change first dimension output from 1D flux integral to (total) amplitudes
      sol[0,*] = sol[0,*] / sqrt(2*!pi*sol[1,*]^2)
      
      RARR = R
   end
endcase

if keyword_set(REFIT) then begin
   if keyword_set(SDSS) then begin
      common sdss_profs, Ie_Exp, q_Exp, Re_Exp, Ie_Dev, q_Dev, Re_Dev  ;; Ie is surface brightness at Re

      case n_sersic OF
         1.0: begin
            Ie_Exp = 1.0
            Re_Exp = Re_sersic
            
            ;; RANGE = [1, Re_Exp*3.95] ;; This range provides decent fits to fExpCut
            RANGE = [1, Re_Exp*4.1] ;; For fits to fExpCut2: Range set to provide better fits near the truncation radius
            R = getrange(RANGE[0],RANGE[1],nstep,/LOG)
            
            ;; profile = fExpCut(R, 0) ;; SDSS-modified truncated profile
            profile = fExpCut2(R, 0) ;; Less drastic truncation -- enables better MGE fits
            
         end
         4.0: begin
            Ie_Dev = 1.0
            Re_Dev = Re_sersic
            
            ;; RANGE = [1, Re_Dev*7.95] ;; This range provides decent fits to fDevCut
            RANGE = [1, Re_Dev*8.1] ;; For fits to fDevCut2: Range set to provide better fits near the truncation radius
            R = getrange(RANGE[0],RANGE[1],nstep,/LOG)
            
            ;;profile = fDevCut(R, 0)  ;; SDSS-modified truncated profile
            profile = fDevCut2(R, 0)  ;; Less drastic truncation -- enables better MGE fits
            
         end
         else: begin
            print, 'Modified SDSS profiles must be either Exp of deV (n=1 or 4)'
            return, -1
         end
      endcase
   endif else begin
      if not(keyword_set(RANGE)) then RANGE = [1,1000]
      R = getrange(RANGE[0],RANGE[1],nstep,/LOG)
      
      ;;model_sersic =  (pixscale)^2*I0_c * exp( -bn_c*(im_circle*pixscale / Re_c)^(1.0/n_c)*(1 + (1.0/(q_c)^2-1)*sin2_theta)^(1.0/(2*n_c)) )
      
      ;;  Determine exact "bn" numerically 
      common sersic2D_common, I0, Re, n, bn, q
      n = n_sersic
      Re = Re_sersic
      bn = newton((2*n - 0.324)>0.1, 'sersic_bn')  ;; Exact solution, interpolated
      
      ;;profile = exp( -bn*( R / Re_sersic)^(1.0/n_sersic) )
      profile = exp( -bn*((R/Re_sersic)^(1.0/n_sersic)-1) )  ;; Makes I0 into Ie (surface brightness at Re)
   endelse
   
   mge_fit_1d, R, profile, NGAUSS=16, SOL=sol
   
   ;; Change first dimension output from 1D flux integral to (total) amplitudes
   sol[0,*] = sol[0,*] / sqrt(2*!pi*sol[1,*]^2)
   
   RARR = R
endif
   
return, sol
end

