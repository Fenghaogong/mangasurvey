;; ----------------------------------------------------------------------------------------------------
;; NAME: MGESTR.PRO
;;
;; PURPOSE: Given an MGEs for a set of objects, MGESTR computes information about the structure of the profile.  This includes the
;; major-axis effective radius, the half-light radius within a circular aperture, and the effective axis ratio.  See Bundy et
;; al. 2012 for further details.
;;
;; CALLING SEQUENCE: 
;;   mgestr, Re_mge, q_mge, re_aper, A_mge = A_mge, sigma_x_mge = sigma_x_mge, sigma_y_mge = sigma_y_mge, mge_struct = mge_struct            
;;
;;
;; OUTPUT:
;;   Re_MGE    On output, returns a vector of effective radii along the major axis
;;   Q_MGE     An (optional) output vector of the effective axis ratios of summed MGE profiles.  
;;             SIGMA_Y must be provided or else Q_MGE = 1.0.
;;   Re_APER   An (optional) output vector of the effective radii of circular apertures containing half the total
;;             flux. SIGMA_Y_MGE must be provided or else this is equivalent to RE_MGE.
;;            
;; INPUT:
;;   A_MGE       An input array of vectors holding the normalization coefficients (amplitudes) of each MGE term.  Dimensions are [N_MGE,
;;               N_terms], where N_MGE is the number of input profiles.
;;   SIGMA_X     An input array of vectors holding the major-axis sigma values of each MGE term.  Dimensions are [N_MGE, N_terms]
;;   SIGMA_Y     An (optional) input array of vectors holding the minor-axis sigma values of each MGE term.  Dimensions are [N_MGE,
;;               N_terms].  
;;   MGE_STRUCT  A structure providing A_MGE, SIGMA_X, and SIGMA_Y.  If provided, MGE_STRUCT overrides the input keywords above.
;;               MGE_STRUCT is provided, e.g., by the SYNMAG routine
;;
;; HISTORY:
;;   Created by Kevin Bundy 5/15/12
;;

;; Define a function whose root occurs when the 2D summed MGE flux within a circular aperture equals half the total flux
function mgeaper_half_root, rr
  common mge_common, A0, sigma, sigma_y, Rlim, F_half

  n_rr = n_elements(rr)
  output = fltarr(n_rr)

  for ii=0, n_rr-1 do begin
     Rlim = rr[ii]
     output[ii] = F_half - 2*QROMO('aperint_1D', 0, rr, K=3, eps=1e-3)          ;; Fastest is K=3 EPS=1e-3   
  endfor

  return, output
end

;; Define the 1D integrand ror 2D Gaussian integral in Cartesian coordinates: 
function aperint_1D, x
  common mge_common, A0, sigma, sigma_y, Rlim, F_half
  n_sum = n_elements(A0)

  if n_sum EQ 1 then integrand = 2*A0 * sigma_y * sqrt(!pi/2.0) * exp( -x^2/(2.0*sigma^2)) * erf(sqrt( (Rlim^2-x^2)/(2.0*sigma_y^2) )) else $
     integrand = total(2*A0 * sigma_y * sqrt(!pi/2.0) * exp( -x^2/(2.0*sigma^2)) * erf(sqrt( (Rlim^2-x^2)/(2.0*sigma_y^2) )))

  return, integrand
end

;; Assuming a fixed axis ratio of 1.0, with sigma's defined along the major-axis, define a function whose root occurs when the
;; resulting 2D summed MGE flux reaches half the total flux (under this assumption)
function mgemaj_half_root, x
  common mge_common, A0, sigma, sigma_y, Rlim, F_half

  n_x = n_elements(x)
  output = fltarr(n_x)

  for ii=0, n_x-1 do output[ii] = F_half - total(2*!pi * A0 * sigma^2 * (1 - exp(-x[ii]^2 / (2 * sigma^2))))
  return, output
end

;; Procedure reads output MGE structure and returns desired MGE amplitudes and major + minor axis sigmas
pro read_mge_struct, mge_struct, A_mge, sigma_x_mge, sigma_y_mge

  if size(mge_struct, /type) NE 8 then begin
     print, 'ERROR: There is a problem with the passed mge_struct'
     return
  endif

  n_mge = n_elements(mge_struct)
  n_terms_max = n_elements(mge_struct[0].A0)

  A_mge = transpose(mge_struct.A0)
  sigma_x_mge = transpose(mge_struct.sigma_x)
  sigma_y_mge = transpose(mge_struct.sigma_y)
end


;; ----------------------------------------------------------------------------------------------------
;; MGEstr - main procedure
;; ----------------------------------------------------------------------------------------------------
pro mgestr, Re_mge, q_mge, re_aper, A_mge = A_mge, sigma_x_mge = sigma_x_mge, sigma_y_mge = sigma_y_mge, mge_struct = mge_struct            

common mge_common, A0, sigma, sigma_y, Rlim, F_half

usage_message1 = 'MGEstr Usage, output arguments: Re_mge, [q_mge], [Re_aper] '
usage_message2 = '              input keywords  : A_mge=, sigma_x=, [sigma_y=] OR mge_struc='

if (n_params() EQ 0 OR n_params() GT 3) then begin
   print, usage_message1
   print, usage_message2
   return
endif


;; Figure out the type and method of input being passed, do some error checking
if n_elements(mge_struct) GT 0 then begin
   print, 'Using the passed structure of MGE values as input'

   read_mge_struct, mge_struct, A_mge, sigma_x_mge, sigma_y_mge
endif else begin
   if not(keyword_set(A_mge)) OR not(keyword_set(sigma_x_mge)) then begin
      print, 'Please supply (at least) A_mge and sigma_x keyword input'
      print, usage_message1
      print, usage_message2
   endif
endelse

size = size(reform(A_mge))
if size[0] EQ 2 then n_mge = size[1] else n_mge = 1  ;; Test how many MGEs were passed

;; Setup output arrays
Re_mge_guess = fltarr(n_mge)
Re_mge_maj = fltarr(n_mge)
Re_mge_aper = fltarr(n_mge)
q_mge = fltarr(n_mge)

;; ----------------------------------------------------------
;; Estimate the major axis half-light radius
;; ----------------------------------------------------------
!except=0
for i=0L, n_mge-1 do begin   
   w_defined = where(A_mge[i,*] GT 0, n_defined)  ;; Ignore terms with zero or negative amplitude

   if n_defined GT 0 then begin
      A0 = A_mge[i,w_defined]
      sigma = sigma_x_mge[i,w_defined] > 1e-3
      ;;sigma_y = sigma_y_mge[i,w_defined]  ;; Don't need this for Re_mge_maj

      F_half = total(2*!pi * A0 * sigma^2)/2  ;; Half-flux sum
      
      Ri = sqrt(1.3863 * sigma^2)
      Re_mge_guess[i] = total(A0 * Ri^3)/total(A0*Ri^2)
      
      Re_mge_maj[i] = newton(Re_mge_guess[i], 'mgemaj_half_root')
   endif
endfor
!except=1


;; ----------------------------------------------------------
;; Estimate the effective axis ratio if desired
;; ----------------------------------------------------------
;;   We compute the mean axis ratio of all MGE terms, weighted by the 2D flux of each term
if n_params() GE 2 then begin
   for i=0L, n_mge-1 do begin
      w_defined = where(A_mge[i,*] GT 0, n_defined)  ;; Ignore terms with zero or negative amplitude

      if n_defined GT 0 then begin
         A0 = A_mge[i,w_defined]
         sigma = sigma_x_mge[i,w_defined] > 1e-3
         sigma_y = sigma_y_mge[i,w_defined]  > 1e-3
         
         q_hat = sigma_y / sigma
         weights = A0 * sigma * sigma_y
         
         q_mge[i] = total(weights * q_hat) / total(weights)
      endif
      
   endfor
endif

;; ----------------------------------------------------------
;; Measure the half-light circular aperture radius if desired
;; ----------------------------------------------------------
;;
!except=0
if n_params() EQ 3 then begin
   for i=0L, n_mge-1 do begin
   w_defined = where(A_mge[i,*] GT 0, n_defined)  ;; Ignore terms with zero or negative amplitude

   if n_defined GT 0 then begin
      A0 = A_mge[i,w_defined]
      sigma = sigma_x_mge[i,w_defined] > 1e-3
      sigma_y = sigma_y_mge[i,w_defined]  > 1e-3
      
      F_half = total(2*!pi * A0 * sigma * sigma_y)/2  ;; Half-flux sum
      
      Re_mge_aper[i] = newton(Re_mge_maj[i], 'mgeaper_half_root')
   endif
     
   endfor
endif
!except=1

;; Store output vectors
Re_MGE = Re_MGE_Maj
re_aper = Re_mge_aper


end
