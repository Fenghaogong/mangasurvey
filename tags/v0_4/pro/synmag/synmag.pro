;; ----------------------------------------------------------------------------------------------------
;; NAME: SYNMAG.PRO
;;
;; PURPOSE: Procedure predicts aperture magnitudes given information about the PSF and the intrinsic profiles of each object.  See
;; Bundy et al. 2012 for details.
;;
;; CALLING SEQUENCE: 
;;   synmag, ModelScale, ModelRad, ModelAB, PSF_FWHM_info, aper_radius, synmag, $
;;           N_SERSIC = N_SERSIC, SDSS = SDSS, ZPT_FOR_I0 = ZPT_FOR_I0, $
;;           MGE_intrinsic = MGE_intrinsic, MGEOut = MGEOut
;;
;; ARGUMENTS: 
;;   ModelScale:        Input array of total profile-fit magnitudes or, if ZPT_FOR_I0, is set, the surface brightness at Re
;;   PSF_FWHM_info:     Input Gaussian FWHM seeing (arcsec) to be convolved with intrinsic profiles. Can be of the following form:
;;                        1) A scalar, in which case a single Gaussian with this FWHM is convolved with all objects;
;;                        2) A vector with size equal to the number of objects, allowing a different FWHM for each object;
;;                        3) An array with dimensions [nmags, 2, n_PSF], where nmags is the number of objects, n_PSF is the number
;;                        of Gaussian PSF components, and the second dimension specifies each component's Gaussian amplitude and sigma.
;;   aper_radius:       Input radii of apertures in arcsec.  All apertures are predicted for all magnitudes
;;   synmag:            Output synthetic aperture magnitude
;;
;; KEYWORDS:
;;   N_SERSIC:        Set to the Sersic-n value of the assumed intrinsic profile shape (if described as a Sersic profile)
;;   SDSS:            Toggle for using the modified SDSS profiles for n=1 (Exponential) or n=4 (de Vaucouleurs)
;;   ZPT_FOR_I0:      If ModelScale is the surface brightness at Re, this required keyword determines the synmag zeropoint
;;   MGE_INTRINSIC:   Array of size [2, n_terms].  The first element of each term is the full amplitude, Amp, defined such that
;;                    the total 2D flux is Amp*2*!pi*sigma^2.  The second element is the Gaussian sigma of each MGE component.  On
;;                    input, MGE_INTRINSIC defines the intrinsic profile assumed.  On output, it returns the MGE that was used.
;;   MGEOUT           Output structure providing the PSF-convolved MGE information.  For easy input to MGESTR.PRO;;   
;;
;; NOTES:
;;   Flux Normalization: If ModelScale is a total magnitude, SYNMAG computes the flux normalization, that is the amplitude of the
;;   fitted profile, by comparing the total flux implied by ModelScale to the 2D integral of the fitted intrinsic profile when its
;;   amplitude is unity.  The integral is taken to infinity which may introduce small biases in the case where ModelScale is
;;   defined by integrating to a finite, but large, radius.  In the case of SDSS, for example, use SYNMAG_SDSS.
;;
;;
;; HISTORY:
;;   Created by KBundy on 8/11/2012: Inherited from sdss_to_aper6.pro
;;
;;
pro synmag, ModelScale, ModelRad, ModelAB, PSF_FWHM_info, aper_radius, synmag, N_SERSIC = N_SERSIC, $
            SDSS = SDSS, ZPT_FOR_I0 = ZPT_FOR_I0, $
            MGE_intrinsic = MGE_intrinsic, MGEOut = MGEOut

if not(keyword_set(FullRad_Multiplier)) then FullRad_multiplier = 10.0 * ModelRad ;; Radius out to which "total flux" is integrated
if not(keyword_set(MGE_scale_sampling)) then MGE_scale_sampling = 100.0  ;; Scale MGE spatial dimension, defines intrinsic profile's radius value at Re; 1 Re = 100 spatial bins

nmags = n_elements(ModelScale)
napers = n_elements(aper_radius)
if n_elements(aper_radius) GT 1 then synmag = fltarr(nmags, napers)-99 else synmag = fltarr(nmags)-99

;; ..................................................
;; Set up PSF information

PSF_FWHM_array_n_dim_orig = size(PSF_FWHM_info, /n_dim)
if PSF_FWHM_array_n_dim_orig LE 1 then begin
   if n_elements(PSF_FWHM_info) EQ 1 then PSF_FWHM_array = fltarr(nmags) + PSF_FWHM_info
   n_PSF = 1
   PSF_sigma = fltarr(nmags) + PSF_FWHM_info / 2.3548 ;; The gaussian sigma of the PSF
   R_PSF = fltarr(nmags) + 1./(2*!pi*PSF_sigma^2)      ;; Amplitude of each PSF Gaussian term
endif

if PSF_FWHM_array_n_dim_orig EQ 2 then begin ;; This is a Multi-Gaussian PSF [2, n_PSF], but same for all sources. 
   n_PSF = n_elements(PSF_FWHM_info[0,*])

;;    PSF_FWHM_full = fltarr(2, n_PSF, nmags)
;;    for k=0,n_PSF-1 do begin                       ;; We need PSF_FWHM_full -> [2, n_PSF, nmags]
;;       replicate_inplace, PSF_FWHM_full, PSF_FWHM_array[0,k], 3, [0,k,0]
;;       replicate_inplace, PSF_FWHM_full, PSF_FWHM_array[1,k], 3, [1,k,0]
;;    endfor

   PSF_FWHM_full = fltarr(nmags, 2, n_PSF)
   for k=0,n_PSF-1 do begin                       ;; We need PSF_FWHM_full -> [nmags, 2, n_PSF]
      replicate_inplace, PSF_FWHM_full, PSF_FWHM_info[0,k], 1, [0,0,k]
      replicate_inplace, PSF_FWHM_full, PSF_FWHM_info[1,k], 1, [0,1,k]
   endfor

   PSF_FWHM_array = PSF_FWHM_full  ;; Replace with full array
   R_PSF = reform(PSF_FWHM_array[*,0,*])            ;; Amplitude of each PSF Gaussian term
   PSF_sigma = reform(PSF_FWHM_array[*,1,*])/2.3548 ;; The gaussian sigma of the PSF
endif


if PSF_FWHM_array_n_dim_orig EQ 3 then begin     ;; We were given a full Multi-Gaussian PSF array, defined for all sources (with dimensions, [nmags, 2, n_PSF]
   R_PSF = reform(PSF_FWHM_info[*,0,*])            ;; Amplitude of each PSF Gaussian term
   PSF_sigma = reform(PSF_FWHM_info[*,1,*])/2.3548 ;; The gaussian sigma of the PSF
endif



;; ..................................................
;; Determine MGE approximation to the intrinsic profile (if needed)
if n_elements(MGE_intrinsic) LT 2 then begin  ;; No MGE for the intrinsic profile was passed

   if n_elements(n_sersic) NE 1 then begin
      print
      print, ' SYNMAG: No MGE Passed.  To use a Sersic intrinsic profile, please set N_SERSIC keyword.'
      print
   endif

   mge = mgsersic(n_sersic, MGE_scale_sampling, NSTEP=NSTEP, RANGE=RANGE, RARR=RARR, SDSS=SDSS)
   mge[1,*] = mge[1,*]/MGE_scale_sampling

   mge_intrinsic = mge

   amp_mge = mge[0,*]    ;; MGE amplitdues
   sigma_mge = mge[1,*]  ;; MGE sigma values

endif else begin
   print
   print, ' SYNMAG: Using passed MGE for the intrinsic profile.'
   print

   mge = mge_intrinsic
   amp_mge = mge[0,*]                  ;; MGE amplitdues
   sigma_mge = mge[1,*]                ;; MGE sigma values
endelse

wgood = where(ModelScale GT 0 AND ModelRad GT 0, ngood)

badflag = bytarr(nmags) + 1  ;; There's a problem if badflag=1
badflag[wgood] = 0

;; The following assumes regular magnitudes -- Luptitudes are different by ~1% at mags fainter than 22.3
if keyword_set(ZPT_FOR_I0) then dummy_zpt = ZPT_FOR_I0 else begin
   ModelMag = ModelScale
   dummy_zpt = 22.5
   total_flux = 10^(0.4*(dummy_zpt-ModelMag))
endelse


;; Structure for holding output convolved MGE information
n_terms_mge = n_elements(mge[0,*])

mgeout = {psf_gauss:fltarr(2, n_psf)-1, $
          A0:fltarr(n_terms_mge*n_psf)-99, sigma_x:fltarr(n_terms_mge*n_psf)-99, sigma_y:fltarr(n_terms_mge*n_psf)-99}
mgeout = replicate(mgeout, nmags)

;; fxbhmake, hdr_mgeout, nmags, /initialize, /date
;; sxaddpar, hdr_mgeout, 'LEN_UNIT', 'Arcsec', 'Units of sigma values'
;; sxaddpar, hdr_mgeout, 'COEF', 'Amplitude', 'Type of MGE normalization coefficient'
;; sxaddpar, hdr_mgeout, 'Pref', '1=Exp, 2=Dev', 'Type of MGE normalization coefficient'
;; sxaddpar, hdr_mgeout, 'N_EXP', n_terms_exp, 'Number of terms in the Exponential MGE'
;; sxaddpar, hdr_mgeout, 'N_DEV', n_terms_dev, 'Number of terms in the deVaucouleurs MGE'

time_integrate = fltarr(nmags)  ;; Time it takes to do MGE flux integrals
time_mock = fltarr(nmags)       ;; Time it takes to do mock image with filtering

A_invft_mge = fltarr(n_terms_mge, n_PSF)  ;; Output PSF-convolved MGE amplitudes
sigma_array = fltarr(n_terms_mge, n_PSF)  ;; Output PSF-convolved major-axis sigmas
sigma_y_array =fltarr(n_terms_mge, n_PSF) ;; Output PSF-convolved minor-axis sigmas

;; ----------------------------------------------------------------------------------------------------
;; Loop over all objects 
for i=0L, nmags-1 do begin  
   s1 = systime(/sec)
   if i GT 0 AND (i - 100*(i/100) EQ 0 OR i EQ nmags-1) then print, ' SYNMAG: Number '+strc(i)+' of '+strc(nmags-1)+'.  Time left: '+dplace(time_left,1)+' min'
   if badflag[i] GE 1 then goto, skip

   ;; ..................................................
   ;; Determine normalization of the flux model
   if keyword_set(ZPT_FOR_I0) then I0 = ModelScale[i] else begin
      ;; Figure out the flux normalization from the total (model) magnitude provided
            
      sigma_y_mge = ModelAB[i] * sigma_mge 

      flux_norm = mgeflux(-1, amp_mge, ModelRad[i]*sigma_mge, ModelRad[i]*sigma_y_mge, /precise)

      ;flux_norm = mgeflux(4.0*ModelRad[i], 1.0, ModelRad[i], qprof=ModelAB[i], /precise, /Exp_SDSS)

      I0 = total_flux[i] / flux_norm 
   endelse

   ;; ----------------------------------------------------------------------------------------------------
   ;; Use the MGE expansion to analytically describe the convolved profile
   mgeout[i].psf_gauss[0,*] = R_PSF[i,*]
   mgeout[i].psf_gauss[1,*] = psf_sigma[i,*]

   q_obj = ModelAB[i]

   sigma_obj = sigma_mge * (ModelRad[i])  ;; Object's instrinsic MGE gaussian sigmas for object in question (arcsec)
   ;A_mge = I0 * (amp_mge / sqrt(2*!pi*sigma_mge^2))            ;; Object's intrinsic MGE gaussian coeffecients of the expansion
   A_mge = I0[0] * amp_mge             ;; Object's intrinsic MGE gaussian coeffecients of the expansion

   ;; Loop over the Gaussian components of the PSF
   for k=0, n_PSF-1 do begin
      A_invft_mge[*,k] = (A_mge * R_PSF[i,k] * 2*!pi * sigma_obj^2 * q_obj * psf_sigma[i,k]^2) / $
                             sqrt( (sigma_obj^2 + psf_sigma[i,k]^2)*(sigma_obj^2*q_obj^2 + psf_sigma[i,k]^2) ) ;; Array of Prefactors (on convolved gaussian term)
      
      sigma_array[*,k] = sqrt(sigma_obj^2 + psf_sigma[i,k]^2)
      sigma_y_array[*,k] = sqrt(sigma_obj^2*q_obj^2 + psf_sigma[i,k]^2)
   endfor
   
   t0_integrate = systime(/sec)
      
   ;; Integrate the summed convolved MGE numerically....    
   flux_sum = mgeflux(aper_radius, A_invft_mge[*], sigma_array[*], sigma_y_array[*])

   ;; Store MGE info output
   mgeout[i].A0 = A_invft_mge[*]
   mgeout[i].sigma_x = sigma_array[*]
   mgeout[i].sigma_y = sigma_y_array[*]
   
   synmag[i,*] = dummy_zpt - 2.5*alog10(flux_Sum)

   tend_integrate = systime(/sec)
   time_integrate[i] = tend_integrate - t0_integrate
   
   ;stop

   skip:
   s2 = systime(/sec)
   time_left = (s2-s1)*(nmags-1-i)/60 ; Minutes   
endfor


end
