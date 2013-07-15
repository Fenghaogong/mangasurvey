;; ----------------------------------------------------------------------------------------------------
;; NAME: SYNMAG_SDSS.PRO
;;
;; PURPOSE: Procedure predicts aperture magnitudes from the SDSS photometry catalog, typically in the r-band, given information
;; about the PSF and the intrinsic profiles of each object.  See Bundy et al. 2012 for details and SYNMAG for a more general
;; routine.
;;
;; CALLING SEQUENCE: 
;;     synmag_sdss, Modelmag, DevMag, DevRad, DevAB, ExpMag, ExpRad, ExpAB, PSF_FWHM_info, aper_radius, synmag, $
;;                  Ie_Dev = Ie_Dev, Ie_Exp = Ie_Exp, mgeout = mgeout, hdr_mgeout = hdr_mgeout
;;
;; ARGUMENTS: 
;;   ModelMag:          Input array of SDSS model magnitudes, typically in the r-band
;;   DevMag:            Input de Vaucouleurs magnitudes
;;   DevRad:            Input de Vaucouleurs radii
;;   DevAB:             Input de Vaucouleurs axis ratios
;;   ExpMag:            Input Exponential magnitudes
;;   ExpRad:            Input Exponential radii
;;   ExpAB:             Input Exponential axis ratios
;;   PSF_FWHM_info:     Input Gaussian FWHM seeing (arcsec) to be convolved with intrinsic profiles. Can be of the following form:
;;                        1) A scalar, in which case a single Gaussian with this FWHM is convolved with all objects;
;;                        2) A vector with size equal to the number of objects, allowing a different FWHM for each object;
;;                        3) An array with dimensions [nmags, 2, n_PSF], where nmags is the number of objects, n_PSF is the number
;;                        of Gaussian PSF components, and the second dimension specifies each component's Gaussian amplitude and sigma.
;;   aper_radius:       Input radii of apertures in arcsec.  All apertures are predicted for all magnitudes
;;   synmag:            Output synthetic aperture magnitude
;;
;; KEYWORDS:
;;   Ie_Dev           Intrinsic de Vaucouleurs profile normalizations, defined as the surface brightness at Re (not calculated if input), default nanomaggies (zpt=22.5)
;;   Ie_Exp           Intrinsic Exponential profile normalizations, defined as the surface brightness at Re (not calculated if input), default nanomaggies (zpt=22.5)
;;   MGEOUT           Output structure providing the PSF-convolved MGE information.  For easy input to MGESTR.PRO
;;   HDR_MGEOUT       Output header for use in saving MGEOUT as a FITS table, e.g., 
;;                      IDL> mwrfits, MGEOUT, 'filename.fits', HDR_MGEOUT, /create
;; 
;; PROCEDURE:
;;   Profile Choice: SYNMAG_SDSS determines whether an Exponential or de Vaucouleurs profile was used to fit each model magnitude
;;   by determining whether the input ModelMag better matches DevMag or ExpMag.
;;
;;
;; HISTORY: v1.0
;;   Created by KBundy on 8/11/2012: Inherited from sdss_to_aper6.pro
;;   v1.1: Changed variable names with I0 to Ie
;;
;;

pro synmag_sdss, Modelmag, DevMag, DevRad, DevAB, ExpMag, ExpRad, ExpAB, PSF_FWHM_info, aper_radius, synmag, $
                   Ie_Dev = Ie_Dev, Ie_Exp = Ie_Exp, mgeout = mgeout, hdr_mgeout = hdr_mgeout

if n_params() EQ 0 then begin
   print
   print, 'Usage: sdss_synmag, Modelmag, DevMag, DevRad, DevAB, ExpMag, ExpRad, ExpAB, PSF_FWHM_array, aper_radius, synmag, $
   print, '       Ie_Dev=, Ie_Exp=, mgeout=, hdr_mgeout=
   return
endif

if not(keyword_set(MGE_scale_sampling)) then MGE_scale_sampling = 100.0  ;; Scale MGE spatial dimension, defines intrinsic profile's radius value at Re; 1 Re = 100 spatial bins

nmags = n_elements(ModelMag)
napers = n_elements(aper_radius)
if n_elements(aper_radius) GT 1 then synmag = fltarr(nmags, napers)-99 else synmag = fltarr(nmags)-99

if not(keyword_set(Ie_Exp)) then Ie_Exp = fltarr(nmags)-1
if not(keyword_set(Ie_Dev)) then Ie_Dev = fltarr(nmags)-1

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
;; Setup MGE approximation to intrinsic SDSS profiles
mge_exp = mgsersic(1.0, MGE_scale_sampling, /SDSS)
mge_dev = mgsersic(4.0, MGE_scale_sampling, /SDSS)

amp_mge_exp = mge_exp[0,*]   ;; MGE amplitdues
sigma_mge_exp = mge_exp[1,*] ;; MGE sigma values

amp_mge_dev = mge_dev[0,*]   ;; MGE amplitdues
sigma_mge_dev = mge_dev[1,*] ;; MGE sigma values

;; ..................................................
;; Determine which objects use which profiles
prof_to_use = bytarr(nmags)

wExp = where(Modelmag GT 0 AND ExpMag GT 0 and ExpRad GT 0 AND $
             abs(ExpMag - ModelMag) LT abs(DevMag - Modelmag), nExp)
wDev = where(Modelmag GT 0 AND DevMag GT 0 and DevRad GT 0 AND $
             abs(DevMag - ModelMag) LT abs(ExpMag - Modelmag), nDev)

if nExp GT 0 then prof_to_use[wExp] = 1.0
if nDev GT 0 then prof_to_use[wDev] = 4.0


;; The following assumes regular magnitudes -- Luptitudes are different by ~1% at mags fainter than 22.3
if keyword_set(ZPT_FOR_Ie) then dummy_zpt = ZPT_FOR_Ie else begin
   dummy_zpt = 22.5
   total_flux = 10^(0.4*(dummy_zpt-ModelMag))
endelse

;; Structure for holding output convolved MGE information
n_terms_exp = n_elements(mge_exp[0,*])
n_terms_dev = n_elements(mge_dev[0,*])
n_terms_max = max([n_terms_exp, n_terms_dev])
mgeout = {psf_gauss:fltarr(2, n_psf)-1, pref:-1, $
          A0:fltarr(n_terms_max*n_psf)-99, sigma_x:fltarr(n_terms_max*n_psf)-99, sigma_y:fltarr(n_terms_max*n_psf)-99}
mgeout = replicate(mgeout, nmags)

mgeout.pref = prof_to_use

fxbhmake, hdr_mgeout, nmags, /initialize, /date
sxaddpar, hdr_mgeout, 'LEN_UNIT', 'Arcsec', 'Units of sigma values'
sxaddpar, hdr_mgeout, 'COEF', 'Amplitude', 'Type of MGE normalization coefficient'
sxaddpar, hdr_mgeout, 'Pref', '1=Exp, 4=Dev', 'Type of MGE normalization coefficient'
sxaddpar, hdr_mgeout, 'N_EXP', n_terms_exp, 'Number of terms in the Exponential MGE'
sxaddpar, hdr_mgeout, 'N_DEV', n_terms_dev, 'Number of terms in the deVaucouleurs MGE'

time_integrate = fltarr(nmags)  ;; Time it takes to do MGE flux integrals
time_mock = fltarr(nmags)       ;; Time it takes to do mock image with filtering

A_invft_mge = fltarr(n_terms_max, n_PSF)  ;; Output PSF-convolved MGE amplitudes
sigma_array = fltarr(n_terms_max, n_PSF)  ;; Output PSF-convolved major-axis sigmas
sigma_y_array =fltarr(n_terms_max, n_PSF) ;; Output PSF-convolved minor-axis sigmas

s1 = systime(/sec)
;; ----------------------------------------------------------------------------------------------------
;; Loop over all objects 
for i=0L, nmags-1 do begin  
   if i GT 0 AND (i - 500*(i/500) EQ 0 OR i EQ nmags-1) then print, ' SYNMAG: Number '+strc(i)+' of '+strc(nmags-1)+'.  Time left: '+dplace(time_left,1)+' min'

   if PSF_sigma[i,0] LT 0 then goto, skip  ;; If first Gaussian seeing FWHM is negative then we skip the object

   case prof_to_use[i] OF
      1.0: begin  ;; ModelMag profile fit is an exponential
         if Ie_Exp[i] LE 0 then begin
            flux_norm = mgeflux(4.0*ExpRad[i], 1.0, ExpRad[i], qprof=ExpAB[i], /Exp_SDSS)
            Ie_Exp[i] = total_flux[i] / flux_norm
         endif
   
         q_obj = ExpAB[i]
         sigma_obj = sigma_mge_exp * (ExpRad[i] / MGE_scale_sampling) 
         A_mge = Ie_Exp[i] * amp_mge_exp

         n_terms_use = n_terms_exp
      end
      4.0: begin  ;; Modelmag profile fit is a de Vaucouleurs
         if Ie_Dev[i] LE 0 then begin
            flux_norm = mgeflux(8.0*DevRad[i], 1.0, DevRad[i], qprof=DevAB[i], /Dev_SDSS)
            Ie_Dev[i] = total_flux[i] / flux_norm
         endif
   
         q_obj = DevAB[i]
         sigma_obj = sigma_mge_dev * (DevRad[i] / MGE_scale_sampling) 
         A_mge = Ie_dev[i] * amp_mge_dev

         n_terms_use = n_terms_dev
      end
      else: goto, skip
   endcase


   ;; ----------------------------------------------------------------------------------------------------
   ;; Use the MGE expansion to analytically describe the convolved profile
   A_invft_mge = fltarr(n_terms_use, n_PSF)  ;; Output PSF-convolved MGE amplitudes
   sigma_array = fltarr(n_terms_use, n_PSF)  ;; Output PSF-convolved major-axis sigmas
   sigma_y_array =fltarr(n_terms_use, n_PSF) ;; Output PSF-convolved minor-axis sigmas

   ;; Loop over the Gaussian components of the PSF
   for k=0, n_PSF-1 do begin
      A_invft_mge[*,k] = (A_mge * R_PSF[i,k] * 2*!pi * sigma_obj^2 * q_obj * psf_sigma[i,k]^2) / $
                             sqrt( (sigma_obj^2 + psf_sigma[i,k]^2)*(sigma_obj^2*q_obj^2 + psf_sigma[i,k]^2) ) ;; Array of Prefactors (on convolved gaussian term)
      
      sigma_array[*,k] = sqrt(sigma_obj^2 + psf_sigma[i,k]^2)
      sigma_y_array[*,k] = sqrt(sigma_obj^2*q_obj^2 + psf_sigma[i,k]^2)
   endfor
   
   t0_integrate = systime(/sec)
      
   ;; Integrate the summed convolved MGE....    
   flux_sum = mgeflux(aper_radius, A_invft_mge[*], sigma_array[*], sigma_y_array[*])

   ;; Store MGE info output
   mgeout[i].psf_gauss[0,*] = R_PSF[i,*]
   mgeout[i].psf_gauss[1,*] = psf_sigma[i,*]
   mgeout[i].A0[0:(n_terms_use*n_PSF)-1] = A_invft_mge[*]
   mgeout[i].sigma_x[0:(n_terms_use*n_PSF)-1] = sigma_array[*]
   mgeout[i].sigma_y[0:(n_terms_use*n_PSF)-1] = sigma_y_array[*]
   
   synmag[i,*] = dummy_zpt - 2.5*alog10(flux_Sum)

   tend_integrate = systime(/sec)
   time_integrate[i] = tend_integrate - t0_integrate

   ;if synmag[i,0] LT ModelMag[i] then stop
   
   skip:
   s2 = systime(/sec)
   ;time_left = (s2-s1)*(nmags-1-i)/60 ; Minutes  
   time_left = (s2-s1)/(i+1.0)*(nmags-1-i)/60 ; Minutes   
endfor


end
