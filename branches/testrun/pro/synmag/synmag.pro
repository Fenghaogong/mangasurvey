;; ----------------------------------------------------------------------------------------------------
;; NAME: SYNMAG.PRO
;;
;; PURPOSE: Procedure predicts aperture magnitudes given information about the PSF and the intrinsic profiles of each object.  See
;; Bundy et al. 2012, AJ, 144, 188 for details.
;;
;; CALLING SEQUENCE: 
;;   synmag, ModelScale, ModelRad, ModelAB, PSF_FWHM_info, aper_radius, synmag, $
;;           N_SERSIC = N_SERSIC, SDSS = SDSS, ZPT_FOR_Ie = ZPT_FOR_Ie, Ie = Ie, $
;;           MGE_intrinsic = MGE_intrinsic, MGEOut = MGEOut
;;
;; ARGUMENTS: 
;;   ModelScale:        Input array of total profile-fit magnitudes or, if ZPT_FOR_Ie, is set, the surface brightness at Re
;;   ModelRad:          Profile model radius, same units as aper_radius
;;   ModelAB:           Profile model axis-ratio, q
;;   PSF_FWHM_info:     Input Gaussian FWHM seeing (arcsec) to be convolved with intrinsic profiles. Can be of the following form:
;;                        1) A scalar, in which case a single Gaussian with this FWHM is convolved with all objects;
;;                        2) A vector with size equal to the number of objects, allowing a different FWHM for each object;

;;                        3) A 2-dimensional array of size [2, n_PSF], [4, n_PSF], or [6, n_PSF] with elements as follows:
;;                             [0, n_PSF] = Gausssian amplitudes for each PSF component
;;                             [1, n_PSF] = Gaussian sigma values if circularly symmetric, or major-axis sigma of each component
;;
;;                           If desired...
;;                             [2, n_PSF] = Gaussian minor-axis sigma values for each component
;;                             [3, n_PSF] = Position Angle offset (degrees counter-clockwise) of the major axis of each PSF component with
;;                                          respect to the major axis of the source profile
;;                           If desired...
;;                             [4, n_PSF] = The x-position offset of each PSF component
;;                             [5, n_PSF] = The y-position offset of each PSF component
;;                             (x,y) offsets are defined in a coordinate frame, centered on the source, with the x-axis oriented
;;                             with the major axis of the source.  
;;
;;                        4) A 3-dimensional array of size [2, n_PSF, nobj], [4, n_PSF, nobj], or [6, n_PSF, nobj], where the
;;                        first two dimensions are as above and the last dimension specifies a different set of PSF Gaussian terms
;;                        for each object.
;;
;;   aper_radius:       Input radii of apertures in arcsec.  All apertures are predicted for all magnitudes
;;   synmag:            Output synthetic aperture magnitude
;;
;; KEYWORDS:
;;   N_SERSIC:        Set to the Sersic-n value of the assumed intrinsic profile shape (if described as a Sersic profile)
;;   SDSS:            Toggle for using the modified SDSS profiles for n=1 (Exponential) or n=4 (de Vaucouleurs)
;;   ZPT_FOR_Ie:      If ModelScale is the surface brightness at Re, this required keyword determines the synmag zeropoint
;;   MGE_INTRINSIC:   Array of size [2, n_terms].  The first element of each term is the full amplitude, Amp, defined such that
;;                    the total 2D flux is Amp*2*!pi*sigma^2.  The second element is the Gaussian sigma of each MGE component.  On
;;                    input, MGE_INTRINSIC defines the intrinsic profile assumed.  On output, it returns the MGE that was used.
;;   MGEOUT           Output structure providing the PSF-convolved MGE information.  For easy input to MGESTR.PRO
;;   Ie               Intrinsic profile's surface brightness as Re   
;;
;; NOTES:
;;   Flux Normalization: If ModelScale is a total magnitude, SYNMAG computes the flux normalization, that is the amplitude of the
;;   fitted profile, by comparing the total flux implied by ModelScale to the 2D integral of the fitted intrinsic profile when its
;;   amplitude is unity.  The integral is taken to infinity which may introduce small biases in the case where ModelScale is
;;   defined by integrating to a finite, but large, radius.  In the case of SDSS, for example, use SYNMAG_SDSS.
;;
;;
;; HISTORY:
;;   Created by K. Bundy on 8/11/2012: Inherited from sdss_to_aper6.pro
;;   Version 2.0 - Add functionality for angle and position offsets between PSF components and profile major axis 9/24/12
;;   Version 2.1 - Implemented consistency of the angle definition to counter-clockwise (12/3/12)
;;   Version 2.2 - Uncommented lines for making output header, and added keyword, hdr_mgeout = hdr_mgeout (2/15/13)
;;   Version 2.3 - Added Ie output and changed I0 to Ie
;;
;;
pro synmag, ModelScale, ModelRad, ModelAB, PSF_FWHM_info, aper_radius, synmag, N_SERSIC = N_SERSIC, $
            SDSS = SDSS, ZPT_FOR_Ie = ZPT_FOR_Ie, Ie = Ie, $
            MGE_intrinsic = MGE_intrinsic, MGEOut = MGEOut, hdr_mgeout = hdr_mgeout

if n_params() LT 6 then begin
   print, '  SYNMAG: Usage:'
   print, '     synmag, ModelScale, ModelRad, ModelAB, PSF_FWHM_info, aper_radius, synmag, $'
   print, '     N_SERSIC = N_SERSIC, SDSS = SDSS, ZPT_FOR_Ie = ZPT_FOR_Ie, Ie = Ie, $'
   print, '     MGE_intrinsic = MGE_intrinsic, MGEOut = MGEOut'
   return
endif



if not(keyword_set(FullRad_Multiplier)) then FullRad_multiplier = 10.0 * ModelRad ;; Radius out to which "total flux" is integrated
if not(keyword_set(MGE_scale_sampling)) then MGE_scale_sampling = 100.0  ;; Scale MGE spatial dimension, defines intrinsic profile's radius value at Re; 1 Re = 100 spatial bins

nmags = n_elements(ModelScale)
napers = n_elements(aper_radius)
if n_elements(aper_radius) GT 1 then synmag = fltarr(nmags, napers)-99 else synmag = fltarr(nmags)-99

;; ..................................................
;; Set up PSF information

PSF_FWHM_array_n_dim_orig = size(PSF_FWHM_info, /n_dim)
if PSF_FWHM_array_n_dim_orig LE 1 then begin  ;; This is a single Gaussian PSF
   if n_elements(PSF_FWHM_info) EQ 1 then PSF_FWHM_array = fltarr(nmags) + PSF_FWHM_info
   n_PSF = 1
   n_params = 2
   psf_smaj = fltarr(nmags) + PSF_FWHM_info / 2.3548  ;; The gaussian sigma of the PSF
   R_PSF = fltarr(nmags) + 1./(2*!pi*psf_smaj^2)      ;; Amplitude of each PSF Gaussian term

   psf_smaj = fltarr(1, nmags) + PSF_FWHM_info / 2.3548  ;; The gaussian sigma of the PSF
   R_PSF = fltarr(1, nmags) + 1./(2*!pi*psf_smaj^2)      ;; Amplitude of each PSF Gaussian term
endif

if PSF_FWHM_array_n_dim_orig EQ 2 then begin ;; This is a Multi-Gaussian PSF, [2, n_PSF] or [4, n_PSF], but the same for all sources. 
   n_PSF = n_elements(PSF_FWHM_info[0,*])
   n_params = n_elements(PSF_FWHM_info[*,0])

   PSF_FWHM_full = reform(fltarr(n_params, n_PSF, nmags), n_params, n_PSF, nmags)
   for k=0,n_PSF-1 do begin                       ;; We need PSF_FWHM_full -> [n_params, n_PSF, nmags]
      for l=0,n_params-1 do begin
         replicate_inplace, PSF_FWHM_full, PSF_FWHM_info[l,k], 3, [l,k,0]
      endfor
   endfor

;;    PSF_FWHM_full = fltarr(nmags, n_params, n_PSF)
;;    for k=0,n_PSF-1 do begin                       ;; We need PSF_FWHM_full -> [nmags, n_params, n_PSF]
;;       replicate_inplace, PSF_FWHM_full, PSF_FWHM_info[0,k], 1, [0,0,k]
;;       replicate_inplace, PSF_FWHM_full, PSF_FWHM_info[1,k], 1, [0,1,k]
;;    endfor

;;    R_PSF = reform(PSF_FWHM_array[*,0,*])            ;; Amplitude of each PSF Gaussian term
;;    psf_smaj = reform(PSF_FWHM_array[*,1,*])/2.3548 ;; The gaussian sigma of the PSF

   PSF_FWHM_array = PSF_FWHM_full  ;; Replace with full array

   R_PSF = reform(PSF_FWHM_array[0,*,*])            ;; Amplitude of each PSF Gaussian term
   psf_smaj = reform(PSF_FWHM_array[1,*,*])/2.3548 ;; The gaussian sigma of the PSF (major-axis if four parameters are given)

   ;; if n_params EQ 4 then begin
   ;;    psf_smin = reform(PSF_FWHM_array[2,*,*])/2.3548
   ;;    PA_offset = reform(PSF_FWHM_array[3,*,*])
   ;; endif else begin
   ;;    psf_smin = psf_smaj
   ;;    PA_offset = 0
   ;; endelse

   case n_params OF
      2: begin
         psf_smin = psf_smaj
         PA_offset = 0
      end
      4: begin
         psf_smin = reform(PSF_FWHM_array[2,*,*])/2.3548
         PA_offset = reform(PSF_FWHM_array[3,*,*])
         PSF_X0_input = 0 * psf_smin
         PSF_Y0_input = 0 * psf_smin
      end
      6: begin
         psf_smin = reform(PSF_FWHM_array[2,*,*])/2.3548
         PA_offset = reform(PSF_FWHM_array[3,*,*])
         PSF_X0_input = reform(PSF_FWHM_array[4,*,*])
         PSF_Y0_input = reform(PSF_FWHM_array[5,*,*])
      end      
   endcase

   ;; Spatial coefficients of the fully-general 2D Gaussian of the form: AA * exp( -(a(x-x0)^2 + 2b(x-x0)(y-y0) + c(y-y0)^2 )
   a_psf = (cos(pa_offset*!dtor)^2.0)/(2.0*psf_smaj^2.0) + (sin(pa_offset*!dtor)^2.0)/(2.0*psf_smin^2.0)
   c_psf = (sin(pa_offset*!dtor)^2.0)/(2.0*psf_smaj^2.0) + (cos(pa_offset*!dtor)^2.0)/(2.0*psf_smin^2.0)
   b_psf = sin(pa_offset*!dtor)*cos(pa_offset*!dtor)*( 1.0/(2.0*psf_smin^2) - 1.0/(2.0*psf_smaj^2) )      ;; (1/(2smin^2) - 1/(2smaj^2)) --> defines clockwise PA of major-axis from x-axis

endif


if PSF_FWHM_array_n_dim_orig EQ 3 then begin     ;; We were given a full Multi-Gaussian PSF array, defined for all sources (with dimensions, [n_params, n_PSF, nmags])
   n_PSF = n_elements(PSF_FWHM_info[0,*,0])
   n_params = n_elements(PSF_FWHM_info[*,0,0])
   
   PSF_FWHM_array = PSF_FWHM_info  ;; Replace with full array

   R_PSF = reform(PSF_FWHM_array[0,*,*])            ;; Amplitude of each PSF Gaussian term
   psf_smaj = reform(PSF_FWHM_array[1,*,*])/2.3548 ;; The gaussian sigma of the PSF (major-axis if four parameters are given)

   ;; if n_params EQ 4 then begin
   ;;    PSF_smin = reform(PSF_FWHM_array[2,*,*])/2.3548  ;; Minor axis sigma of the PSF
   ;;    PA_offset = reform(PSF_FWHM_array[3,*,*])           ;; Position angle offset, counterclockwise in degrees, between the profile major axis and PSF component major axis
   ;; endif  else begin
   ;;    psf_smin = psf_smaj
   ;;    PA_offset = 0
   ;; endelse

   case n_params OF
      2: begin
         psf_smin = psf_smaj
         PA_offset = 0
      end
      4: begin
         psf_smin = reform(PSF_FWHM_array[2,*,*])/2.3548
         PA_offset = reform(PSF_FWHM_array[3,*,*])
         PSF_X0_input = 0 * psf_smin
         PSF_Y0_input = 0 * psf_smin
      end
      6: begin
         psf_smin = reform(PSF_FWHM_array[2,*,*])/2.3548
         PA_offset = reform(PSF_FWHM_array[3,*,*])
         PSF_X0_input = reform(PSF_FWHM_array[4,*,*])
         PSF_Y0_input = reform(PSF_FWHM_array[5,*,*])
      end      
   endcase

   
   ;; Spatial coefficients of the fully-general 2D Gaussian of the form: AA * exp( -(a(x-x0)^2 + 2b(x-x0)(y-y0) + c(y-y0)^2 )
   a_psf = (cos(pa_offset*!dtor)^2.0)/(2.0*psf_smaj^2.0) + (sin(pa_offset*!dtor)^2.0)/(2.0*psf_smin^2.0)
   c_psf = (sin(pa_offset*!dtor)^2.0)/(2.0*psf_smaj^2.0) + (cos(pa_offset*!dtor)^2.0)/(2.0*psf_smin^2.0)
   ;b_psf = sin(pa_offset*!dtor)*cos(pa_offset*!dtor)*( 1.0/(2.0*psf_smin^2) - 1.0/(2.0*psf_smaj^2) )      ;; (1/(2smin^2) - 1/(2smaj^2)) --> defines clockwise PA of major-axis from x-axis
   b_psf = sin(pa_offset*!dtor)*cos(pa_offset*!dtor)*( 1.0/(2.0*psf_smaj^2) - 1.0/(2.0*psf_smin^2) )      ;; (1/(2smaj^2) - 1/(2smin^2)) --> defines counter-clockwise PA of major-axis from x-axis

endif


;; ..................................................
;; Determine MGE approximation to the intrinsic profile (if needed)
if n_elements(MGE_intrinsic) LT 2 then begin  ;; No MGE for the intrinsic profile was passed

   if n_elements(n_sersic) NE 1 then begin
      print
      print, ' SYNMAG: No MGE Passed.  To use a Sersic intrinsic profile, please set N_SERSIC keyword to a single value.'
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
if keyword_set(ZPT_FOR_Ie) then dummy_zpt = ZPT_FOR_Ie else begin
   ModelMag = ModelScale
   dummy_zpt = 22.5
   total_flux = 10^(0.4*(dummy_zpt-ModelMag))
endelse

Ie = fltarr(nmags)-1

;; Structure for holding output convolved MGE information
n_terms_mge = n_elements(mge[0,*])

mgeout = {psf_gauss:fltarr(n_params, n_psf)-1, $
          A0:fltarr(n_terms_mge*n_psf)-99, sigma_x:fltarr(n_terms_mge*n_psf)-99, sigma_y:fltarr(n_terms_mge*n_psf)-99}
mgeout = replicate(mgeout, nmags)

fxbhmake, hdr_mgeout, nmags, /initialize, /date
sxaddpar, hdr_mgeout, 'LEN_UNIT', 'Arcsec', 'Units of sigma values'
sxaddpar, hdr_mgeout, 'COEF', 'Amplitude', 'Type of MGE normalization coefficient'
;; sxaddpar, hdr_mgeout, 'Pref', '1=Exp, 2=Dev', 'Type of MGE normalization coefficient'
;; sxaddpar, hdr_mgeout, 'N_EXP', n_terms_exp, 'Number of terms in the Exponential MGE'
;; sxaddpar, hdr_mgeout, 'N_DEV', n_terms_dev, 'Number of terms in the deVaucouleurs MGE'

time_integrate = fltarr(nmags)  ;; Time it takes to do MGE flux integrals
time_mock = fltarr(nmags)       ;; Time it takes to do mock image with filtering

A_invft_mge = fltarr(n_terms_mge, n_PSF)  ;; Output PSF-convolved MGE amplitudes
sigma_array = fltarr(n_terms_mge, n_PSF)  ;; Output PSF-convolved major-axis sigmas
sigma_y_array =fltarr(n_terms_mge, n_PSF) ;; Output PSF-convolved minor-axis sigmas

delta = fltarr(n_terms_mge, n_PSF)
a_out = fltarr(n_terms_mge, n_PSF)
b_out = fltarr(n_terms_mge, n_PSF)
c_out = fltarr(n_terms_mge, n_PSF)

;; We prepare X/Y offset arrays at this point -- this is reformating them from one offset per PSF term to an offset for each PSF term crossed with each MGE term
if n_params GE 4 then begin
   psf_x0 = fltarr(n_terms_mge, n_PSF, nmags)
   psf_y0 = fltarr(n_terms_mge, n_PSF, nmags)

   for j=0, n_terms_mge-1 do begin
      psf_x0[j,*,*] = psf_x0_input
      psf_y0[j,*,*] = psf_y0_input
   endfor
endif

;; ----------------------------------------------------------------------------------------------------
;; Loop over all objects 
for i=0L, nmags-1 do begin  
   s1 = systime(/sec)


   if i GT 0 AND (i - 100*(i/100) EQ 0 OR i EQ nmags-1) then begin
      print, format='($, A)', strjoin(strarr(100)+string(8B))  ;; Erase 100 characters of previous line
      print, format='($, A)', ' SYNMAG: Number '+strc(i)+' of '+strc(nmags-1)+'.  Time left: '+dplace(time_left,1)+' min'
   endif
   if badflag[i] GE 1 then goto, skip

   ;; ..................................................
   ;; Determine normalization of the flux model
   if keyword_set(ZPT_FOR_Ie) then Ie[i] = ModelScale[i] else begin
      ;; Figure out the flux normalization from the total (model) magnitude provided
            
      sigma_y_mge = ModelAB[i] * sigma_mge 

      flux_norm = mgeflux(-1, amp_mge, ModelRad[i]*sigma_mge, ModelRad[i]*sigma_y_mge, /precise)

      ;flux_norm = mgeflux(4.0*ModelRad[i], 1.0, ModelRad[i], qprof=ModelAB[i], /precise, /Exp_SDSS)

      Ie[i] = total_flux[i] / flux_norm 
   endelse

   ;; ----------------------------------------------------------------------------------------------------
   ;; Use the MGE expansion to analytically describe the convolved profile
   mgeout[i].psf_gauss[0,*] = R_PSF[*,i]
   mgeout[i].psf_gauss[1,*] = psf_smaj[*,i]

   if n_params GE 4 then begin      
      mgeout[i].psf_gauss[2,*] = psf_smin[*,i]
      mgeout[i].psf_gauss[3,*] = pa_offset[*,i]
   endif

   if n_params EQ 6 then begin
      mgeout[i].psf_gauss[4,*] = psf_x0_input[*,i]
      mgeout[i].psf_gauss[5,*] = psf_y0_input[*,i]
   endif

   q_obj = ModelAB[i]

   sigma_obj = sigma_mge * (ModelRad[i])  ;; Object's instrinsic MGE gaussian sigmas for object in question (arcsec)
   ;A_mge = Ie * (amp_mge / sqrt(2*!pi*sigma_mge^2))            ;; Object's intrinsic MGE gaussian coeffecients of the expansion
   A_mge = Ie[i] * amp_mge             ;; Object's intrinsic MGE gaussian coeffecients of the expansion

   if n_params EQ 2 then begin
      ;; Loop over the Gaussian components of the PSF
      for k=0, n_PSF-1 do begin
         A_invft_mge[*,k] = (A_mge * R_PSF[k,i] * 2*!pi * sigma_obj^2 * q_obj * psf_smaj[k,i]^2) / $
                            sqrt( (sigma_obj^2 + psf_smaj[k,i]^2)*(sigma_obj^2*q_obj^2 + psf_smaj[k,i]^2) ) ;; Array of Prefactors (on convolved gaussian term)
         
         sigma_array[*,k] = sqrt(sigma_obj^2 + psf_smaj[k,i]^2)
         sigma_y_array[*,k] = sqrt(sigma_obj^2*q_obj^2 + psf_smaj[k,i]^2)
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
   endif

   if n_params GE 4 then begin

      ;; Loop over the Gaussian components of the PSF
      for k=0, n_PSF-1 do begin
         ;; Useful PSF-convolved output constant: delta = (sigx^2 * sigy^2)/(2*smaj^2*smin^2) + a*sigx^2 + c*sigy^2 + 0.5
         delta[*,k] = (sigma_obj^2.0 * (sigma_obj^2.0 * q_obj^2))/(2.0*psf_smaj[k,i]^2.0*psf_smin[k,i]^2.0) + a_psf[k,i]*sigma_obj^2 + c_psf[k,i]*(sigma_obj*q_obj)^2 + 0.5

         ;; Spatial coefficients of the PSF-convolved profile
         a_out[*,k] = (1.0/(4.0*delta[*,k]))*( (sigma_obj^2.0*q_obj^2)/(psf_smaj[k,i]^2.0*psf_smin[k,i]^2.0) + 2*a_psf[k,i] )  ;; so-called, alpha
         c_out[*,k] = (1.0/(4.0*delta[*,k]))*( (sigma_obj^2.0)/(psf_smaj[k,i]^2.0*psf_smin[k,i]^2.0) + 2*c_psf[k,i] )          ;; so-called, gamma
         b_out[*,k] = b_psf[k,i] / (2.0*delta[*,k])                                                                            ;; so-called, beta

         A_invft_mge[*,k] = A_mge * R_PSF[k,i] * sigma_obj^2 * q_obj * !pi * sqrt(2.0/delta[*,k])  ;; Amplitudes

      endfor
      
      ;; Determine convolved output's position angle: 
      ;;   Clockwise: tan(2*theta) = 2*beta / (gamma - alpha)      
      ;;   Counter-clockwise:  tan(2*theta) = 2*beta / (alpha - gamma)  
      ;; pa_out = !radeg*atan(2.0*b_out, (c_out-a_out))/2.0  ;; Degrees clockwise, Note the use of 2 arguments in atan(y,x) to avoid undetermined angle degeneracies
      pa_out = !radeg*atan(2.0*b_out, (a_out-c_out))/2.0  ;; Degrees counter-clockwise, Note the use of 2 arguments in atan(y,x) to avoid undetermined angle degeneracies
      
      ;; Determine convolved output's major- and minor-axis Gaussian widths
      ;;   2 smaj_out^2 = ( (alpha + gamma) + sqrt((alpha+gamma)2 - 4*(alpha*gamma - beta^2)) ) / (2*(alpha*gamma - beta^2))
      ;;   2 smaj_out^2 = ( (alpha + gamma) - sqrt((alpha+gamma)2 - 4*(alpha*gamma - beta^2)) ) / (2*(alpha*gamma - beta^2))
      smaj_out = sqrt(0.5*(((a_out + c_out) + sqrt( (a_out+c_out)^2.0 - 4.0*(a_out*c_out - b_out^2) ))/(2.0*(a_out*c_out - b_out^2))))
      smin_out = sqrt(0.5*(((a_out + c_out) - sqrt( (a_out+c_out)^2.0 - 4.0*(a_out*c_out - b_out^2) ))/(2.0*(a_out*c_out - b_out^2))))

      t0_integrate = systime(/sec)
      
      ;; Integrate the summed convolved MGE numerically....    
      flux_sum = mgeflux(aper_radius, A_invft_mge[*], smaj_out[*], smin_out[*], xc = -(psf_x0[*,*,i])[*], yc = -(psf_y0[*,*,i])[*], PA = pa_out[*])  ;; 
      
      ;; Store MGE info output
      mgeout[i].A0 = A_invft_mge[*]
      mgeout[i].sigma_x = smaj_out[*]
      mgeout[i].sigma_y = smin_out[*]
      
      synmag[i,*] = dummy_zpt - 2.5*alog10(flux_Sum)
      
      tend_integrate = systime(/sec)
      time_integrate[i] = tend_integrate - t0_integrate
   endif

      
   skip:
   s2 = systime(/sec)
   time_left = (s2-s1)*(nmags-1-i)/60 ; Minutes   
endfor
print

end
