;; This file presents some IDL scripts that demonstrate how to use the SYNMAG package.  These scripts can be copy-and-pasted into
;; your IDL terminal.  Please refer to Bundy et al. 2012, AJ, 144, 188, and kindly cite this reference when using the SYNMAG software.
;;
;; Kevin Bundy: Aug 11, 2012

;; --------------------------------------------------
;; The "all-purpose" SYNMAG.PRO
;; --------------------------------------------------
;;   SYNMAG.PRO is designed to be general, so that it can be applied to a variety of datasets.  In addition to this routine,
;;   SDSS_TO_SYNMAG.PRO is specific to the purpose of measuring SYNMAGs using SDSS data.  Specific information about the calling
;;   sequence, arguments, and keywords are provided in the SYNMAG.PRO header.
;;
;; The following is an example of how to use SYNMAG.PRO...  We will apply it to a subset of an SDSS catalog, even though a
;; SYNMAG_SDSS is a better choice when working with SDSS data.

redcol = 255

;; Load the SDSS test catalog (from the Stripe 82 Coadd Photometry)
s82demo = mrdfits('s82demo.fits',1)

wequiv_Dev = where(abs(s82demo.DevMag_R - s82demo.R) LT abs(s82demo.ExpMag_R - s82demo.R), nequiv_Dev) ;; Index on S82demo
wequiv_Exp = where(abs(s82demo.ExpMag_R - s82demo.R) LT abs(s82demo.DevMag_R - s82demo.R), nequiv_Exp) 

;; Define a single-component gaussian PSF and assume it is the same for each object. This PSF will be applied before measuring
;; SYNMAG photometry
PSF_FWHM_constant = 2.1 ;; arcsec

aperture_radii  =[2, 3, 4, 5]/2.0  ;; arcsec, aperture radii for the synthetic aperture photometry

;; Calculate r-band SYNMAGs for the single-component, constant PSF and the SDSS (truncated) Exponential fits
synmag, s82demo.ExpMag_R, s82demo.ExpRad_R, s82demo.ExpAB_R, PSF_FWHM_constant, aperture_radii, synmag_cnst_PSF, n_sersic = 1.0, /SDSS

ww = wequiv_Exp

;; Make some plots
wset, 0 & !p.multi=[0,1,3] & chars=2 & ii=1 & yr = 1.5*[-0.1, 0.1] ;; Compare Exp FIBERMAG versus analytic aperture mags
plot, s82demo[ww].ExpAB_R, (s82demo[ww].fibermag_r - synmag_cnst_PSF[ww,ii]), ps=3, xr=[0,1], yr=yr, xtit='q', chars=chars & oplot, [0,10], [0,0], col=redcol ;; by q
plot, s82demo[ww].ExpRad_R, (s82demo[ww].fibermag_r - synmag_cnst_PSF[ww,ii]), ps=3, xr=[0,4], yr=yr, xtit='Re', chars=chars & oplot, [0,10], [0,0], col=redcol ;; by Radius
plot, s82demo[ww].FracDev_R, (s82demo[ww].fibermag_r - synmag_cnst_PSF[ww,ii]), ps=3, xr=[0,1], yr=yr, xtit='FracDev', chars=chars & oplot, [0,10], [0,0], col=redcol ;; by FracExp
!p.multi=0 & wshow


;; Define a 2-component Gaussian PSF
FWHM1 = 2.07
R_PSF1 = (9/10.)/(2*!pi*(FWHM1/2.3548)^2)
PSF_FWHM_dbl = [ [R_PSF1, FWHM1], [R_PSF1/36, 2*FWHM1]]

;; Calculate r-band synmags for the 2-component PSF and the SDSS (truncated) Exponential fits
;;   Note: The double PSF description could be made more a, introducing 0.01-0.02 mag offsets
synmag, s82demo.ExpMag_R, s82demo.ExpRad_R, s82demo.ExpAB_R, PSF_FWHM_dbl, aperture_radii, synmag_dbl_PSF, n_sersic = 1.0, /SDSS

ww = wequiv_Exp

wset, 0 & !p.multi=[0,1,3] & chars=2 & ii=1 & yr = 1.5*[-0.1, 0.1] ;; Compare Exp FIBERMAG versus analytic aperture mags
plot, s82demo[ww].ExpAB_R, (s82demo[ww].fibermag_r - synmag_dbl_PSF[ww,ii]), ps=3, xr=[0,1], yr=yr, xtit='q', chars=chars & oplot, [0,10], [0,0], col=redcol ;; by q
plot, s82demo[ww].ExpRad_R, (s82demo[ww].fibermag_r - synmag_dbl_PSF[ww,ii]), ps=3, xr=[0,4], yr=yr, xtit='Re', chars=chars & oplot, [0,10], [0,0], col=redcol ;; by Radius
plot, s82demo[ww].FracDev_R, (s82demo[ww].fibermag_r - synmag_dbl_PSF[ww,ii]), ps=3, xr=[0,1], yr=yr, xtit='FracDev', chars=chars & oplot, [0,10], [0,0], col=redcol ;; by FracExp
!p.multi=0 & wshow


;; --------------------------------------------------
;; SYNMAG_SDSS.PRO: 
;; --------------------------------------------------
;; Define a 2-component Gaussian PSF
FWHM1 = 2.07
R_PSF1 = (9/10.)/(2*!pi*(FWHM1/2.3548)^2)
PSF_FWHM_dbl = [ [R_PSF1, FWHM1], [R_PSF1/36, 2*FWHM1]]

;; Calculate synmags, specifically for SDSS
synmag_sdss, s82demo.r, s82demo.DevMag_r, s82demo.DevRad_r, s82demo.DevAB_r, s82demo.ExpMag_r, s82demo.ExpRad_r, s82demo.ExpAB_r, PSF_FWHM_dbl, aperture_radii, synmag_sdss, $
             mgeout = mgeout, hdr_mgeout = hdr_mgeout

ww = indgen(n_elements(s82demo))  ;; Plot all objects
wset, 0 & !p.multi=[0,1,3] & chars=2 & ii=1 & yr = 1.5*[-0.1, 0.1] ;; Compare Exp FIBERMAG versus analytic aperture mags
plot, s82demo[ww].ExpAB_R, (s82demo[ww].fibermag_r - synmag_sdss[ww,ii]), ps=3, xr=[0,1], yr=yr, xtit='q', chars=chars & oplot, [0,10], [0,0], col=redcol ;; by q
plot, s82demo[ww].ExpRad_R, (s82demo[ww].fibermag_r - synmag_sdss[ww,ii]), ps=3, xr=[0,4], yr=yr, xtit='Re', chars=chars & oplot, [0,10], [0,0], col=redcol ;; by Radius
plot, s82demo[ww].FracDev_R, (s82demo[ww].fibermag_r - synmag_sdss[ww,ii]), ps=3, xr=[0,1], yr=yr, xtit='FracDev', chars=chars & oplot, [0,10], [0,0], col=redcol ;; by FracExp
!p.multi=0 & wshow


;; --------------------------------------------------
;; MGESTR.PRO: 
;; --------------------------------------------------
;; Measure the PSF-convolved structural parameters of the output ("observed") profile
;;   1) The major axis half-light radius
;;   2) The circular aperture containing half the light
;;   3) The effective axis ratio
mgestr, Re_mge, q_mge, Re_aper, mge_struct = mgeout  ;; mgeout is delivered by the synmag_sdss call above

ww = wequiv_Exp
!p.multi = [0,3,1] & chars=2
plot, s82demo[ww].ExpRad_r, Re_mge[ww], ps=3, chars=chars, xtit='Intrinsic Re', ytit='Re_MGE' & oplot, [0,100], [0,100], col=redcol
plot, s82demo[ww].ExpRad_r, Re_aper[ww], ps=3, chars=chars, xtit='Intrinsic Re', ytit='Re_Aper'  & oplot, [0,100], [0,100], col=redcol
plot, s82demo[ww].ExpAB_r, q_mge[ww], ps=3, chars=chars, xtit='Intrinsic q', ytit='q_MGE'  & oplot, [0,100], [0,100], col=redcol

ww = wequiv_Dev
!p.multi = [0,3,1] & chars=2
plot, s82demo[ww].DevRad_r, Re_mge[ww], ps=3, chars=chars, xtit='Intrinsic Re', ytit='Re_MGE' & oplot, [0,100], [0,100], col=redcol
plot, s82demo[ww].DevRad_r, Re_aper[ww], ps=3, chars=chars, xtit='Intrinsic Re', ytit='Re_Aper'  & oplot, [0,100], [0,100], col=redcol
plot, s82demo[ww].DevAB_r, q_mge[ww], ps=3, chars=chars, xtit='Intrinsic q', ytit='q_MGE'  & oplot, [0,100], [0,100], col=redcol



