;+
; function mlskysubtract
;
; This is the high-level sky subtraction routine, analagous to what is
; done in BOSS extract_object.pro.  It calls mlcalcsky to actually
; figure out the sky values.;
; Returns 0 if everything ok, returns an error code if there was a problem.
;
; Author:
;   David R. Law (Dunlap Institute; drlaw@di.utoronto.ca)
;   First written: 01/28/2013
;   Last modified: 02/27/2013
;
; REVISION HISTORY:
;   v1.0: 28-Jan-2013  D. Law
;       Imported code from BOSS extract_object.pro, start integrating
;       with MaNGA algorithms and calling
;   v1.1: 05-Feb-2013  D. Law
;       Revised call to require output filename
;   v1.2: 15-Feb-2013  D. Law
;       Revised to require keyword 'wave', and handle 2'', 3'', 5''
;       fibers seperately.  Outputs extra extensions to spSFrame
;       where everything is rectified to the same wave scale.
;   v1.3: 27-Feb-2013 D. Law
;       Added ximg and superflat to output extensions again.  Now
;       mimics spFrame except with slitmap in place of plugmap.
;   v2: 23-Jun-2013 D. Law
;       Revised to include possibility of using edge fibers for sky
;       fibers, do multi-pass subtraction for continuum and masked
;       line regions, create output plot quality files.
;       Create output poisson comparison extension.
;-

function mlskysubtract, spframefile, spsframefile, obsparam, slitmap, camera, wave, visual=visual, edgesky=edgesky, plotfile=plotfile, flatfile=flatfile
; NB flatfile only matters for poisson calculation

;;;;;;;;;;;;;;;;;;;;;;;;;;

; Read in files
mlframe_read,spframefile,objflux=flux,objivar=fluxivar,wset=vacset,mask=pixelmask, $
   dispset=dispset,ximg=ximg,superflat=superflat,hdr=objhdr

; Read in flatfield, if specified
if (keyword_set(flatfile)) then flat=mrdfits(flatfile,0)

; Set up plot file if specified, otherwise use display
if (keyword_set(plotfile)) then begin
  set_plot,'ps'
  device,filename=plotfile,/color
  splog,'Printing to ',plotfile
endif else set_plot,'x'


nx = (size(flux,/dim))[0] 
ny = (size(flux,/dim))[1] 

; Set up the blue/red frame parameters
if (camera eq 'b1') then nbkpt = 3*nx/4
if (camera eq 'r1') then nbkpt = nx

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; Define which fibers are sky fibers of 2, 3, 5'' flavor

ftype2=where((slitmap.fsize eq 2.))
ftype3=where((slitmap.fsize eq 3.))
ftype5=where((slitmap.fsize eq 5.))

; Split up structures accordingly
nbkpt2=nbkpt
nbkpt3=nbkpt
nbkpt5=nbkpt

flux2=flux[*,ftype2]
flux3=flux[*,ftype3]
flux5=flux[*,ftype5]

ny2 = (size(flux2,/dim))[1] 
ny3 = (size(flux3,/dim))[1] 
ny5 = (size(flux5,/dim))[1] 

fluxivar2=fluxivar[*,ftype2]
fluxivar3=fluxivar[*,ftype3]
fluxivar5=fluxivar[*,ftype5]

pixelmask2=pixelmask[*,ftype2]
pixelmask3=pixelmask[*,ftype3]
pixelmask5=pixelmask[*,ftype5]

traceset2xy,vacset,junk,loglam
loglam2=loglam[*,ftype2]
loglam3=loglam[*,ftype3]
loglam5=loglam[*,ftype5]

vacset2=traceset_trim(vacset,ftype2)
vacset3=traceset_trim(vacset,ftype3)
vacset5=traceset_trim(vacset,ftype5)

traceset2xy,dispset,junk,dispval
dispval2=dispval[*,ftype2]
dispval3=dispval[*,ftype3]
dispval5=dispval[*,ftype5]

slitmap2=slitmap[ftype2]
slitmap3=slitmap[ftype3]
slitmap5=slitmap[ftype5]

superflat2=superflat[*,ftype2]

; Define where the skies are within each subgroup
if (keyword_set(edgesky)) then begin
  ; All block edges
  edgefibid=[1,30,31,49,65,101,132,161,162,180,182,210,211,241,242,260,261,290,321,350,351,380,381,399,400,429,445,481,482,500,501,530,531,560]

  ; Block 'edges' appropriate for 6653 masking out object flux
  ;edgefibid=[1,30,47,65,101,132,159,164,173,182,210,214,241,254,259,261,290,351,380,400,429,445,481,531,560]
  ; Block 'edges' for 6653, 3 per block
  ;edgefibid=[1,17,30,47,65,80,101,132,144,159,164,173,182,198,210,214,228,241,254,259,261,276,290,351,365,380,400,429,445,464,481,531,546,560]

  ; Match to find where the edge fiber indices are, and how many
  ; Match across entire slit
  match,slitmap.fiberid,edgefibid,iskies,junk,count=nskies
  ; Match across subset of 2'' fibers only
  match,slitmap2.fiberid,edgefibid,iskies2,junk,count=nskies2
endif else begin
  ; Match across entire slit
  iskies=where((slitmap.ifuname eq 'SKY2') AND (slitmap.plugstatus eq 1), nskies)
  ; Match across subset of 2'' fibers only
  iskies2=where((slitmap2.ifuname eq 'SKY2') AND (slitmap2.plugstatus eq 1), nskies2)
endelse

iskies3=where((slitmap3.ifuname eq 'SKY3') AND (slitmap3.plugstatus eq 1), nskies3)
iskies5=where((slitmap5.ifuname eq 'SKY5') AND (slitmap5.plugstatus eq 1), nskies5)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; First pass
tai_mid=double(obsparam.taistart)+obsparam.exptime/2.

; Do 2'' fibers
skystruct2 = mlcalcsky(flux2, fluxivar2, loglam2, slitmap2, $
  skysub2, skysubivar2, iskies=iskies2, pixelmask=pixelmask2, $
  fibermask=fibermask2, upper=3.0, lower=3.0, tai=tai_mid, nbkpt=nbkpt2)

if (NOT keyword_set(skystruct2)) then begin
  splog,'Problem with skysubtract2- quit!'
  mlquitmanga3d,-200L 
endif

; Do 3'' fibers
skystruct3 = mlcalcsky(flux3, fluxivar3, loglam3, slitmap3, $
  skysub3, skysubivar3, iskies=iskies3, pixelmask=pixelmask3, $
  fibermask=fibermask3, upper=3.0, lower=3.0, tai=tai_mid, nbkpt=nbkpt3)

if (NOT keyword_set(skystruct3)) then begin
  splog,'Problem with skysubtract3- quit!'
  mlquitmanga3d,-200L 
endif

; Do 5'' fibers
skystruct5 = mlcalcsky(flux5, fluxivar5, loglam5, slitmap5, $
  skysub5, skysubivar5, iskies=iskies5, pixelmask=pixelmask5, $
  fibermask=fibermask5, upper=3.0, lower=3.0, tai=tai_mid, nbkpt=nbkpt5)

if (NOT keyword_set(skystruct5)) then begin
  splog,'Problem with skysubtract5- quit!'
  mlquitmanga3d,-200L 
endif

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; If any of the sky-fibers are bad, then re-do sky-subtraction.

; Do 2'' fibers
ibadfib2 = where(djs_median(skysub2[*,iskies2]^2 * $
  skysubivar2[*,iskies2], 1) GT 2.0)               
if (ibadfib2[0] NE -1) then begin
  fibermask2[iskies2[ibadfib2]] = fibermask2[iskies2[ibadfib2]] OR fibermask_bits('BADSKYFIBER')

  splog, 'Calling skysubtract2 again; masked skyfibers',string(iskies2[ibadfib2])
  skystruct2 = mlcalcsky(flux2, fluxivar2, loglam2, slitmap2, $
    skysub2, skysubivar2, iskies=iskies2, pixelmask=pixelmask2, $
     fibermask=fibermask2, upper=10.0, lower=10.0, tai=tai_mid, nbkpt=nbkpt2)

  if (NOT keyword_set(skystruct2)) then begin
   splog,'Problem with skysubtract2- quit!'
   mlquitmanga3d,-200L 
  endif
endif

; Do 3'' fibers
ibadfib3 = where(djs_median(skysub3[*,iskies3]^2 * $
  skysubivar3[*,iskies3], 1) GT 2.0)               
if (ibadfib3[0] NE -1) then begin
  fibermask3[iskies3[ibadfib3]] = fibermask3[iskies3[ibadfib3]] OR fibermask_bits('BADSKYFIBER')

  splog, 'Calling skysubtract3 again; masked skyfibers',string(iskies3[ibadfib3])
  skystruct3 = mlcalcsky(flux3, fluxivar3, loglam3, slitmap3, $
    skysub3, skysubivar3, iskies=iskies3, pixelmask=pixelmask3, $
     fibermask=fibermask3, upper=10.0, lower=10.0, tai=tai_mid, nbkpt=nbkpt3)

  if (NOT keyword_set(skystruct3)) then begin
   splog,'Problem with skysubtract3- quit!'
   mlquitmanga3d,-200L 
  endif
endif

; Do 5'' fibers
ibadfib5 = where(djs_median(skysub5[*,iskies5]^2 * $
  skysubivar5[*,iskies5], 1) GT 2.0)               
if (ibadfib5[0] NE -1) then begin
  fibermask5[iskies5[ibadfib5]] = fibermask5[iskies5[ibadfib5]] OR fibermask_bits('BADSKYFIBER')

  splog, 'Calling skysubtract5 again; masked skyfibers',string(iskies5[ibadfib5])
  skystruct5 = mlcalcsky(flux5, fluxivar5, loglam5, slitmap5, $
    skysub5, skysubivar5, iskies=iskies5, pixelmask=pixelmask5, $
     fibermask=fibermask5, upper=10.0, lower=10.0, tai=tai_mid, nbkpt=nbkpt5)

  if (NOT keyword_set(skystruct5)) then begin
   splog,'Problem with skysubtract5- quit!'
   mlquitmanga3d,-200L 
  endif
endif

; QA plots for chi^2 from 1D sky-subtraction.
;if keyword_set(VISUAL) then getwindow,/open
qaplot_skysub, flux2, fluxivar2, skysub2, skysubivar2, vacset2, iskies2, title=' 1D Sky-subtraction'
;qaplot_skysub, flux3, fluxivar3, skysub3, skysubivar3, vacset3, iskies3, title=' 1D Sky-subtraction'
;qaplot_skysub, flux5, fluxivar5, skysub5, skysubivar5, vacset5, iskies5, title=' 1D Sky-subtraction'

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;----------
; Sky-subtract one final time, this time with dispset (PSF subtraction)
; (rejected sky fibers from above remain rejected).
; Modify pixelmask in this call.
; Sets skysub2b and skysubivar2b
nskypoly = 2L
skystruct = mlcalcsky(flux2, fluxivar2, loglam2, slitmap2, $
  skysub2b, skysubivar2b, iskies=iskies2, pixelmask=pixelmask2, $
  fibermask=fibermask2, upper=10.0, lower=10.0, tai=tai_mid, $
  ; dispset=dispset, $ ; Why is this commented out??
  dispval=dispval2, $
  npoly=nskypoly, nbkpt=nbkpt2, $
  relchi2set=relchi2set, newmask=newmask)
pixelmask = newmask

;if (NOT keyword_set(skystruct)) then begin
;  splog,'Problem with skysubtract- quit!'
;  mlquitmanga3d,-200L 
;endif

; QA plots for chi^2 from 2D sky-subtraction.
;if keyword_set(VISUAL) then getwindow,/open
;qaplot_skysub, flux, fluxivar, skysub, skysubivar, $
;  vacset, iskies, title=' 2D Sky-subtraction'

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; DRL second pass sky subtraction for 2'' fibers
; Subtract off low-order continuum fit to all sky fibers within a block
; This sets skysub2c and skysubivar2c (same as 2b by default)
skysub2c=skysub2b
skysubivar2c=skysubivar2b

; Loop through the v-groove blocks
firstblock=min(slitmap2.blockid)
lastblock=max(slitmap2.blockid)
for i=firstblock,lastblock do begin
  ; Identify the fibers in this block
  dofibers=where(slitmap2.blockid eq i,ndo)
  ; Match where the fibers in block are sky fibers
  match,dofibers,iskies2,temp1,temp2,count=nsky

  ; If there was at least 1 sky fiber, use them to do subtraction
  if (nsky ge 1) then begin
    ; Identify the sky fibers in this block
    skyfibers=dofibers[temp1]
    ; And assign their flux, ivar, and wave solutions to temporary arrays
    tempflux=skysub2b[*,skyfibers]
    tempivar=skysubivar2b[*,skyfibers]
    tempwave=loglam2[*,skyfibers]

    ; Sort according to wavelength
    isort=sort(tempwave)
    tempwave=tempwave[isort]
    tempflux=tempflux[isort]
    tempivar=tempivar[isort]

    ; Bspline fit to the sky spectra
    ; Use lowest order with minimal number of breakpoints
    ; to get a straight line fit
    ;everyn=3000
    ;nord=2
    ; Use everyn=500*nsky, nord=3 to get a good-looking fit in r?
    nord=3
    everyn=500*nsky
    groupsize=nsky
    bkpt=0
    sset=bspline_iterfit(tempwave,double(tempflux),invvar=double(tempivar),nord=nord,everyn=everyn,bkpt=bkpt,maxrej=0,yfit=skyfit) 

    tempskyfit=interpol(skyfit,tempwave,loglam2[*,dofibers])
    skysub2c[*,dofibers]=skysub2b[*,dofibers]-tempskyfit

    ; Quality control plot
    ; Set display yrange to about 5sigma at the midpoint
    temp=mlmeanclip(tempflux[nsky*nx/2.:nsky*nx/2.+nx/4.],ymean,ysig)
    plot,tempflux,yrange=[-5*ysig,5*ysig]
    oplot,skyfit,color=250
    xyouts, nsky*nx/2.,4*ysig,strcompress('Block #: ' +string(i)),color=250
  endif
endfor

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; DRL third pass sky subtraction for 2'' fibers
; Subtract off residuals of edge fibers for wavelengths where
; subtraction was bad.
; This sets skysub2d and skysubivar2d.

skysub2d=skysub2c
skysubivar2d=skysubivar2c

; Set typical inverse gain
; Ideally will calculate somewhere, for now just use what I estimate
; are typical values.  Multiply fluxes by this number to get e- units.
if (camera eq 'b1') then igain=0.913 ; e-/ADU
if (camera eq 'r1') then igain=0.827 ; e-/ADU

; Calculate ratio of noise to poisson
nphot_skysub=skysub2c*superflat2*igain
nphot=flux2*superflat2*igain
pratio=nphot_skysub/sqrt(nphot)
; Make sure no NaN or Inf values, zero them out
pratio[where(finite(pratio) eq 0)]=0.

; Interpolate all sky fiber vectors to a common wavelength grid
; (pick the wavelength solution of the first sky fiber, as good as any)
skyblock=pratio[*,iskies2]
ivarblock=skysubivar2c[*,iskies2]
for i=1,nskies2-1 do begin
  skyblock[*,i]=interpol(skyblock[*,i],loglam2[*,iskies2[i]],loglam2[*,iskies2[0]])
  ivarblock[*,i]=interpol(ivarblock[*,i],loglam2[*,iskies2[i]],loglam2[*,iskies2[0]])
endfor

; At each wavelength, compare sky fibers to determine deviations
; from Poisson across the slit.  It's only the wavelengths where
; rms across blocks isn't consistent with poisson that we want
; to do another subtraction.
; NB: This isn't *quite* right, since putting things on the same
; wave grid introduces correlations, but close enough for
; the purpose of flagging bad lines.
pvector=fltarr(nx)
qvector=fltarr(nx)
ivarvector=fltarr(nx)
for i=0,nx-1 do begin
  a=mlmeanclip(skyblock[i,*],a1,a2,clipsig=5.)
  pvector[i]=a2
  qvector[i]=a1
  b=mlmeanclip(ivarblock[i,*],b1,b2,clipsig=5.)
  ivarvector[i]=b1
endfor

; Fit the overall trend of pvector with a low-order bspline
everyn=300
bkpt=0

sset=bspline_iterfit(loglam2[*,iskies2[0]],double(pvector),invvar=double(ivarvector),nord=3,everyn=everyn,bkpt=bkpt,maxrej=0,yfit=pfit) 
; Subtract off the bspline fit
pvector=pvector-pfit
; Zero out end caps
pvector[where(ivarvector eq 0.)]=0.
; Figure out typical rms of vector
temp=mlmeanclip(pvector[where(ivarvector ne 0.)],temp1,temp2)
; Identify where peaks are n-sigma above zero
; This is where skyline subtraction wasn't great between blocks

nsiglinefail=5.
linefail=where(abs(pvector) gt nsiglinefail*temp2,nfail)
if (nfail ne 0) then begin
  ; Define a box 1 pixel either side
  linefailmin=linefail-1
  linefailmax=linefail+1
  ; Ensure these don't go beyond ends of spectrum
  temp=where(linefailmin lt 0,ntemp)
  if (ntemp ne 0) then linefailmin[temp]=0
  temp=where(linefailmax ge nx,ntemp)
  if (ntemp ne 0) then linefailmax[temp]=nx-1
endif

; Define a skyline mask in 1d vector form.  This identifies
; where there was high rms in residuals.
skylinevec=intarr(nx)
if (nfail ne 0) then begin
  for i=0,nfail-1 do skylinevec[linefailmin[i]:linefailmax[i]]=1
endif
; Kludge for b1 camera: don't flag very ends of wave regime
if (camera eq 'b1') then skylinevec[0:1000]=0
if (camera eq 'b1') then skylinevec[3000:nx-1]=0

; Now also identify where there were high mean values in residuals.
a=mlmeanclip(qvector[where(qvector ne 0.)],a1,a2)
qvector=(qvector-a1)/a2
nsiglinefail=3
linefail=where(abs(qvector) gt nsiglinefail,nfail)
failvec=intarr(nx)
failvec[linefail]=1
; If flagged pixel i and i+2, also flag pixel i+1
for i=0,nfail-1 do if ((failvec[linefail[i]] eq 1)and(failvec[linefail[i]+2] eq 1)) then failvec[linefail[i]+1]=1
linefail=where(failvec eq 1,nfail)
skylinevecb=intarr(nx)
; Look for where 3 pixels flagged in a row
for i=0,nfail-1 do begin
  if ((failvec[linefail[i]] eq 1)and(failvec[linefail[i]+1] eq 1)and(failvec[linefail[i]+2] eq 1)) then begin
    skylinevecb[linefail[i]:linefail[i]+2]=1
  endif
endfor
; Kludge for b1 camera: don't flag very ends of wave regime
if (camera eq 'b1') then skylinevecb[0:1000]=0
if (camera eq 'b1') then skylinevecb[3000:nx-1]=0

; Combine the two flagging methods
overlap=where(skylinevecb eq 1,noverlap)
if (noverlap ne 0) then skylinevec[overlap]=1
; Interpolate the mask to the wavelength solution of every fiber
skylinemask2=interpol(skylinevec,loglam2[*,iskies2[0]],loglam2[*,0:ny2-1])
skylinemask2[where(skylinemask2 ne 1)]=0

; Loop through the v-groove blocks to do the resubtraction
for i=firstblock,lastblock do begin
  ; Identify the fibers in this block
  dofibers=where(slitmap2.blockid eq i,ndo)
  ; Match where the fibers in block are sky fibers
  match,dofibers,iskies2,temp1,temp2,count=nsky

  ; If there was at least 1 sky fiber, use them to do subtraction
  if (nsky ge 1) then begin
    ; Identify the sky fibers in this block
    skyfibers=dofibers[temp1]
    ; And assign their flux, ivar, and wave solutions to temporary arrays
    tempflux=skysub2c[*,skyfibers]
    tempivar=skysubivar2c[*,skyfibers]
    tempwave=loglam2[*,skyfibers]

    ; If only 1 sky, then it IS the reference fit vector
    if (nsky eq 1) then begin
      skyfit=tempflux
    ; If more than 1 sky, bspline them together
    endif else begin
      ; Sort according to wavelength
      isort=sort(tempwave)
      tempwave=tempwave[isort]
      tempflux=tempflux[isort]
      tempivar=tempivar[isort]
      ; Kill any bad values using ivar as a mask
      tempflux[where(tempivar eq 0.)]=0.

      nord=2
      everyn=1
      groupsize=nsky
      bkpt=0
      sset=bspline_iterfit(tempwave,double(tempflux),nord=nord,everyn=everyn,bkpt=bkpt,maxrej=0,yfit=skyfit) 
    endelse

    ; Interpolate the fit to the wavelength of each fiber in the block
    tempskyfit=interpol(skyfit,tempwave,loglam2[*,dofibers])
    ; And subtract it off
    skysub2d[*,dofibers]=skysub2c[*,dofibers]-tempskyfit*skylinemask2[*,dofibers]
    ; Technically I should probably change the ivar channel too...
  endif
endfor

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; Re-assemble the different fiber sizes back into the same structures
skysub=fltarr(nx,ny)
skysubivar=fltarr(nx,ny)
pixelmask=fltarr(nx,ny)
fibermask=fltarr(ny)
skylinemask=fltarr(nx,ny)

skysub[*,ftype2]=skysub2d[*,*]
; Modification June 23 2013
; For convenience of test-run analysis, zero out 3'' and 5'' fibers
;skysub[*,ftype3]=skysub3[*,*]
;skysub[*,ftype5]=skysub5[*,*]
skysub[*,ftype3]=0.
skysub[*,ftype5]=0.

skylinemask[*,ftype2]=skylinemask2[*,*]

; Define full-size version of the different sky passes
skysub_pass1=fltarr(nx,ny)
skysub_pass1[*,ftype2]=skysub2
skysub_pass2=fltarr(nx,ny)
skysub_pass2[*,ftype2]=skysub2b
skysub_pass3=fltarr(nx,ny)
skysub_pass3[*,ftype2]=skysub2c
skysub_pass4=fltarr(nx,ny)
skysub_pass4[*,ftype2]=skysub2d

skysubivar[*,ftype2]=skysubivar2d[*,*]
skysubivar[*,ftype3]=skysubivar3[*,*]
skysubivar[*,ftype5]=skysubivar5[*,*]

pixelmask[*,ftype2]=pixelmask2[*,*]
pixelmask[*,ftype3]=pixelmask3[*,*]
pixelmask[*,ftype5]=pixelmask5[*,*]

fibermask[ftype2]=fibermask2[*]
fibermask[ftype3]=fibermask3[*]
fibermask[ftype5]=fibermask5[*]

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; Save the sky-subtracted flux values as is, and now modify flambda.
flambda = skysub
flambdaivar = skysubivar
skyimg = flux - flambda

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; DRL- telluric correction would be here, but not done yet



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;----------
; Interpolate over masked pixels, just for aesthetic purposes

flambda = djs_maskinterp(flambda, flambdaivar EQ 0, /const, iaxis=0 )
skysub_pass1=djs_maskinterp(skysub_pass1, flambdaivar EQ 0, /const, iaxis=0 )
skysub_pass2=djs_maskinterp(skysub_pass2, flambdaivar EQ 0, /const, iaxis=0 )
skysub_pass3=djs_maskinterp(skysub_pass3, flambdaivar EQ 0, /const, iaxis=0 )
skysub_pass4=djs_maskinterp(skysub_pass4, flambdaivar EQ 0, /const, iaxis=0 )

;----------
; Combine FIBERMASK and PIXELMASK to FINALMASK
finalmask = pixelmask
ntrace=ny ;???
for itrace=0, ntrace-1 do $
  finalmask[*,itrace] = finalmask[*,itrace] OR fibermask[itrace]

;----------
; Disable some mask bits in regions where 'NODATA' is set
q_nodata = (finalmask AND sdss_flagval('SPPIXMASK','NODATA')) NE 0
discards = ['NEARBADPIXEL','LOWFLAT','SCATTEREDLIGHT','NOSKY']
for j=0, n_elements(discards)-1 do $
  finalmask = finalmask - q_nodata $
  * (finalmask AND sdss_flagval('SPPIXMASK',discards[j]))

;----------
; Get an estimate of the relative chi^2 at each pixel.
; Do this with a simple linear interpolation.
if (keyword_set(relchi2set)) then begin
  xx = 0
  traceset2xy, vacset, xx, loglam
;   rchi2img = interpol(relchi2struct.chi2, relchi2struct.wave, loglam)
  rchi2img = bspline_valu(loglam, relchi2set)
   ; Compute the mean relative chi2 of sky-subtraction, after masking
   ; bad regions of the CCD
  fval = sdss_flagval('SPPIXMASK','NOPLUG') $
    + sdss_flagval('SPPIXMASK','BADTRACE') $
    + sdss_flagval('SPPIXMASK','BADFLAT') $
    + sdss_flagval('SPPIXMASK','BADARC') $
    + sdss_flagval('SPPIXMASK','LOWFLAT') $
    + sdss_flagval('SPPIXMASK','NOSKY') $
    + sdss_flagval('SPPIXMASK','NODATA') $
    + sdss_flagval('SPPIXMASK','BADFLUXFACTOR')
  indx = where((finalmask AND fval) EQ 0 AND flambdaivar NE 0, ct)
  if (ct EQ 0) then skychi2 = 0. $
  else skychi2 = mean(rchi2img[indx])
endif else begin
  rchi2img = 0 * flambda + 1.
  skychi2 = 0.
endelse
sxaddpar, objhdr, 'SKYCHI2', skychi2, ' Mean chi^2 of sky-subtraction'

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


sxaddpar, objhdr, 'MANGAVER', getenv('MANGAVER')

; Determine output details from provided spsframefile name
outname=spsframefile
len=strlen(outname)
suffix=strmid(outname,len-3,3)
; If filename ended in .gz suffix, remove it.
; Will be added back when we gzip later
if (suffix eq '.gz') then outname=strmid(outname,0,len-3)

mwrfits, flambda, outname, objhdr, /create;sky subtracted flux
mwrfits, flambdaivar, outname   ; sky subtracted inverse variance
mwrfits, finalmask, outname ; final pixel mask
mwrfits, vacset, outname    ;trace-set for wavelength sol; wset
mwrfits, dispset, outname   ;trace-set for dispersion sol
mwrfits, slitmap, outname  ;slitmap used
mwrfits, skyimg, outname    ;sky flux
; Not sure what these last two are, just spit out same as was read in
mwrfits, ximg, outname      ;x pos on CCD
mwrfits, superflat, outname  ;superflat vector from quartz lamps

; Calculate ratio of noise to poisson
nphot_skysub=flambda*superflat*igain
nphot=flux*superflat*igain

; If a flatfield was specified, xply by flat values of each fiber
if (keyword_set(flatfile)) then begin
  nphot=nphot*flat
  nphot_skysub=nphot_skysub*flat
endif

pratio=nphot_skysub/sqrt(nphot)
; Make sure no NaN or Inf values, zero them out
pratio[where(finite(pratio) eq 0)]=0.
mwrfits, pratio, outname


; Put everything on the same wavelength solution just to be helpful
; for Bershady analysis (not to be used for science)
lambdaber=10.^loglam
flamber=fltarr((size(wave))[1],ny)
flamivarber=fltarr((size(wave))[1],ny)
maskber=fltarr((size(wave))[1],ny)
wavearray=fltarr((size(wave))[1],ny)

for bershady=0,ny-1 do begin
  flamber[*,bershady]=interpol(flambda[*,bershady],lambdaber[*,bershady],wave)
  flamivarber[*,bershady]=interpol(flambdaivar[*,bershady],lambdaber[*,bershady],wave)
  maskber[*,bershady]=interpol(finalmask[*,bershady],lambdaber[*,bershady],wave)
  wavearray[*,bershady]=wave[*]
endfor
mwrfits, flamber, outname
mwrfits, flamivarber, outname
mwrfits, maskber, outname
mwrfits, wavearray, outname

; gzip output file
spawn, ['gzip', '-f', outname], /noshell

;;;;;;;;;;;;;;;;
; Make an extra output file for sky subtraction testing June 23 2013
outname2=outname+'_skypasses.fits'
;;;;;;;;;;;;;;;;

; Put everything on a common wavelength grid
ss1=fltarr((size(wave))[1],ny)
ss2=fltarr((size(wave))[1],ny)
ss3=fltarr((size(wave))[1],ny)
ss4=fltarr((size(wave))[1],ny)
smask=fltarr((size(wave))[1],ny)
for i=0,ny-1 do begin
  ss1[*,i]=interpol(skysub_pass1[*,i],lambdaber[*,i],wave)
  ss2[*,i]=interpol(skysub_pass2[*,i],lambdaber[*,i],wave)
  ss3[*,i]=interpol(skysub_pass3[*,i],lambdaber[*,i],wave)
  ss4[*,i]=interpol(skysub_pass4[*,i],lambdaber[*,i],wave)
  smask[*,i]=interpol(skylinemask[*,i],lambdaber[*,i],wave)
endfor

mwrfits,ss1,outname2,objhdr, /create
mwrfits,ss2,outname2
mwrfits,ss3,outname2
mwrfits,ss4,outname2
mwrfits,smask,outname2
mwrfits,wavearray,outname2

; Close out the plot file
if (keyword_set(plotfile)) then device,/close

return,0
end
