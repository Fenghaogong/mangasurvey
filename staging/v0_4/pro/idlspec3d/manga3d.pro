;+
; function manga3d.pro
;
; This is the MAIN pipeline reduction program to do the 3d cube work.
; See README manga3d.readme for more information.
;
; Author:
;   David R. Law (Dunlap Institute; drlaw@di.utoronto.ca)
;   First written: 11/15/2012
;   Last modified: 04/03/2013
;
; REVISION HISTORY:
;   v1: 23-Nov-2012  D. Law
;       First written.
;   v1.1: 15-Jan-2013 D. Law
;       Revision from output log to use splog
;   v1.2: 05-Feb-2013 D. Law
;       Revised output format to condense products into single file
;       with multiple extensions.  Changed directory structure of
;       intermediate products to include versioning.
;   v1.3: 07-Mar-2013 D. Law
;       Added some S/N calculations for individual fibers.
;       Fixed inverse variance channel combination.
;   v1.4: 03-Apr-2013 D. Law
;       Modifications to combination routine to properly handle
;       inverse variance in making the data cubes.
;       Additional tweak in tags/v0_4 to fix minor bug in ivar usage.
;
; Program to do cube-stage reduction of MaNGA data
;
; The present version of the code loops over all exposures, for each
; one it reads in the bundlemap and slitmap.  If these are the same
; for each exposure some small increase in speed could be achieved
; by only reading these once.  However, it is possible that different
; bundles or pluggings might need to be combined together, so in
; interests of this flexibility read each time.
;-

pro manga3d, input=input, outname=outname, inspect=inspect, doall=doall

  ; 32-bit integers
  Compile_Opt DEFINT32

  ; Shared use variables
  common MANGA_SHARE, pixscale, platescale, tstart
  tstart=systime()
  pixscale=0.5; 0.5 arcsec output pixels
  platescale=60.; 60 microns per arcsec

  ; Default input settings to read m3dplan.par from the
  ; current working directory
  if (keyword_set(input)) then input=input $
  else input='m3dplan.par'

  ; Read which files to process into 'redx'
  status=mlreadrdx(input,redx,rhdr)
  ; If read failed, exit
  if status ne 0 then mlquitmanga3d,status

  ; Read flavor of reduction (sos or full)
  ; Default in case of bad entry is full
  flavor=yanny_par(rhdr,'flavor')
  if flavor ne 'sos' then flavor='full'

  ; Default output name is 'DataCube'
  if (keyword_set(outname)) then outname=outname $
  else outname='DataCube'

  ; Start logging
  splog,filename=strcompress(strmid(input,0,strpos(input,'.par'))+'.log')
  splog,'Reduction started at ',tstart
  splog,'Using reduction pipeline version ',getenv('MANGADRP_VER')
  splog,strcompress('Processing using flavor '+flavor)

  ; Pipeline status starts off as 0
  status=0

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Read in basic info on which files to process,
; how many exposures.  Determine total number of fibers when
; all exposures are combined, and set up some data structures.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ; Read calibration file to set default wavelength grid 'wave'
  ; NOTE that this isn't the same as what is set by reductions...
  ; Currently interpolating reduced spectra onto this grid
  status=mlsetwcalib(wave)
  ; If read failed, exit
  if status ne 0 then mlquitmanga3d,status

  ; Define basic parameters
  nexp=(size(redx))[1]; Number of exposures to combine
  nftotal=mlcounttotalfibers(redx.ifuname,nexp); Total number of fibers between all exposures
  nwave=(size(wave))[1]; Total number of spectral elements in default wavelength grid

  ; Set up vectors and variables
  objname=''; Target names (e.g., 'NGC4333')
  objra=0.D; Target RA
  objdec=0.D; Target DEC
  xrel=fltarr(nwave,nftotal)
  yrel=fltarr(nwave,nftotal)
  dataarray=fltarr(nwave,nftotal)
  ivararray=fltarr(nwave,nftotal)
  fiber_status=intarr(nftotal)
  bigmask=intarr(nwave,nftotal)

  ; Other things
  totalexptime=0.

  filepath=''

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Loop over input exposures: Read in all of the spectra from files,
; figure out x,y offsets and fiber status
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

runningindex=0; Index of what fiber has been filled

for i=0,nexp-1 do begin
  splog,''
  splog,strcompress('Processing exposure number '+string(i+1))

  ; Figure out which IFU is used for this exposure, and
  ; how many fibers it has
  nfiber=mlgetbundlesize(redx[i].ifuname)

  ; Read the appropriate bundlemap 'bmap'
  status=mlreadbm(redx[i].ifuname,redx[i].mjd,bmap)
  ; If read failed, exit
  if status ne 0 then mlquitmanga3d,status

  ; Inspect the bundle map?
  if (keyword_set(inspect)) then begin
    bundleok='y'
    read,'Bundlemap ok? (y/n): ',bundleok
    if bundleok ne 'y' then begin
      splog,strcompress('User flagged bundle '+ifuname+' for exposure '+string(i)+' as bad.')
      mlquitmanga3d,-17L
    endif
  endif

  ; Read the appropriate slitmap 'slitmap'
  status=mlreadslm(redx[i].plate,redx[i].mjd,slitmap)
  ; If read failed, exit
  if status ne 0 then mlquitmanga3d,status

  ; Filepath to data files
  if flavor eq 'full' then $
    filepath=strcompress(getenv('MANGA_SPECTRO_REDUX')+'/'+getenv('RUN2D')+'/'+string(redx[i].plate)+'/'+string(redx[i].mjd)+'/',/remove_all)
  if flavor eq 'sos' then $
    filepath=strcompress(getenv('MANGA_QUICK')+'/'+string(redx[i].mjd)+'/',/remove_all)

  ; If there is no exp*.par summary file of exposure details in the data directory
  ; then create it from FITS headers.  Also force its creation if the
  ; /doall flag was set
  expfile=strcompress('exp-'+string(redx[i].plate)+'-'+string(redx[i].mjd)+'.par',/remove_all)
  junk=findfile(filepath+expfile,count=ct)
  if ((ct eq 0)or(keyword_set(doall))) then begin
    status=mlmakeexp(filepath,expfile,flavor)
    ; If make failed, exit
    if status ne 0 then mlquitmanga3d,status
  endif

  ; Read key parameters from exp*.par summary file.  We could just
  ; calculate that stuff here instead of reading it from a .par file,
  ; but logging, and more importantly modifying things to test, are
  ; easier this way.
  status=mlreadexp(filepath+expfile,obsparam)
  ; If read failed, exit
  if status ne 0 then mlquitmanga3d,status
  ; Select only the line appropriate to chosen exposure number
  obsparam=obsparam[where(obsparam.expnum eq redx[i].expnum)]
  totalexptime+=obsparam.exptime

  ; As of 01/22/13, dither positions not in FITS headers so not in obsparam
  ; structure.  Kludge them from the redx structure into obsparam here
  obsparam.dposn=redx[i].dposn

  ; When multiple carts, need to have some kind of check here on whether
  ; the specified IFU name was actually installed on the plate at the time...

  ; Target information
  subslit=slitmap[where(slitmap.ifuname eq redx[i].ifuname)]
  if (i eq 0) then begin ; If first exposure, define object coordinates
    objra=subslit[0].ra
    objdec=subslit[0].dec
  endif
  ; If not first exposure, check objects name and coords match the first!
  if (i ne 0) then begin
    if (subslit[0].ra ne objra) then begin
      status=-51L
      splog,strcompress('RA mismatch between exposure '+string(redx[0].expnum) +' and exposure '+ string(redx[i].expnum) +': '+objra+' vs '+subslit[0].ra)
      mlquitmanga3d,status
    endif
    if (subslit[0].dec ne objdec) then begin
      status=-52L
      splog,strcompress('DEC mismatch between exposure '+string(redx[0].expnum) +' and exposure '+ string(redx[i].expnum) +': '+objdec+' vs '+subslit[0].dec)
      mlquitmanga3d,status
    endif
  endif

  ; Print name and coordinates to log
  splog,''
  splog,strcompress('Target RA: '+string(objra))
  splog,strcompress('Target DEC: '+string(objdec))
  splog,''

  ; Set up temporary data vectors
  spec=fltarr(nwave,nfiber)
  specivar=fltarr(nwave,nfiber)
  xinbundle=fltarr(nfiber)
  yinbundle=fltarr(nfiber)
  fstat=intarr(nfiber)

  ; Call mlimportspec. This reads the actual spectra from the CCDs, performs
  ; sky subtraction, flux calibration, combines blue+red channels, and
  ; interpolates everything to a common wavelength vector.
  ; Also populates xinbundle, yinbundle, and fthrough with positioning
  ; and throughput information for the appropriate bundle.
  status=mlimportspec(wave,redx[i],obsparam,slitmap,bmap,spec,specivar,xinbundle,yinbundle,fstat,filepath,flavor,doall=doall)

  ; Move spectra from the temporary vector to permanent vector
  dataarray[0,runningindex]=spec[*,*]
  ivararray[0,runningindex]=specivar[*,*]

  ; Collate dead fibers in bundle.  Matt used 1=good, I use 0=good...
  ; This leads to some awkwards code inconsistencies...
  for j=0,nfiber-1 do begin
    if fstat[j] eq 1 then fiber_status[runningindex+j]=0 $
    else fiber_status[runningindex+j]=1
  endfor

  ; Determine relative positioning of the fibers.
  ; This is COMPLICATED and multi-step, put it in a subroutine.
  tempx=fltarr(nwave,nfiber)
  tempy=fltarr(nwave,nfiber)

  status=mlsetifupositions(obsparam,objra,objdec,xinbundle,yinbundle,wave,tempx,tempy)

  ;  Given final grid of x,y positions, make a map showing S/N ratio
  ; in individual fibers for each exposure
  status=mlfibersn(tempx/pixscale,tempy/pixscale,fstat,wave,spec,specivar,FiberSNgMap,FiberSNrMap,FiberSNiMap,FiberSNTable)
  snfile=strcompress(filepath+getenv('MANGADRP_VER')+'/snmap-'+redx[i].expnum+'-'+redx[i].ifuname+'.fits')
  mwrfits,FiberSNgMap,snfile,/create
  mwrfits,FiberSNrMap,snfile
  mwrfits,FiberSNiMap,snfile
  FiberSNTable[*,0]+=objra
  FiberSNTable[*,1]+=objdec
  FiberSNTable[*,3]+=objra
  FiberSNTable[*,4]+=objdec
  FiberSNTable[*,6]+=objra
  FiberSNTable[*,7]+=objdec
  mwrfits,FiberSNTable,snfile

  ; Move positions from the temporary vector to permanent vector
  ; and rescale to units of output pixels
  xrel[0,runningindex]=tempx[*,*]/pixscale
  yrel[0,runningindex]=tempy[*,*]/pixscale

;  mldtest,xrel[2000,runningindex:runningindex+nfiber-1],yrel[2000,runningindex:runningindex+nfiber-1],dataarray[2000,runningindex:runningindex+nfiber-1]

  runningindex+=nfiber; Bump running index by the number of fibers just added
endfor

; TEST- measure position in an exp
;mldtest,xrel[2000,0:18],yrel[2000,0:18],dataarray[2000,0:18]

; end test

; Set up output paths, create directories if needed
outpath=getenv('MANGA_DIR')+'/3dredux/'
spawn,strcompress('mkdir -p '+outpath)
if (flavor eq 'sos') then outpath=outpath+'sos/' $
else outpath=outpath+'full/'
spawn,strcompress('mkdir -p '+outpath)
outpath=outpath+getenv('MANGADRP_VER')+'/'
spawn,strcompress('mkdir -p '+outpath)
outpath=outpath+redx[0].plate+'/'
spawn,strcompress('mkdir -p '+outpath)
outpath=outpath+redx[0].ifuname+'/'
spawn,strcompress('mkdir -p '+outpath)

; Now run the S/N ratio computation for ALL fibers together to make a
; pseudo-map of S/N^2
; DRL- shouldn't need to do this EVERY time??  Existence check?
status=mlfibersn(xrel,yrel,fiber_status,wave,dataarray,ivararray,FiberSNgMap,FiberSNrMap,FiberSNiMap,FiberSNTable,/flipstat)
snfile=strcompress(outpath+'/snmapSTACK-'+redx[0].ifuname+'.fits')
mwrfits,FiberSNgMap,snfile,/create
mwrfits,FiberSNrMap,snfile
mwrfits,FiberSNiMap,snfile

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Given final grid of x,y positions, make a map showing effective coverage.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; Make a high resln coverage map every 1000 Angstroms from 3600 to 9600
; for visual analysis purposes
splog,'Making coverage map'
status=mlhirescovmap(CovMap,xrel,yrel,fiber_status,wave)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Arrays are populated: now combine together into data cube
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; Array size slightly bigger than max offset plus a bit
;Xsize=fix(2*(max(abs(xrel)) > max(abs(yrel)))+5)
;Ysize=Xsize
; Tweak, set these to be constant for a given nfiber
nfiber=mlgetbundlesize(redx[0].ifuname)
if nfiber eq 127 then Xsize=72
if nfiber eq 19 then Xsize=33
if nfiber eq 61 then Xsize=53
Ysize=Xsize
xrel=xrel+Xsize/2.; Don't want negative locations, want 0-Xsize, and 0-Ysize
yrel=yrel+Ysize/2.
DataCube=fltarr(Xsize,Ysize,nwave)
IvarCube=fltarr(Xsize,Ysize,nwave)
interpimage=fltarr(Xsize,Ysize)

x=fltarr(nftotal)
y=fltarr(nftotal)
dat=fltarr(nftotal)
ivar_in=fltarr(nftotal)
mask=intarr(nftotal)
for k=0,nwave-1 do begin
  if (k eq 0) then T=systime()
  Tnow=systime()
  tperk=mltimer(T,Tnow)/(k > 1)
  if ((k mod 200)eq 0) then splog,'Looping cube over wavelength space: ',(k*100.)/nwave,'% complete, ', tperk*(nwave-k),' seconds remaining'

  x[*]=xrel[k,*]
  y[*]=yrel[k,*]
  dat[*]=dataarray[k,*]
  ivar_in[*]=ivararray[k,*]
  mask[*]=fiber_status[*]

  ; Tweak mask to include bad pixels at a given wavelength
  badloc=where((dat eq 0.) or (fiber_status eq 1) or (ivar_in le 0.))
  if (size(badloc))[0] ne 0 then mask[badloc]=1
  bigmask[k,*]=mask[*]

  goodloc=where(mask eq 0)

  ; Do the interpolation
  ; If there were good values, combine those values
  if (size(goodloc))[0] ne 0 then begin

    xgood=x[goodloc]
    ygood=y[goodloc]
    datgood=dat[goodloc]
    maskgood=mask[goodloc]
    ivar_ingood=ivar_in[goodloc]

    interpimage=mlcalifainterp_ivar(xgood,ygood,datgood,maskgood,[Xsize,Ysize],1.6/pixscale,0.7/pixscale,ivar_ingood,ivar_out)
    DataCube[*,*,k]=interpimage
    IvarCube[*,*,k]=ivar_out
  ; Otherwise default to zeroes
  endif else begin
    DataCube[*,*,i]=0.
    IvarCube[*,*,i]=0.
  endelse
endfor

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Make header info for the cube
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; Look at the bottom of spcoadd_v5.pro for info on more useful
; keywords and things
mkhdr,head,DataCube

fxaddpar,head,'CDELT1',pixscale
fxaddpar,head,'CDELT2',pixscale
fxaddpar,head,'CUNIT1','arcsec'
fxaddpar,head,'CUNIT2','arcsec'
fxaddpar,head,'CTYPE1','LINEAR'
fxaddpar,head,'CTYPE2','LINEAR'

fxaddpar,head,'CUNIT3','Angstroms'
fxaddpar,head,'CRVAL3',alog10(wave[nwave-1])
fxaddpar,head,'CD1_1',1
fxaddpar,head,'CD2_2',2
fxaddpar,head,'CD3_3',alog10(wave[1])-alog10(wave[0])
fxaddpar,head,'CRPIX3',nwave
fxaddpar,head,'DC-FLAG',1
fxaddpar,head,'CTYPE3','LOG10'
fxaddpar,head,'DISPAXIS',3

fxaddpar,head,'UNITS','1e-17 erg/s/cm^2/Angstrom'

fxaddpar,head,'SURVEY','MANGA'
fxaddpar,head,'MANGADRP_VER',getenv('MANGADRP_VER')


; Total integration time in decimal hours
fxaddpar,head,'EXPTIME',totalexptime/60.

outfile=strcompress(outpath+outname+'.fits',/remove_all)
   ; HDU #0 is cube
mwrfits, DataCube, outfile, head, /create
   ; HDU #1 is inverse variance
mwrfits, IvarCube, outfile
   ; HDU #2 is wavelength
mwrfits, wave, outfile
   ; HDU #3 is RSS
mwrfits, dataarray, outfile
   ; HDU #4 is coverage map
mwrfits, CovMap, outfile

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Clean up and quit
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  mlquitmanga3d,status

return
end
