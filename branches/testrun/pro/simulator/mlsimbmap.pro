; function mlsimbmap
;
; Old pre-test run code to simulate locations of fibers within
; a bundle of given size.

; Calculations individual fiber locations for MaNGA fiber bundles
; Heavily revised from early simulation code and downselected
; to focus on only MaNGA-style bundles as of October 2012.
; Fibers are assumed to be polymicro 120/132/150 micron
; fibers arranged in regular hexagonal bundles.

; Required input arguments:
; BType: Integer bundle type (-1, 7, 19, 37, 61, 91, 127, 169)
;        Note that BType=-1 will read fiber locations from a
;        .bm format input file
;
; Optional input arguments:
; BMfile: Specify .bm format input file for an as-built bundle
;   Code assumes that this file is in the
;   $MANGADB/metrology/
;   directory
;
; Returns 0 if everything worked, -1 if error
;
; Output arguments:
; nfiber: Number of fibers
; fiber_address: 2-component fiber address
; fiber_xcen and fiber_ycen: vectors of fiber centers
;  computed by this function, returned in arcseconds
; fiber_status: Fiber status (0=live, 1=dead)
; corediam: Fiber core diameter in arcsec
; dither_rad: radius offset for dithers in arcsec
;  (i.e., quantifies space between 'holes')
; bundle_rmax: Maximum fiber distance from center in arcsec
;  (useful for defining a simulation box size)
; /spiral: Flag to use spiral pattern if available
;
; Author:
;   David R. Law (Dunlap Institute; drlaw@di.utoronto.ca)
;   First written: 07/15/2011
;   Last modified: 01/17/2013
;
; REVISION HISTORY:
;   v1: 17-Oct-2012  D. Law
;       Adapted previous bundle mapping code to pipeline development
;       version.  Removed archaic bundle options.  Added 169-fiber bundles
;       and 2D addressing.
;   v1.2: 22-Dec-2012  D. Law
;       Adapted previous version to deal with Yanny-style bundle maps
;   v1.3: 17-Jan-2013  D. Law
;       Tweaked filepath to metrology data, implementation of dead
;       fibers, added spiral indexing option to as-built bundles.

; mlxyaddress converts from fiber address to x,y position
; in microns, given ideal fiber spacing values
pro mlxyaddress,fiber_address,fiber_xcen,fiber_ycen,xoff,yoff
  fiber_xcen[*]=(fiber_address[*,0]*2+fiber_address[*,1])*xoff
  fiber_ycen[*]=fiber_address[*,1]*yoff
return
end

function mlsimbmap,BType,nfiber,fiber_xcen,fiber_ycen,fiber_status,corediam,dither_rad,bundle_rmax,BMfile=BMfile,spiral=spiral

  ; Define fundamental geometric parameters common to all bundles
  pltscl=60.;microns/arcsecond
  corediam=120.; Fiber core diameter in microns
  fiberdiam=150.; Fiber buffer diameter in microns
  ;dither_rad=fiberdiam*0.577; Dither offset radius in microns (This was
  ; old offset distance, when position 1 was the zero-point)
  dither_rad=fiberdiam/3.; Dither offset radius from zero-position in microns
  ; (Note that this is the distance assuming all 3 positions are offset
  ; from a central point.  Division by 3.0 is EXACT: 3.0=(2*cos(30))**2)
  ; Define ideal fiber-fiber spacing
  xoff=fiberdiam*cos(60.*3.1415926536/180.)
  yoff=fiberdiam*sin(60.*3.1415926536/180.)

  ; If BType is unknown, exit
  if ((BType ne -1)and(BType ne 7)and(BType ne 19)and(BType ne 37)and(BType ne 61)and(BType ne 91)and(BType ne 127)and(BType ne 169)) then begin
     print,'ERROR: Unknown BType in mlbundlemap.  Exit!'
     return,1
  endif

  ; If BType is -1, read in fiber locations from the specified .bm file
  if (BType eq -1) then begin
    BMfile=strcompress(getenv('MANGAROOT')+'/mangadb/metrology/'+BMfile)

    bmdata=yanny_readone(BMfile,'BUNDLEMAP',hdr=hdr)
    fiber_xcen=-bmdata.xpmm*1000.
    fiber_ycen=bmdata.ypmm*1000.

    nfiber=(size(fiber_xcen))[1]
    ; Extract fiber status information
    fiber_status=intarr(nfiber)
    for i=0,nfiber-1 do begin
      if bmdata[i].gbu eq 1 then fiber_status[i]=0 $
      else fiber_status[i]=1
    endfor

   ; If /spiral keyword set, reorder to spiral indexing
   if keyword_set(spiral) then begin
     sporder=bmdata.isp+1; For some reason Matt uses 0-indexing on spiral order...
                         ; Convert to 1-indexing for consistency
     reorder=sort(sporder)
     fiber_xcen=fiber_xcen(reorder)
     fiber_ycen=fiber_ycen(reorder)
     fiber_status=fiber_status(reorder)
   endif

  endif

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;; 7-fiber minibundles (white)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  if (BType eq 7) then begin
    nfiber=7
    fiber_xcen=fltarr(nfiber)
    fiber_ycen=fltarr(nfiber)
    fiber_address=intarr(nfiber,2)
    fiber_status=intarr(nfiber)

    ; Fiber 0 is the lower-left corner, serpentine layout
    ; Row 1 (bottom row)
    for i=0,1 do begin
       fiber_address[i,0]=i
       fiber_address[i,1]=-1
    endfor
    ; Row 2 (middle row)
    for i=2,4 do begin
       fiber_address[i,0]=1-(i-2)
       fiber_address[i,1]=0
    endfor
    ; Row 2 (top row)
    for i=5,6 do begin
       fiber_address[i,0]=i-6
       fiber_address[i,1]=1
    endfor
  mlxyaddress,fiber_address,fiber_xcen,fiber_ycen,xoff,yoff
  endif

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;; 19-fiber bundles (red)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  if (BType eq 19) then begin
    nfiber=19
    fiber_xcen=fltarr(nfiber)
    fiber_ycen=fltarr(nfiber)
    fiber_address=intarr(nfiber,2)
    fiber_status=intarr(nfiber)

    ; Fiber 0 is the lower-left corner, serpentine layout
    ; Row 1 (bottom row)
    for i=0,2 do begin
       fiber_address[i,0]=i
       fiber_address[i,1]=-2
    endfor
    ; Row 2
    for i=3,6 do begin
       fiber_address[i,0]=2-(i-3)
       fiber_address[i,1]=-1
    endfor
    ; Row 3 (middle row)
    for i=7,11 do begin
       fiber_address[i,0]=(i-7)-2
       fiber_address[i,1]=0
    endfor
    ; Row 4
    for i=12,15 do begin
       fiber_address[i,0]=1-(i-12)
       fiber_address[i,1]=1
    endfor
    ; Row 5 (top row)
    for i=16,18 do begin
       fiber_address[i,0]=(i-16)-2
       fiber_address[i,1]=2
    endfor
  mlxyaddress,fiber_address,fiber_xcen,fiber_ycen,xoff,yoff
endif

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;; 37-fiber bundles (green)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  if (BType eq 37) then begin
    nfiber=37
    fiber_xcen=fltarr(nfiber)
    fiber_ycen=fltarr(nfiber)
    fiber_address=intarr(nfiber,2)
    fiber_status=intarr(nfiber)

    ; Fiber 0 is the lower-left corner, serpentine layout
    ; Row 1 (bottom row)
    for i=0,3 do begin
       fiber_address[i,0]=i
       fiber_address[i,1]=-3
    endfor
    ; Row 2
    for i=4,8 do begin
       fiber_address[i,0]=3-(i-4)
       fiber_address[i,1]=-2
    endfor
    ; Row 3
    for i=9,14 do begin
       fiber_address[i,0]=(i-9)-2
       fiber_address[i,1]=-1
    endfor
    ; Row 4 (middle row)
    for i=15,21 do begin
       fiber_address[i,0]=3-(i-15)
       fiber_address[i,1]=0
    endfor
    ; Row 5
    for i=22,27 do begin
       fiber_address[i,0]=(i-22)-3
       fiber_address[i,1]=1
    endfor
    ; Row 6
    for i=28,32 do begin
       fiber_address[i,0]=1-(i-28)
       fiber_address[i,1]=2
    endfor
    ; Row 7 (top row)
    for i=33,36 do begin
       fiber_address[i,0]=(i-33)-3
       fiber_address[i,1]=3
    endfor
  mlxyaddress,fiber_address,fiber_xcen,fiber_ycen,xoff,yoff
  endif

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;; 61-fiber bundles (orange)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  if (BType eq 61) then begin
    nfiber=61
    fiber_xcen=fltarr(nfiber)
    fiber_ycen=fltarr(nfiber)
    fiber_address=intarr(nfiber,2)
    fiber_status=intarr(nfiber)

    ; Fiber 0 is the lower-left corner, serpentine layout
    ; Row 1 (bottom row)
    for i=0,4 do begin
       fiber_address[i,0]=i
       fiber_address[i,1]=-4
    endfor
    ; Row 2
    for i=5,10 do begin
       fiber_address[i,0]=4-(i-5)
       fiber_address[i,1]=-3
    endfor
    ; Row 3
    for i=11,17 do begin
       fiber_address[i,0]=(i-11)-2
       fiber_address[i,1]=-2
    endfor
    ; Row 4
    for i=18,25 do begin
       fiber_address[i,0]=4-(i-18)
       fiber_address[i,1]=-1
    endfor
    ; Row 5 (middle row)
    for i=26,34 do begin
       fiber_address[i,0]=(i-26)-4
       fiber_address[i,1]=0
    endfor
    ; Row 6
    for i=35,42 do begin
       fiber_address[i,0]=3-(i-35)
       fiber_address[i,1]=1
    endfor
    ; Row 7
    for i=43,49 do begin
       fiber_address[i,0]=(i-43)-4
       fiber_address[i,1]=2
    endfor
    ; Row 8
    for i=50,55 do begin
       fiber_address[i,0]=1-(i-50)
       fiber_address[i,1]=3
    endfor
    ; Row 9 (top row)
    for i=56,60 do begin
       fiber_address[i,0]=(i-56)-4
       fiber_address[i,1]=4
    endfor
  mlxyaddress,fiber_address,fiber_xcen,fiber_ycen,xoff,yoff
endif

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;; 91-fiber bundles (blue)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  if (BType eq 91) then begin
    nfiber=91
    fiber_xcen=fltarr(nfiber)
    fiber_ycen=fltarr(nfiber)
    fiber_address=intarr(nfiber,2)
    fiber_status=intarr(nfiber)

    ; Fiber 0 is the lower-left corner, serpentine layout
    ; Row 1 (bottom row)
    for i=0,5 do begin
       fiber_address[i,0]=i
       fiber_address[i,1]=-5
    endfor
    ; Row 2
    for i=6,12 do begin
       fiber_address[i,0]=5-(i-6)
       fiber_address[i,1]=-4
    endfor
    ; Row 3
    for i=13,20 do begin
       fiber_address[i,0]=(i-13)-2
       fiber_address[i,1]=-3
    endfor
    ; Row 4
    for i=21,29 do begin
       fiber_address[i,0]=5-(i-21)
       fiber_address[i,1]=-2
    endfor
    ; Row 5
    for i=30,39 do begin
       fiber_address[i,0]=(i-30)-4
       fiber_address[i,1]=-1
    endfor
    ; Row 6 (middle row)
    for i=40,50 do begin
       fiber_address[i,0]=5-(i-40)
       fiber_address[i,1]=0
    endfor
    ; Row 7
    for i=51,60 do begin
       fiber_address[i,0]=(i-51)-5
       fiber_address[i,1]=1
    endfor
    ; Row 8
    for i=61,69 do begin
       fiber_address[i,0]=3-(i-61)
       fiber_address[i,1]=2
    endfor
    ; Row 9
    for i=70,77 do begin
       fiber_address[i,0]=(i-70)-5
       fiber_address[i,1]=3
    endfor
    ; Row 10
    for i=78,84 do begin
       fiber_address[i,0]=1-(i-78)
       fiber_address[i,1]=4
    endfor
    ; Row 11 (top row)
    for i=85,90 do begin
       fiber_address[i,0]=(i-85)-5
       fiber_address[i,1]=5
    endfor
  mlxyaddress,fiber_address,fiber_xcen,fiber_ycen,xoff,yoff
  endif

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;; 127-fiber bundles (purple)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  if (BType eq 127) then begin
    nfiber=127
    fiber_xcen=fltarr(nfiber)
    fiber_ycen=fltarr(nfiber)
    fiber_address=intarr(nfiber,2)
    fiber_status=intarr(nfiber)

    ; Fiber 0 is the lower-left corner, serpentine layout
    ; Row 1 (bottom row)
    for i=0,6 do begin
       fiber_address[i,0]=i
       fiber_address[i,1]=-6
    endfor
    ; Row 2
    for i=7,14 do begin
       fiber_address[i,0]=6-(i-7)
       fiber_address[i,1]=-5
    endfor
    ; Row 3
    for i=15,23 do begin
       fiber_address[i,0]=(i-15)-2
       fiber_address[i,1]=-4
    endfor
    ; Row 4
    for i=24,33 do begin
       fiber_address[i,0]=6-(i-24)
       fiber_address[i,1]=-3
    endfor
    ; Row 5
    for i=34,44 do begin
       fiber_address[i,0]=(i-34)-4
       fiber_address[i,1]=-2
    endfor
    ; Row 6
    for i=45,56 do begin
       fiber_address[i,0]=6-(i-45)
       fiber_address[i,1]=-1
    endfor
    ; Row 7 (middle row)
    for i=57,69 do begin
       fiber_address[i,0]=(i-57)-6
       fiber_address[i,1]=0
    endfor
    ; Row 8
    for i=70,81 do begin
       fiber_address[i,0]=5-(i-70)
       fiber_address[i,1]=1
    endfor
    ; Row 9
    for i=82,92 do begin
       fiber_address[i,0]=(i-82)-6
       fiber_address[i,1]=2
    endfor
    ; Row 10
    for i=93,102 do begin
       fiber_address[i,0]=3-(i-93)
       fiber_address[i,1]=3
    endfor
    ; Row 11
    for i=103,111 do begin
       fiber_address[i,0]=(i-103)-6
       fiber_address[i,1]=4
    endfor
    ; Row 12
    for i=112,119 do begin
       fiber_address[i,0]=1-(i-112)
       fiber_address[i,1]=5
    endfor
    ; Row 13 (top row)
    for i=120,126 do begin
       fiber_address[i,0]=(i-120)-6
       fiber_address[i,1]=6
    endfor
  mlxyaddress,fiber_address,fiber_xcen,fiber_ycen,xoff,yoff
  endif

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;; 169-fiber bundles (grey)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  if (BType eq 169) then begin
    nfiber=169
    fiber_xcen=fltarr(nfiber)
    fiber_ycen=fltarr(nfiber)
    fiber_address=intarr(nfiber,2)
    fiber_status=intarr(nfiber)

    ; Fiber 0 is the lower-left corner, serpentine layout
    ; Row 1 (bottom row)
    for i=0,7 do begin
       fiber_address[i,0]=i
       fiber_address[i,1]=-7
    endfor
    ; Row 2
    for i=8,16 do begin
       fiber_address[i,0]=7-(i-8)
       fiber_address[i,1]=-6
    endfor
    ; Row 3
    for i=17,26 do begin
       fiber_address[i,0]=(i-17)-2
       fiber_address[i,1]=-5
    endfor
    ; Row 4
    for i=27,37 do begin
       fiber_address[i,0]=7-(i-27)
       fiber_address[i,1]=-4
    endfor
    ; Row 5
    for i=38,49 do begin
       fiber_address[i,0]=(i-38)-4
       fiber_address[i,1]=-3
    endfor
    ; Row 6
    for i=50,62 do begin
       fiber_address[i,0]=7-(i-50)
       fiber_address[i,1]=-2
    endfor
    ; Row 7
    for i=63,76 do begin
       fiber_address[i,0]=(i-63)-6
       fiber_address[i,1]=-1
    endfor
    ; Row 8 (middle row)
    for i=77,91 do begin
       fiber_address[i,0]=7-(i-77)
       fiber_address[i,1]=0
    endfor
    ; Row 9
    for i=92,105 do begin
       fiber_address[i,0]=(i-92)-7
       fiber_address[i,1]=1
    endfor
    ; Row 10
    for i=106,118 do begin
       fiber_address[i,0]=5-(i-106)
       fiber_address[i,1]=2
    endfor
    ; Row 11
    for i=119,130 do begin
       fiber_address[i,0]=(i-119)-7
       fiber_address[i,1]=3
    endfor
    ; Row 12
    for i=131,141 do begin
       fiber_address[i,0]=3-(i-131)
       fiber_address[i,1]=4
    endfor
    ; Row 13
    for i=142,151 do begin
       fiber_address[i,0]=(i-142)-7
       fiber_address[i,1]=5
    endfor
    ; Row 14
    for i=152,160 do begin
       fiber_address[i,0]=1-(i-152)
       fiber_address[i,1]=6
    endfor
    ; Row 15 (top row)
    for i=161,168 do begin
       fiber_address[i,0]=(i-161)-7
       fiber_address[i,1]=7
    endfor
  mlxyaddress,fiber_address,fiber_xcen,fiber_ycen,xoff,yoff
endif

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;; All bundles: convert to arcsec
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; Convert from microns to arcsec using standard platescale
  fiber_xcen=fiber_xcen/pltscl
  fiber_ycen=fiber_ycen/pltscl
  dither_rad=dither_rad/pltscl
  corediam=corediam/pltscl
  bundle_rmax=max(sqrt(fiber_xcen*fiber_xcen+fiber_ycen*fiber_ycen))+fiberdiam/2./pltscl

return,0
end
