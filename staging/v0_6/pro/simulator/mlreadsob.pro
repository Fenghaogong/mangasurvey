; function mlreadsob
;
; Reads a simulated observing (.sob) setup file
; This file does all of the figuring out which observations and dither
; positions to use, so that these details are transparent in the
; simulation code itself.
;
; The only required input parameter is the sample observation .sob file
;
; This code is designed to fill the following parameters with
; pertinent single values:
;
; This code also defines the dimensions of the following parameters
; and fills them with pertinent values:
;
;
; Author:
;   David R. Law (Dunlap Institute; drlaw@di.utoronto.ca)
;   First written: 10/03/2012
;   Last modified: 01/17/2013
;
; REVISION HISTORY:
;   v1: 09-Oct-2012  D. Law
;       Adapted previous bundle mapping code to pipeline development
;       version.  Removed archaic bundle options.  Added 169-fiber bundles
;       and 2D addressing.
;  v1.1: 17-Jan-2013  D. Law
;       Tweaked to read in crosstalk value, fix read bug

forward_function mlradecsex
forward_function mlhasex

function mlreadsob,sobfile,pixscaleW,pixscaleO,boxxplier,ditherfrac,rseed,xtalk,guidelam,tc,rh,atmp,BType,BMfile,decl,knownxyerr,unknownxyerr,deadf,nexp,PA,seeing,exptime,ha,dposn

  openr,lunsob,sobfile,/get_lun

  junk='';String for reading comment blocks

  ; Read first comment block
  for i=0,2 do readf,lunsob,junk

  ; Read simulator parameters block
  pixscaleW=0.0
  pixscaleO=0.0
  boxxplier=0.0
  ndither=0
  ditherfrac=0.0
  rseed=0
  xtalk=0.
  guidelam=0.
  tc=0.
  rh=0.
  atmp=0.
  for i=0,1 do readf,lunsob,junk
  readf,lunsob,pixscaleW,junk
  readf,lunsob,pixscaleO,junk
  readf,lunsob,boxxplier,junk
  readf,lunsob,ndither,junk
  readf,lunsob,ditherfrac,junk
  readf,lunsob,rseed,junk
  readf,lunsob,xtalk,junk
  readf,lunsob,guidelam,junk
  readf,lunsob,tc,junk
  readf,lunsob,rh,junk
  readf,lunsob,atmp,junk
  readf,lunsob,junk

  ; Read target, bundle, and basic setup block
  BType=0
  BMfile=''
  dec=''
  knownxyerr=0.
  unknownxyerr=0.
  deadf1=0
  deadf2=0
  deadf3=0
  for i=0,2 do readf,lunsob,junk
  readf,lunsob,BType,junk
  readf,lunsob,BMfile
  BMfile=strmid(BMfile,0,strpos(BMfile,' ')); Strip out filename
  readf,lunsob,dec
  dec=strmid(dec,0,strpos(dec,' ')); Strip out declination
  mlradecsex,'00:00:00',dec,junk1,decl,/inpsex; Convert dec to double precision and decimal degrees
  readf,lunsob,knownxyerr,junk
  readf,lunsob,unknownxyerr,junk
  readf,lunsob,deadf1
  readf,lunsob,deadf2
  readf,lunsob,deadf3
  deadf=intarr(3)
  deadf[0]=deadf1
  deadf[1]=deadf2
  deadf[2]=deadf3

  ; Read observation details
  nnight=0
  nbad=0
  for i=0,5 do readf,lunsob,junk
  readf,lunsob,nnight,junk
  readf,lunsob,nbad,junk
  readf,lunsob,junk

  ; Calculate total number of good exposures to deal with
  nexp=ndither*nnight-nbad
  ; Define counting index to assign values with
  j=0

  ; Set matrix sizes for position angle, seeing, exposure time, hour angle
  ; These will have values for each exposure
  PA=dblarr(nexp)
  dposn=intarr(nexp)
  seeing=dblarr(nexp)
  exptime=dblarr(nexp)
  ha=dblarr(nexp)

  ; set temporary variable for read buffer
  patemp=0.
  qtemp=0
  dtemp=0
  seetemp=0.
  timetemp=0.
  hatemp=''
  hatemp2=0.

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; Read details for Night 1
  for i=0,1 do readf,lunsob,junk
  readf,lunsob,patemp,junk
  readf,lunsob,junk
  ; Exposure 1
  readf,lunsob,qtemp,junk
  readf,lunsob,dtemp,junk
  readf,lunsob,seetemp,junk
  readf,lunsob,timetemp,junk
  readf,lunsob,hatemp
  hatemp=strmid(hatemp,0,strpos(hatemp,' ')); Strip out hour angle
  mlhasex,hatemp,hatemp2,/inpsex; Convert hour angle from string to decimal degrees
  if (qtemp eq 0) then begin
    PA[j]=patemp
    dposn[j]=dtemp
    seeing[j]=seetemp
    exptime[j]=timetemp
    ha[j]=hatemp2
    j=j+1
  endif
  readf,lunsob,junk
  if (j eq nexp) then begin; If all exposures accounted for, stop scanning
    close,lunsob
    free_lun,lunsob
    return,0
  endif
  ; Exposure 2
  readf,lunsob,qtemp,junk
  readf,lunsob,dtemp,junk
  readf,lunsob,seetemp,junk
  readf,lunsob,timetemp,junk
  readf,lunsob,hatemp
  hatemp=strmid(hatemp,0,strpos(hatemp,' ')); Strip out hour angle
  mlhasex,hatemp,hatemp2,/inpsex; Convert hour angle from string to decimal degrees
  if (qtemp eq 0) then begin
    PA[j]=patemp
    dposn[j]=dtemp
    seeing[j]=seetemp
    exptime[j]=timetemp
    ha[j]=hatemp2
    j=j+1
  endif
  readf,lunsob,junk
  if (j eq nexp) then begin; If all exposures accounted for, stop scanning
    close,lunsob
    free_lun,lunsob
    return,0
  endif
  ; Exposure 3
  readf,lunsob,qtemp,junk
  readf,lunsob,dtemp,junk
  readf,lunsob,seetemp,junk
  readf,lunsob,timetemp,junk
  readf,lunsob,hatemp
  hatemp=strmid(hatemp,0,strpos(hatemp,' ')); Strip out hour angle
  mlhasex,hatemp,hatemp2,/inpsex; Convert hour angle from string to decimal degrees
  if (qtemp eq 0) then begin
    PA[j]=patemp
    dposn[j]=dtemp
    seeing[j]=seetemp
    exptime[j]=timetemp
    ha[j]=hatemp2
    j=j+1
  endif
  readf,lunsob,junk
  if (j eq nexp) then begin; If all exposures accounted for, stop scanning
    close,lunsob
    free_lun,lunsob
    return,0
  endif

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; Read details for Night 2
  for i=0,1 do readf,lunsob,junk
  readf,lunsob,patemp,junk
  readf,lunsob,junk
  ; Exposure 1
  readf,lunsob,qtemp,junk
  readf,lunsob,dtemp,junk
  readf,lunsob,seetemp,junk
  readf,lunsob,timetemp,junk
  readf,lunsob,hatemp
  hatemp=strmid(hatemp,0,strpos(hatemp,' ')); Strip out hour angle
  mlhasex,hatemp,hatemp2,/inpsex; Convert hour angle from string to decimal degrees
  if (qtemp eq 0) then begin
    PA[j]=patemp
    dposn[j]=dtemp
    seeing[j]=seetemp
    exptime[j]=timetemp
    ha[j]=hatemp2
    j=j+1
  endif
  readf,lunsob,junk
  if (j eq nexp) then begin; If all exposures accounted for, stop scanning
    close,lunsob
    free_lun,lunsob
    return,0
  endif
  ; Exposure 2
  readf,lunsob,qtemp,junk
  readf,lunsob,dtemp,junk
  readf,lunsob,seetemp,junk
  readf,lunsob,timetemp,junk
  readf,lunsob,hatemp
  hatemp=strmid(hatemp,0,strpos(hatemp,' ')); Strip out hour angle
  mlhasex,hatemp,hatemp2,/inpsex; Convert hour angle from string to decimal degrees
  if (qtemp eq 0) then begin
    PA[j]=patemp
    dposn[j]=dtemp
    seeing[j]=seetemp
    exptime[j]=timetemp
    ha[j]=hatemp2
    j=j+1
  endif
  readf,lunsob,junk
  if (j eq nexp) then begin; If all exposures accounted for, stop scanning
    close,lunsob
    free_lun,lunsob
    return,0
  endif
  ; Exposure 3
  readf,lunsob,qtemp,junk
  readf,lunsob,dtemp,junk
  readf,lunsob,seetemp,junk
  readf,lunsob,timetemp,junk
  readf,lunsob,hatemp
  hatemp=strmid(hatemp,0,strpos(hatemp,' ')); Strip out hour angle
  mlhasex,hatemp,hatemp2,/inpsex; Convert hour angle from string to decimal degrees
  if (qtemp eq 0) then begin
    PA[j]=patemp
    dposn[j]=dtemp
    seeing[j]=seetemp
    exptime[j]=timetemp
    ha[j]=hatemp2
    j=j+1
  endif
  readf,lunsob,junk
  if (j eq nexp) then begin; If all exposures accounted for, stop scanning
    close,lunsob
    free_lun,lunsob
    return,0
  endif

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; Read details for Night 3
  for i=0,1 do readf,lunsob,junk
  readf,lunsob,patemp,junk
  readf,lunsob,junk
  ; Exposure 1
  readf,lunsob,qtemp,junk
  readf,lunsob,dtemp,junk
  readf,lunsob,seetemp,junk
  readf,lunsob,timetemp,junk
  readf,lunsob,hatemp
  hatemp=strmid(hatemp,0,strpos(hatemp,' ')); Strip out hour angle
  mlhasex,hatemp,hatemp2,/inpsex; Convert hour angle from string to decimal degrees
  if (qtemp eq 0) then begin
    PA[j]=patemp
    dposn[j]=dtemp
    seeing[j]=seetemp
    exptime[j]=timetemp
    ha[j]=hatemp2
    j=j+1
  endif
  readf,lunsob,junk
  if (j eq nexp) then begin; If all exposures accounted for, stop scanning
    close,lunsob
    free_lun,lunsob
    return,0
  endif
  ; Exposure 2
  readf,lunsob,qtemp,junk
  readf,lunsob,dtemp,junk
  readf,lunsob,seetemp,junk
  readf,lunsob,timetemp,junk
  readf,lunsob,hatemp
  hatemp=strmid(hatemp,0,strpos(hatemp,' ')); Strip out hour angle
  mlhasex,hatemp,hatemp2,/inpsex; Convert hour angle from string to decimal degrees
  if (qtemp eq 0) then begin
    PA[j]=patemp
    dposn[j]=dtemp
    seeing[j]=seetemp
    exptime[j]=timetemp
    ha[j]=hatemp2
    j=j+1
  endif
  readf,lunsob,junk
  if (j eq nexp) then begin; If all exposures accounted for, stop scanning
    close,lunsob
    free_lun,lunsob
    return,0
  endif
  ; Exposure 3
  readf,lunsob,qtemp,junk
  readf,lunsob,dtemp,junk
  readf,lunsob,seetemp,junk
  readf,lunsob,timetemp,junk
  readf,lunsob,hatemp
  hatemp=strmid(hatemp,0,strpos(hatemp,' ')); Strip out hour angle
  mlhasex,hatemp,hatemp2,/inpsex; Convert hour angle from string to decimal degrees
  if (qtemp eq 0) then begin
    PA[j]=patemp
    dposn[j]=dtemp
    seeing[j]=seetemp
    exptime[j]=timetemp
    ha[j]=hatemp2
    j=j+1
  endif
  readf,lunsob,junk
  if (j eq nexp) then begin; If all exposures accounted for, stop scanning
    close,lunsob
    free_lun,lunsob
    return,0
  endif

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; Read details for Night 4
  for i=0,1 do readf,lunsob,junk
  readf,lunsob,patemp,junk
  readf,lunsob,junk
  ; Exposure 1
  readf,lunsob,qtemp,junk
  readf,lunsob,dtemp,junk
  readf,lunsob,seetemp,junk
  readf,lunsob,timetemp,junk
  readf,lunsob,hatemp
  hatemp=strmid(hatemp,0,strpos(hatemp,' ')); Strip out hour angle
  mlhasex,hatemp,hatemp2,/inpsex; Convert hour angle from string to decimal degrees
  if (qtemp eq 0) then begin
    PA[j]=patemp
    dposn[j]=dtemp
    seeing[j]=seetemp
    exptime[j]=timetemp
    ha[j]=hatemp2
    j=j+1
  endif
  readf,lunsob,junk
  if (j eq nexp) then begin; If all exposures accounted for, stop scanning
    close,lunsob
    free_lun,lunsob
    return,0
  endif
  ; Exposure 2
  readf,lunsob,qtemp,junk
  readf,lunsob,dtemp,junk
  readf,lunsob,seetemp,junk
  readf,lunsob,timetemp,junk
  readf,lunsob,hatemp
  hatemp=strmid(hatemp,0,strpos(hatemp,' ')); Strip out hour angle
  mlhasex,hatemp,hatemp2,/inpsex; Convert hour angle from string to decimal degrees
  if (qtemp eq 0) then begin
    PA[j]=patemp
    dposn[j]=dtemp
    seeing[j]=seetemp
    exptime[j]=timetemp
    ha[j]=hatemp2
    j=j+1
  endif
  readf,lunsob,junk
  if (j eq nexp) then begin; If all exposures accounted for, stop scanning
    close,lunsob
    free_lun,lunsob
    return,0
  endif
  ; Exposure 3
  readf,lunsob,qtemp,junk
  readf,lunsob,dtemp,junk
  readf,lunsob,seetemp,junk
  readf,lunsob,timetemp,junk
  readf,lunsob,hatemp
  hatemp=strmid(hatemp,0,strpos(hatemp,' ')); Strip out hour angle
  mlhasex,hatemp,hatemp2,/inpsex; Convert hour angle from string to decimal degrees
  if (qtemp eq 0) then begin
    PA[j]=patemp
    dposn[j]=dtemp
    seeing[j]=seetemp
    exptime[j]=timetemp
    ha[j]=hatemp2
    j=j+1
  endif
  readf,lunsob,junk
  if (j eq nexp) then begin; If all exposures accounted for, stop scanning
    close,lunsob
    free_lun,lunsob
    return,0
  endif

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; Read details for Night 5
  for i=0,1 do readf,lunsob,junk
  readf,lunsob,patemp,junk
  readf,lunsob,junk
  ; Exposure 1
  readf,lunsob,qtemp,junk
  readf,lunsob,dtemp,junk
  readf,lunsob,seetemp,junk
  readf,lunsob,timetemp,junk
  readf,lunsob,hatemp
  hatemp=strmid(hatemp,0,strpos(hatemp,' ')); Strip out hour angle
  mlhasex,hatemp,hatemp2,/inpsex; Convert hour angle from string to decimal degrees
  if (qtemp eq 0) then begin
    PA[j]=patemp
    dposn[j]=dtemp
    seeing[j]=seetemp
    exptime[j]=timetemp
    ha[j]=hatemp2
    j=j+1
  endif
  readf,lunsob,junk
  if (j eq nexp) then begin; If all exposures accounted for, stop scanning
    close,lunsob
    free_lun,lunsob
    return,0
  endif
  ; Exposure 2
  readf,lunsob,qtemp,junk
  readf,lunsob,dtemp,junk
  readf,lunsob,seetemp,junk
  readf,lunsob,timetemp,junk
  readf,lunsob,hatemp
  hatemp=strmid(hatemp,0,strpos(hatemp,' ')); Strip out hour angle
  mlhasex,hatemp,hatemp2,/inpsex; Convert hour angle from string to decimal degrees
  if (qtemp eq 0) then begin
    PA[j]=patemp
    dposn[j]=dtemp
    seeing[j]=seetemp
    exptime[j]=timetemp
    ha[j]=hatemp2
    j=j+1
  endif
  readf,lunsob,junk
  if (j eq nexp) then begin; If all exposures accounted for, stop scanning
    close,lunsob
    free_lun,lunsob
    return,0
  endif
  ; Exposure 3
  readf,lunsob,qtemp,junk
  readf,lunsob,dtemp,junk
  readf,lunsob,seetemp,junk
  readf,lunsob,timetemp,junk
  readf,lunsob,hatemp
  hatemp=strmid(hatemp,0,strpos(hatemp,' ')); Strip out hour angle
  mlhasex,hatemp,hatemp2,/inpsex; Convert hour angle from string to decimal degrees
  if (qtemp eq 0) then begin
    PA[j]=patemp
    dposn[j]=dtemp
    seeing[j]=seetemp
    exptime[j]=timetemp
    ha[j]=hatemp2
    j=j+1
  endif
  readf,lunsob,junk
  if (j eq nexp) then begin; If all exposures accounted for, stop scanning
    close,lunsob
    free_lun,lunsob
    return,0
  endif

  close,lunsob
  free_lun,lunsob
return,1; If we made it this far, something failed
end
