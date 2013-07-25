;+
; NAME:
;   runManga
;
; PURPOSE:
;     Runs the Manga pipeline, assuming all dependencies are set up
;
; CALLING SEQUENCE:
;   runManga, [ indir=, outdir=, sos=, full=, ]
;
; INPUTS:
;   none
;
; OPTIONAL INPUTS:
;   cam   - set this only if the input data does not contain data for all cameras ;
;           array of string names e.g. ['b1','r1'] ; default is all 4 cameras = ['b1','b2','r1','r2']
;   inop  - set this to a string to indicate an optional sub-folder inside the input data directory where the data is actually located
;   
;   outop - set this to a string to indicate an optional sub-solder inside the output directory where the resulting data files are to be stored
;   
;   plate - set this to indicate what plate is being reduced, used to set the output directory for the FULL pipeline
;
;   survey - set this to indicate what survey you're running [default: MANGA] ;useful for BOSS diagnostic through MANGA pipeline         
;             
; KEYWORDS:
;   sos        - set this keyword to run the quick reduction pipeline (Son Of Spectra BOSS pipeline) ; if no keywords set, SOS is default
;   
;   full       - set this keyword to run the full data reduction pipeline
;   
;   visual     - set this keyword to display plots during the reduction process of various steps & outputs
;   
;   peaks      - set this keyword to use the GETPEAKS method for tracing initial fiber centers instead of original TRACE_CEN
;   
;   skipproc   - set this keyword to skip the intial processing steps with SDSSPROC ;
;                this keyword only used if input data has already been processed ; implemented for the simulated post-processed MANGA input
;                for testing purposes only; most likely will never be used with real data
;                
;   stopgap    - set this keyword to implement various stop-gaps in the code to get simulated input to run to completion, for testing purposes only
;
;   multiple   - set this keyword to indicate to run the reduction pipeline on all files in multiple sub-directories within your data directory
;                keywords INOP and OUTOP must not be set
;                 
;   SOS-ONLY:
;   
;   red/blue   - set any of these keywords to process only those files from the RED or BLUE cameras ; e.g. /red will reduce only cameras r1 and r2
;   
;   flat/arc/sci - set any of these keywords to process only those flavors ; e.g. /flat will reduce only the flats
;
;   FULL-ONLY:
;
;   calib      - set this keyword to indicate to only reduce calibration files
;

pro runManga, SOS=SOS, FULL=FULL, CAMS=cams, VISUAL=VISUAL, PEAKS=PEAKS, SKIPPROC=SKIPPROC, STOPGAP=STOPGAP, SURVEY=SURVEY, $
                RED=RED, BLUE=BLUE, FLAT=FLAT, SCI=SCI, ARC=ARC, INOP=INOP, OUTOP=OUTOP, PLATE=PLATE, MULTIPLE=MULTIPLE, CALIB=CALIB

on_error,0
compile_opt idl2

starttime = systime(/seconds)
thismem = memory()
maxmem = 0
   
inop = (n_elements(inop) eq 0) ? '' : inop
outop = (n_elements(outop) eq 0) ? '' : outop
plate = (n_elements(plate) eq 0) ? '' : string(plate)

;to run on multiple sub-directories inside your data directory
if keyword_set(MULTIPLE) then begin
  paths = file_search(getenv2('MANGA_SPECTRO_DATA')+'*',/test_directory)
  subdirs = strmid(paths,strlen(getenv2('MANGA_SPECTRO_DATA')))
  
  ;select only MJD subdirs
  subdirs = subdirs[where(valid_num(subdirs) eq 1)]
  paths = paths[where(valid_num(subdirs) eq 1)]
  plate = strarr(n_elements(subdirs))
  
  ;for full pipeline
  if keyword_set(FULL) then begin
    plug = file_search(paths,'plPlugMapM*')
    plate = strmid(plug,46,4)
  endif    
  
  for numsub=0,n_elements(subdirs)-1 do $
    runManga, SOS=SOS, FULL=FULL, CAMS=cams, VISUAL=VISUAL, PEAKS=PEAKS, SKIPPROC=SKIPPROC, STOPGAP=STOPGAP, SURVEY=SURVEY, $
                RED=RED, BLUE=BLUE, FLAT=FLAT, SCI=SCI, ARC=ARC, INOP=subdirs[numsub], OUTOP=subdirs[numsub], PLATE=plate[numsub], CALIB=CALIB
  print, 'FINISHED ALL REDUCTIONS'
  print, 'Completetion Time: '+strtrim((systime(/seconds)-starttime)/60.,2)+' in minutes'
  return              
endif


core = getenv2('MANGA_DIR')
indir = getenv2(getenv2('MANGA_SPECTRO_DATA')+inop,/add)
redux = getenv2(getenv2('MANGA_SPECTRO_REDUX')+getenv2('RUN2D')+strtrim(plate,2)+'/'+outop,/add) ; full pipeline output dir
outdir = getenv2(getenv2('MANGA_QUICK')+outop,/add)                         ; sos pipeline outut dir
;plugdir = getenv('SPECLOG_DIR')+'55538/'
plugdir = indir
plugfile = file_search(plugdir, 'plPlugMapM*')
if file_test(plugfile[0]) eq 0 then begin
  print, 'ERROR: NO PLUGMAP FOUND IN DATA DIRECTORY!'
  return
endif
plugfile = strmid(plugfile,strpos(plugfile,'/',/reverse_search)+1,strlen(plugfile))

;check if dirs exist, if not then create
if file_test(outdir,/directory) eq 0 then spawn, '\mkdir -p '+outdir 
if strlen(plate) gt 1 then if file_test(redux, /directory) eq 0 then spawn, '\mkdir -p '+redux

;grab data files
data = file_search(indir+'sdR*')
if data[0] eq '' then begin
  print, "ERROR: No DATA found in input directory"
  return
endif
names = strmid(data,strpos(data[0],'sdR'),strlen(data[0]))

;set some default stuff
cams = (n_elements(cams) eq 0) ? ['b1','b2','r1','r2'] : cams
survey = (n_elements(survey) eq 0) ? 'manga' : survey
if not keyword_set(SOS) and not keyword_set(FULL) then SOS=1

; set for MANGA December sky-test
cams = ['b1','r1'] 
sone = (strupcase(survey) eq 'MANGA') ? 1 : 0
peaks =  (strupcase(survey) eq 'MANGA') ? 1 : 0

setcolors, /system_variable, /silent

;run SOS pipeline
if keyword_set(SOS) then $
  aporeduce, names, indir=indir, outdir=outdir, plugdir=plugdir, plugfile=plugfile, /no_diskcheck, /no_lock, skipproc=skipproc, stopgap=stopgap, $
    peaks=peaks, survey=survey, camnames=cams, red=red, blue=blue, flat=flat, sci=sci, arc=arc, visual=visual,sone=sone, data='real'


;run FULL pipeline
if keyword_SET(FULL) then begin
  if keyword_set(VISUAL) then xdisplay=1

  ;outop not MJD
  if NOT valid_num(outop) then begin
    print, 'ERROR: OUTOP MUST BE OF THE FORM MJD NUMBER!'
    return
  endif

  ;no plate given
  if strlen(plate) eq 0 then begin
    print, 'ERROR: PLATE ID MUST BE SPECIFIED!'
    return
  endif

  planfile = file_search(redux,'*Plan2d*'+outop+'*')
  ;write plan file
  if file_test(planfile) eq 0 then begin
    spplan2d, topdir=getenv2('MANGA_SPECTRO_DATA'), mjd=long(outop), outdir=redux, /nosplog
    planfile = file_search(redux,'*Plan2d*'+outop+'*')
  endif

  ; Figure out slitmap, and pass it to other routines
  temp=mlreadslm(plate,inop,slitmap)

  ;plan = strmid(planfile,strpos(planfile,redux)+strlen(redux), strlen(planfile))
  mlreduce2d, planfile,survey=survey,skipproc=skipproc, stopgap=stopgap,outdir=redux, indir=indir,docams=cams, $
    visual=visual, xdisplay=xdisplay, sone=sone, calib=calib,slitmap=slitmap,/bershady,peaks=peaks
endif


print, 'FINISHED REDUCTIONS'
print, 'Completetion Time: '+strtrim((systime(/seconds)-starttime)/60.,2)+' in minutes'

end
