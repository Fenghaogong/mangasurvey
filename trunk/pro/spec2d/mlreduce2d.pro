;+
; NAME:
;   mlreduce2d
;
; PURPOSE:
;   Calling script for SPREDUCE that reduces a night of data according
;   to a plan file.
;
; CALLING SEQUENCE:
;   mlreduce2d, [ planfile, docams=, /do_telluric, /xdisplay, $
;    /writeflatmodel, /writearcmodel, /bbspec ]
;
; INPUTS:
;
; OPTIONAL INPUTS:
;   planfile   - Name(s) of output plan file; default to reducing all
;                plan files matching 'spPlan2d*.par'
;   docams     - Cameras to reduce; default to ['b1', 'b2', 'r1', 'r2']
;   do_telluric- Passed to EXTRACT_OBJECT
;   xdisplay   - Send plots to X display rather than to plot file
;   writeflatmodel - passed to SPCALIB via SPREDUCE to trigger writing
;                    out of flat model info to file.
;   writearcmodel  - passed to SPCALIB via SPREDUCE to trigger writing
;                    out of arc model info to file.
;   bbspec         - use bbspec extraction code
;
; OUTPUT:
;
; COMMENTS:
;   The following environment variables must be set:
;      BOSS_SPECTRO_DATA
;      SPECLOG_DIR
;      SPECFLAT_DIR
;   Look for raw FITS data files in BOSS_SPECTRO_DATA/MJD.
;   Look for plug map files in SPECLOG_DIR/MJD.
;   Look for spectroscopic flat files in SPECFLAT_DIR.
;
; EXAMPLES:
;
; BUGS:
;   This routine spawns the Unix command 'mkdir'.
;
; PROCEDURES CALLED:
;   cpbackup
;   idlspec2d_version()
;   idlutils_version()
;   splog
;   spreduce
;   yanny_free
;   yanny_par()
;   yanny_read
;
; INTERNAL SUPPORT ROUTINES:
;
; REVISION HISTORY:
;   02-Nov-1999  Written by David Schlegel, Princeton.
;      Apr-2010  Added "write[flat,arc]model" pass-through (A. Bolton, Utah)
;   15-Aug-2011  Added pass-through for spatial split of sky model (A. Bolton, Utah)
;-
;------------------------------------------------------------------------------

pro mlreduce2d, planfile, DOCAMS=docams, DO_TELLURIC=DO_TELLURIC, $
                   XDISPLAY=XDISPLAY, WRITEFLATMODEL=WRITEFLATMODEL, WRITEARCMODEL=WRITEARCMODEL, $
                   BBSPEC=BBSPEC, VISUAL=VISUAL, SURVEY=survey, OUTDIR=outdir, INDIR=indir, PLUGDIR=plugdir, $
                   SKIPPROC=SKIPPROC, RED=RED, BLUE=BLUE, CALIB=CALIB, SCI=SCI, ARC=ARC, FLAT=FLAT
                   
   on_error, 0
   compile_opt idl2
   
   ;set some colors to use for plots
   setcolors, /system_variables, /silent

   ;check for environment variables
   rawdata_dir = getenv2('MANGA_SPECTRO_DATA')
   specflat_dir = getenv2('SPECFLAT_DIR')
   if (~keyword_set(rawdata_dir)) then message, 'Must set environment variable MANGA_SPECTRO_DATA' 
   if (~keyword_set(specflat_dir)) then message, 'Must set environment variable SPECFLAT_DIR' 
      
   ;set default directories
   if n_elements(indir) eq 0 then indir='./' 
   if n_elements(outdir) eq 0 then outdir='./'
   if n_elements(plugdir) eq 0 then plugdir=getenv2('MANGACORE_DIR')+getenv2('MANGACORE_VER')+'plugmaps/'

   ;check if spPlan2d file exists
   if (n_elements(planfile) eq 0) then begin
      planfile = file_search(outdir, 'spPlan2d*.par')
      message, 'ERROR: No spPlan2d file exists!'
   endif
   
;--------------------LOOP OVER MANY PLAN FILES
   ;if multiple plan files exist, then call this script recursively for each such plan file.
   if (N_elements(planfile) GT 1) then begin
      for i=0, N_elements(planfile)-1 do $
       mlreduce2d, planfile[i], docams=docams, do_telluric=do_telluric, $
        xdisplay=xdisplay, writeflatmodel=writeflatmodel, writearcmodel=writearcmodel, bbspec=bbspec, $
        visual=visual, survey=survey, outdir=outdir, indir=indir, plugdir=plugdir, skipproc=skipproc, $
        red=red, blue=blue, calib=calib, sci=sci, arc=arc, flat=flat
      return
   endif

   ;set the cameras if not set
   if (~keyword_set(docams)) then docams = ['b1', 'r1', 'b2', 'r2']

   thismem = memory()
   maxmem = 0

;---------------------- CONSTRUCT INFO FROM PLANFILE

   ; Strip path from plan file name, and change to that directory
   thisplan = (strmid(planfile[0],0,2) eq 'sp') ? fileandpath(outdir+planfile[0], path=thispath) : fileandpath(planfile[0],path=thispath)
   cd, thispath, current=origdir
   if (~keyword_set(thispath)) then cd, origdir

   ; Find the SPEXP structure
   allseq = yanny_readone(thisplan, 'SPEXP', hdr=hdr, /anon)
   if (N_elements(allseq) EQ 0) then begin
      message, 'ABORT: No SPEXP structures in plan file ' + thisplan
      cd, origdir
      return
   endif
    
   ; Find keywords from the header
   run2d = strtrim(string(yanny_par(hdr,'RUN2D')),2)
   mjd = long(yanny_par(hdr, 'MJD'))
   mjdstr = string(mjd, format='(i05.5)')

;----------------------
   ;Finalize some directory paths
   inputdir = (n_elements(indir) eq 0) ? concat_dir(rawdata_dir, mjdstr) : indir
   plugdir = (n_elements(plugdir) eq 0) ? getenv2('MANGACORE_DIR')+getenv2('MANGACORE_VER')+'plugmaps/' : plugdir

;----------------------- START LOGGING
   stime0 = systime(1)
   platemjd = string(yanny_par(hdr,'plateid'), format='(i04.4)') + '-' + mjdstr
   logfile = 'spDiag2d-' + platemjd + '.log'
   plotfile = 'spDiag2d-' + platemjd + '.ps' 
 
   ;Open and start a log file
   if (keyword_set(logfile)) then begin
      cpbackup, logfile
      splog, filename=logfile
      splog, 'Log file ' + logfile + ' opened ' + systime()
   endif
   
   ;Open and start a plot file for output
   if (keyword_set(plotfile) AND ~keyword_set(xdisplay)) then begin
      cpbackup, plotfile
      set_plot, 'ps'
      device, filename=plotfile, /color
      splog, 'Plot file ' + plotfile
   endif else set_plot, 'x'
   
   ;Initial log
   splog, 'IDL version: ' + string(!version,format='(99(a," "))')
   spawn, 'uname -a', uname
   splog, 'UNAME: ' + uname[0]
   splog, 'DISPLAY=' + getenv('DISPLAY')

   splog, 'idlspec2d version ' + idlspec2d_version()
   splog, 'idlutils version ' + idlutils_version()
   spawn, 'specflat_version', flatvers, /noshell
   splog, 'specflat version ' + flatvers[0]
   spawn, 'speclog_version', slogvers, /noshell
   splog, 'speclog version ' + slogvers[0]

   splog, 'Plan file ' + thisplan
   splog, 'DOCAMS = ', docams

   ;set camera names
   camnames = docams
   ncam = N_elements(camnames) 
 
;----------------FILTER FILES 
    ;select out specific flavors
   if keyword_set(calib) then allseq = allseq[where(allseq.flavor ne 'science')]
   if keyword_set(arc) then allseq = allseq[where(allseq.flavor eq 'arc')]
   if keyword_set(flat) then allseq = allseq[where(allseq.flavor eq 'flat')]
   if keyword_set(sci) then allseq = allseq[where(allseq.flavor eq 'science')] 

   ;select out red or blue camera  [b1, r1, b2, r2]
   if keyword_set(red) then allseq.name[[0,2]] = 'UNKNOWN'
   if keyword_set(blue) then allseq.name[[1,3]] = 'UNKNOWN'
 
;-----------------START REDUCTIONS
   ; Find all the unique plate plugging names
   allnames = allseq[ sort(allseq.mapname) ].mapname
   allnames = allnames[ uniq(allnames) ]
   
   ;loop over all plugging names   
   for imap=0, N_elements(allnames)-1 do begin
      ; Get the plate ID number from any (e.g., the first) exposure with this sequence ID number
      thismap = allnames[imap]
      j = where(allseq.mapname EQ thismap)
      plateid = allseq[j[0]].plateid
      platestr = string(plateid, format='(i4.4)')

      ;get start time of all reductions on plate, with plugging name imap
      stime1 = systime(1)
      splog, 'Begin plate ' + platestr + ' at ' + systime()

      ; Find the corresponding plug map file
      plugfile = 'plPlugMapM-' + thismap + '.par'
      splog, 'Plug map file = ', plugfile

      ;loop over all camera names [b1,r1,b2,r2]
      for ido=0, n_elements(docams)-1 do begin

         icam = (where(camnames EQ docams[ido], camct))[0]
         splog, prelog=camnames[icam]
         if (camct NE 1) then message, 'Non-unique camera ID: ' + docams[ido]

         ; Set the flag for splitting the sky model between spatial CCD halves:
         ; (Re: ticket #1388: strange r2 amplifier-boundary break)
         if ((camnames[icam] eq 'r2') and (mjd ge 55300)) then splitsky = 1B else splitsky = 0B

         ; Find science frames ; objname = String array with all science exposures at this sequence + camera
         objind = where(allseq.mapname EQ thismap AND (allseq.flavor EQ 'science' OR allseq.flavor EQ 'smear') AND allseq.name[icam] NE 'UNKNOWN' )
         if (objind[0] eq -1) then objname = '' else objname = allseq[objind].name[icam]


         if ~ifall(allseq.name[icam] eq 'UNKNOWN') then begin

             ; Select **all** flat exposures at this sequence + camera
             flatind = where(allseq.mapname EQ thismap AND allseq.flavor EQ 'flat' AND allseq.name[icam] NE 'UNKNOWN', nflat )
             if (nflat GT 0) then begin
                flatname = allseq[flatind].name[icam]
             endif else begin
                flatname = ''
                if (keyword_set(SCI) or keyword_set(ARC)) then splog, 'SCIENCE Keyword set, thus skipping flats.  Assuming flats exist for MAPNAME= ' + thismap + ', PLATEID= ' + platestr + ', CAMERA= ' + camnames[icam] else $
                   splog, 'ABORT: No flat for MAPNAME= ' + thismap + ', PLATEID= ' + platestr + ', CAMERA= ' + camnames[icam]
             endelse
    
             ; Select **all** arc exposures at this sequence + camera
             arcind = where(allseq.mapname EQ thismap AND allseq.flavor EQ 'arc' AND allseq.name[icam] NE 'UNKNOWN', narc )
             if (narc GT 0) then begin
                arcname = allseq[arcind].name[icam]
             endif else begin
                arcname = ''
                if (keyword_set(SCI) or keyword_set(FLAT)) then splog, 'SCIENCE Keyword set, thus skipping arcs.  Assuming arcs exist for MAPNAME= ' + thismap + ', PLATEID= ' + platestr + ', CAMERA= ' + camnames[icam] else $ 
                   splog, 'ABORT: No arc for MAPNAME= ' + thismap + ', PLATEID= ' + platestr + ', CAMERA= ' + camnames[icam]
             endelse
    
             ;Get start time before reductions of camera ido
             stime2 = systime(1)
    
             print, 'flats: ', flatname
             print, 'arcs: ', arcname
             print, 'sci: ', objname
                
                
             ;FIX THIS HERE   
                
             ; Reduce this set of frames (all objects w/same plate + camera)
             ;if (keyword_set(arcname) AND keyword_set(flatname)) then begin
             ; if ~ifall(allseq.name[icam] eq 'UNKNOWN') then begin
                 plottitle = 'PLATE='+platestr + ' MJD='+strtrim(string(mjd),2)+' '
                mlreduce, flatname, arcname, objname, run2d=run2d, plugfile=plugfile, lampfile=lampfile, indir=inputdir, plugdir=plugdir, outdir=outdir, $
                 plottitle=plottitle, do_telluric=do_telluric, writeflatmodel=writeflatmodel, writearcmodel=writearcmodel,$
                 bbspec=bbspec, splitsky=splitsky, visual=visual, red=red, blue=blue, calib=calib, sci=sci, arc=arc, flat=flat
             ;endif
    
             splog, 'Time to reduce camera ', camnames[icam], ' = ', systime(1)-stime2, ' seconds', format='(a,a,a,f6.0,a)'
    
             heap_gc   ; garbage collection
         endif ;end of "if (objind[0] ne -1)" original

         splog, prelog=''
      endfor ; End loop for camera number

      splog, 'Time to reduce all cameras = ', systime(1)-stime1, ' seconds', format='(a,f6.0,a)'

   endfor ; End loop for plugging name

;---------------- END REDUCTIONS
   ; Track memory usage
   thismem = memory()
   maxmem = maxmem > thismem[3]
   splog, 'Max memory usage = ', string(maxmem/1e6,format='(f7.1)'), ' MB'

   splog, 'Total time for SPREDUCE2D = ', systime(1)-stime0, ' seconds', format='(a,f6.0,a)'
   splog, 'Successful completion of SPREDUCE2D at ' + systime()

   ; Close log and plot files
   if (keyword_set(plotfile) AND ~keyword_set(xdisplay)) then begin
      device, /close
      set_plot, 'x'
   endif

   if (keyword_set(logfile)) then splog, /close

   ; Change back to original directory
   cd, origdir
   return

end

  