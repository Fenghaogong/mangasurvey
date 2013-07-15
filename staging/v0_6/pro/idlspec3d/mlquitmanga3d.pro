; program mlquitmanga3d
;
; Exit the manga3d program cleanly and return everything to base IDL level.
; Can be called early with a specified error code
;
; Author:
;   David R. Law (Dunlap Institute; drlaw@di.utoronto.ca)
;   First written: 11/19/2012
;   Last modified: 01/15/2013
;
; REVISION HISTORY:
;   v1: 19-Nov-2012  D. Law
;       First written
;   v1.1: 15-Jan-2013 D. Law
;       Tweaked logging
pro mlquitmanga3d,status
  ; Shared use variables
  common MANGA_SHARE, pixscale,platescale, tstart

  ; If status ne 0, there was an error, output a message about it
  if (status ne 0) then begin
    splog,strcompress('Error code '+string(status)+', quit!')
    splog,'Check the output log or error codes for more information'
  endif

  tstop=systime(); Completion time
  splog,''
  splog,'Reduction stopped at ',tstop
  splog,'Total time: ',mltimer(tstart,tstop),' seconds'

  ; Close down output log gracefully
  splog,/close

  ; Close all active programs and functions and return command
  ; to the main IDL prompt
  retall

end
