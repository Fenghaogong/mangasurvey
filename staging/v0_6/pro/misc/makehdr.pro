;; Procedure creates a FITS table header based on an existing idl structure and then helps the user to add column information and
;; comments to it.
;;
;; KBundy 6/1/11

pro makehdr, hdr, structure, info_file, INITIALIZE = INITIALIZE

if n_elements(hdr) LE 1 OR keyword_set(INITIALIZE) then begin  ;; The hdr has not yet been created

   fxbhmake, hdr, n_elements(structure), /date, /initialize
   mwrfits, structure[0], 'test.fits', hdr, /create  ;; We write out a dummy FITS table to get column info stored in hdr
   spawn, 'rm test.fits'                             ;; Delete it

   wttype = where(strmatch(hdr, '*TTYPE*'), nttype)

   if n_elements(info_file) EQ 0 then begin
      print, '  ***  Copy/Paste and then edit the following into your information file               ***'
      print, '  ***  Make sure the number of columns is 80 characters                                ***'
      print, '  ***  Add comments on a newline, beginning with COMMENT and 1st character on column 9 ***'
      print
      print, '# Column information for FITS table...'
      print, '# KBundy'
      print
      forprint, hdr[wttype], textout=2
      print
   endif else begin
      openw, 1, info_file, width=80 & !textunit=1
      printf, 1, '#  ***  Make sure the number of columns is 80 characters                     ***'
      printf, 1, '#  ***  Add comments on a newline, with COMMENT and 1st char on column 9     ***'

      printf, 1
      printf, 1, '# Column information for FITS table...'
      printf, 1, '#   '+systime()
      printf, 1
      printf, 1, '# COMMENT  Comment starts here'
      printf, 1
      forprint, hdr[wttype], textout=5, /nocomment
      printf, 1
      close, 1
   endelse
    
   return
endif else begin

   if n_elements(info_file) EQ 0 then begin
      print, 'Please specify a column information file'
      return
   endif

   ttype_vals = sxpar(hdr, 'TTYPE*', count = n_ttype, comment = ttype_comment)

   
   string = ''
   i = 0

   openr, unit, info_file, /get_lun

   while (eof(unit) EQ 0) do begin
      readf, unit, string
      if strmid(string, 0, 5) EQ 'TTYPE' then begin ;; TTYPE 
         if n_elements(descriptions) EQ 0 then descriptions = strmid(string, 32) else descriptions = [descriptions, strmid(string, 32)]
      endif 

      if strmid(string, 0, 7) EQ 'COMMENT' then begin  ;; Comments
         if n_elements(comments) EQ 0 then comments = strmid(string, 9) else comments = [comments, strmid(string, 9)]
      endif 

      i = i+1
   endwhile
   free_lun, unit

   if n_elements(descriptions) EQ n_ttype then begin
      for i=0, n_ttype-1 do sxaddpar, hdr, 'TTYPE'+strc(i+1), ttype_vals[i], descriptions[i]
   endif else print, 'Number of TTYPE entries in the hdr does not match the number in '+info_file

   if n_elements(comments) GT 0 then begin
      for i=0, n_elements(comments)-1 do sxaddpar, hdr, 'COMMENT', comments[i], BEFORE = 'DATE'
   endif

endelse

end
