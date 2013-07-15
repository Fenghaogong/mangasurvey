;+
; Purpose: compute the focus offset according to the shape of the focal plane 
;          at different wavelengths and different position on the plate
; Usage:
;   focusoffset,xfocal,yfocal, wavearr, focusarr
; Input:
;   xfocal: x-position on plate in millimeter; it can be a scalar or 
;           a 1-d array of Npos elements but have to match the size of
;           yfocal.
;   yfocal: y-position on plate in millimeter; it can be a scalar or
;           a 1-d array of Npos elements but have to match the size of 
;           xfocal
;   wavearr: desired wavelengths for which the focus offset is to be computed.
;            an array of [Nwave] elements 
;  
; Output: 
;   focusarr: array of focus offset correponding to the wavelength array.
;             If xfocal and yfocal are scalars, focusarr will be a 1-d array 
;             with size Nwave. 
;             If xfocal and yfocal are arrays, then focusarr will be a 2-d 
;             array with the dimension (Nwave, Npos) with each row 
;             corresponding to the offsets for one position.
;-
pro focusoffset, xfocal, yfocal, wavearr, focusarr
common focustable, fldangle, lambda ,surface1,fldht,surface2,foffset,rms

   if n_elements(fldangle) eq 0 then begin
      dir = getenv('MANGAROOT')+'/manga3d/'+getenv('MANGAVER')+'/etc/'
      readcol,dir+'focalplane.txt',fldangle,lambda,surface1,fldht,surface2,foffset,rms
   endif

   nlambda = 5 
   nangle = 7
   fldangle=reform(fldangle,nangle,nlambda)
   lambda=reform(lambda,nangle,nlambda)
   foffset=reform(foffset,nangle,nlambda)


   npos = n_elements(xfocal)
   if n_elements(yfocal) ne n_elements(xfocal) then $
     message,'Error: the size of array xfocal does not match the size of array yfocal'

;  the scale on the focal plane is 3.62730mm/arcmin (Gunn et al. 2006)
   xfocal=float(xfocal)
   yfocal=float(yfocal)
   angle = sqrt(xfocal*xfocal+yfocal*yfocal)/3.62730

   ss = sort(lambda[0,*])
   lambda = lambda[*,ss]
   fldangle = fldangle[*,ss]
   foffset = foffset[*,ss]
   df = fltarr(5,npos) ; focus offset array for 4000, 4600, 5300, 6500, 9000A
   for il=0,4 do begin
     df[il,*] = interpol(foffset[*,il],fldangle[*,il], angle,/quadratic)
   endfor

   nwave = n_elements(wavearr)
   focusarr = fltarr(nwave,npos)
   for i=0,npos-1 do begin
      focusarr[*,i] = interpol(df[*,i],reform(lambda[0,*]),wavearr,/quadratic)
   endfor
end


