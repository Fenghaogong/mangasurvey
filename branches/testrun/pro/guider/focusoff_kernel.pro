;+
;  build a 2-d kernel to be convolved to compute the actual PSF at a given 
;    focus offset
; Input:  offset --- focus offset in micron
;          pixel  --- pixel scale of the kernel in micron
; Output:  2-d kernel 
;-
function focusoff_kernel,offset,pixel=pixel,mkernel=mkernel

   absoff = abs(offset)
   siz = ceil(absoff/5./pixel+1)
   siz = ceil(siz/2.0)*2+1
   kernel = fltarr(siz,siz)

   magnification = 20 > ceil(pixel/(absoff/20.))
   msiz = fix(siz*magnification)
   mkernel = fltarr(msiz,msiz)
   magkernel = fltarr(msiz,msiz)
   xm = findgen(msiz)-msiz/2
   x2 = xm#(fltarr(msiz)+1)
   y2 = (fltarr(msiz)+1)#xm
   rr = sqrt(x2*x2+y2*y2)*pixel/magnification
   t = where(rr ge (absoff/10.)/2. and rr le (absoff/5.)/2.)
   mkernel[t] = 1.0/(magnification^2)
   kernel = rebin(mkernel,siz,siz)

   return,kernel/total(kernel)
end
