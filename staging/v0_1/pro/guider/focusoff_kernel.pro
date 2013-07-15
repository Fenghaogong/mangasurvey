function focusoff_kernel,offset

   siz = fix(offset/5*1.2)
   siz = (siz/2)*2+1
   kernel = fltarr(siz,siz)
   xx = findgen(siz)-siz/2
   x2 = xx#(fltarr(siz)+1)
   y2 = (fltarr(siz)+1)#xx
   rr = sqrt(x2*x2+y2*y2)
   t = where(rr ge (offset/10.)/2. and rr le (offset/5.)/2.)
   kernel[t] = 1
   return,kernel
end
