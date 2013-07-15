pro pickoption,files

   for i=0,n_elements(files)-1 do begin
      obj=mrdfits(files[i],1,/silent)
;      if plate eq 6870 
      nobj=n_elements(obj)
      meanthrupt = total(obj.thrupt,2)/nobj
      diff = obj.thrupt-meanthrupt#(fltarr(nobj)+1)
      rms  = sqrt(total(diff^2,2)/nobj)
      err = rms/sqrt(nobj-1.)

      fracrms = rms/meanthrupt
      fracerr = err/meanthrupt
      t1=where(obj[0].loglam gt alog10(3700.) and obj[0].loglam lt alog10(3800.))
      t2=where(obj[0].loglam gt alog10(5000.) and obj[0].loglam lt alog10(5100.))
      t3=where(obj[0].loglam gt alog10(6600.) and obj[0].loglam lt alog10(6700.))
      t4=where(obj[0].loglam gt alog10(8900.) and obj[0].loglam lt alog10(9000.))

      print,files[i]
;      print,median(fracrms[t1]),median(fracrms[t2]),median(fracrms[t3]),median(fracrms[t4])
      f1 = median(meanthrupt[t1])
      f2 = median(meanthrupt[t2])
      f3 = median(meanthrupt[t3])
      f4 = median(meanthrupt[t4])
      print,f1/f2,f3/f2,f4/f2,format='(3f7.3)'
stop
      
      print,median(fracerr[t1]),median(fracerr[t2]),median(fracerr[t3]),median(fracerr[t4])
   endfor

end
