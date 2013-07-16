function svdfcn,dataindex, N
common svdinput, x_svd
  return,x_svd[*,dataindex]
end

function doublegau, coords, para,amps=amps
common subimageblock, stamp,errstamp,nx
common svdinput, x_svd 
   xc = para[0]
   yc = para[1]
   sig1sq = para[2]
   sig2sq = para[3]
   xx = (coords mod nx)-xc
   yy = (fix(coords/nx))-yc
   ny = n_elements(xx)/nx
   xx = reform(xx,nx,ny)
   yy = reform(yy,nx,ny)

;   multiply = 10
   x1 = xx[0:nx-1,0]
   y1 = reform(yy[0,0:ny-1],ny)

   xb = [x1,x1[nx-1]+1]
   yb = [y1,y1[ny-1]+1]

   tx1 = (xb-0.5)/sqrt(2.*sig1sq)
   ty1 = (yb-0.5)/sqrt(2.*sig1sq)
   ex1 = sqrt(!pi)/2.*sqrt(2.*sig1sq)*erf(tx1)
   ey1 = sqrt(!pi)/2.*sqrt(2.*sig1sq)*erf(ty1)
   exp1_x = ex1[1:nx]-ex1[0:nx-1]
   exp1_y = ey1[1:ny]-ey1[0:ny-1]

   tx2 = (xb-0.5)/sqrt(2.*sig2sq)
   ty2 = (yb-0.5)/sqrt(2.*sig2sq)
   ex2 = sqrt(!pi)/2.*sqrt(2.*sig2sq)*erf(tx2)
   ey2 = sqrt(!pi)/2.*sqrt(2.*sig2sq)*erf(ty2)
   exp2_x = ex2[1:nx]-ex2[0:nx-1]
   exp2_y = ey2[1:ny]-ey2[0:ny-1]

   templ1 = exp1_x # exp1_y
   templ2 = exp2_x # exp2_y
   templ = [[reform(templ1,nx*ny)],[reform(templ2,nx*ny)]]
   
   valid = where(errstamp lt 1.d29,nvalid)

   amps = myregress(transpose(templ[valid,*]),stamp[valid],measure_errors=errstamp[valid],sigma=error,const=const,chisq=chisq,status=regress_status)
   if regress_status ne 0 then begin 
      print,'Regress fail, using svdfit.'
      x_svd = transpose(templ[valid,*])
      amps = svdfit(indgen(nvalid),stamp[valid],2,measure_errors=errstamp[valid],function_name='svdfcn')
   endif
   amps = reform(amps)
   fit = templ#amps
   return,fit
end

