;+
; function mlcalifainterp_ivar
;
; As mlcalifainterp, but properly handles inverse variance as well.
;
; Function interpolates irregularly-sampled data to a regular grid
; Assumes input is in the form of 3 1-dimensional vectors containing
; the x coordinates, y coordinates, and fluxes respectively.
; Algorithm used is a flux-conserving variation of Shepard's method
; adopted from the CALIFA algorithm described by Sanchez et al.
;
; Input:
;      coordX = 1-d vector containing X coordinates of centres of fibers
;              size: [# of fibers * # of dither positions]
;              IN OUTPUT PIXEL COORDINATES
;
;      coordY = 1-d vector containing Y coordinates of centres of fibers
;              size: [# of fibers * # of dither positions]
;
;      F_in = 1-d vector containing flux through each fiber
;             size: [# of fibers * # of dither positions]
;
;      fiber_status = 1-d vector defining fiber status (0=live, 1=dead)
;
;      dim_out = array containing [x, y] dimensions of output image
;
;      rlim = scalar containing boundary limit radius
;             IN OUTPUT PIXEL COORDINATES
;
;      ivar_in = inverse variance corresponding to F_in
;
;      ivar_out = inverse variance corresponding to output image
;
;      SIGMA = scalar defining sigma of gaussian for weight
;              IN OUTPUT PIXEL COORDINATES
;
; Output:
;      2-D double precision array of given dimensions
;
; Author:
;   David R. Law (Dunlap Institute; drlaw@di.utoronto.ca)
;   First written: 05/15/2012
;   Last modified: 04/03/2013
;
; Modification history:
;   v1: 05-May-2012 Maryna Tsybulska
;       Original coding
;   v2: 09-Oct-2012 D. Law
;       Modified function to conform with ml library file format,
;       some minor modifications and bug fixes.
;   V3: 12-Feb-2013 D. Law
;       Temporary fix of flux scaling issue.  STILL IMPERFECT.
;   V4: 03-Apr_2013 D. Law
;       Added parameters to figure out inverse variance as well, subject
;       to same flux scaling uncertainties as the flux itself.
;-

function mlcalifainterp_ivar, coordX, coordY, F_in, fiber_status, dim_out, rlim, sigma, ivar_in, ivar_out

; Dimensions
ntot = (size(F_in, /dimensions))[0]; Number of total samples
all_dim = [dim_out, ntot]; Output X size, Y size, total samples

; Output image:
r_image = dblarr(dim_out)
ivar_out=dblarr(dim_out)

; X and Y output pixel coordinate arrays

arr_xcoord = dindgen(dim_out[0])
arr_ycoord = dindgen(dim_out[1])

; Calculade distances from each fiber to each pixel and calculate weights
; Note: need only r < r_fib
; Store weights in a series of slices for each sample
arr_weights = dblarr(all_dim)

for i=0,ntot-1 do begin
   ; Array defining radii away from an output fiber location
   ; Set by default to rlim+1, or anything greater than rlim
   ; (since r<=rlim is the criterion for including in the calculation)
   arr_radius = replicate(double(rlim+1), dim_out)
   ; Figure out the region of influence in which we actually need to calculate radii
   xmin = coordX[i] - rlim > 0
   xmax = coordX[i] + rlim < (dim_out[0]-1)
   ymin = coordY[i] - rlim > 0
   ymax = coordY[i] + rlim < (dim_out[1]-1)

   ; weird but this works faster than to store size in separate variables
  ; Calculate radii in this region of influence
   dim_c = [n_elements(arr_xcoord[xmin:xmax]),n_elements(arr_ycoord[ymin:ymax])]
   arr_radius[xmin:xmax,ymin:ymax] = sqrt(rebin((arr_xcoord[xmin:xmax]-coordX[i])^2,dim_c) $
                              + rebin(transpose((arr_ycoord[ymin:ymax]-coordY[i])^2),dim_c))
   tocalc = where(arr_radius le rlim)
   arr_weights[tocalc+i*n_elements(arr_radius)] = exp(-0.5/sigma^2*arr_radius[tocalc]^2)
endfor

; Delete dead fibers
if ((size(where(fiber_status ne 0)))[0] ne 0) then arr_weights[*, *, where(fiber_status ne 0)] = 0.

; Figure out the normalization matrix
matr_norm = total(arr_weights, 3)
matr_norm += (matr_norm eq 0) ; we don't want to divide by zero

; Calculate weighted sum.  Must also normalize by a factor quantifying
; how big the droplet size is in terms of output pixels.  This is
; roughly the size of where matr_norm is non-zero, divided by the
; number of fibers.  It's also roughly rlim*rlim*PI, which is
; roughly the number of pixels info gets dropped over (though
; this doesn't account for square edge in round area effects).
; This depends on provided rlim though, but using the actual fiber
; radius gives an incorrect answer???  Approach below is good enough
; for now, but systematics remain at factor of 2 level.
;
; Should I use the max of rlim or rfiber??
for i=0,ntot-1 do begin
   alpha=arr_weights[*,*,i] / matr_norm / (rlim*rlim*!DPI)
   r_image += F_in[i] * alpha 
   ; Combine errors by using formula
   ; sigma^2=Sum(alpha^2*sigma_i^2)
   ; where sigma_i = 1./sqrt(ivar_in)
   ivar_out += alpha*alpha/ivar_in[i]
endfor

; Change ivar_out from variance to inverse variance
ivar_out=1./ivar_out

return, r_image
end
