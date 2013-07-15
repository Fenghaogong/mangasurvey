;; ----------------------------------------------------------------------------------------------------
;; MLFADDRESS.PRO
;;   Convert between 2-D and spiral fiber addressing in hexagonal bundles for MaNGA
;;
;; USAGE:
;;   faddress, input, output
;;
;; ARGUMENTS:
;;   input:   Either a spiral index number (scalar), or a 2D hexagonal coordinate position (a,b)
;;   output:  Returns whichever addressing system was not supplied on input
;;
;; KEYWORDS (OUTPUT):
;;   RING:    Returns the outer most ring number of the given fiber (where RING=1 is for 7-fiber bundle)
;;   CLOCK:   Clocking position of the fiber in the outer ring (0 = x-axis); counter-clockwise
;;   PHI:     Position angle in degrees, counter-clockwise, with respect to x-axis (a-axis)
;;
;; K. Bundy 9/2/12

;; ----------------------------------------------------------------------------------------------------
;; Define AB_RING function:
;;   For a given ring number, n -- where n=1 corresponds to a 7-fiber bundle -- determine (a,b) coordinates around the ring
;;   Returns a vector of [a,b] coordinates
;;   Fails for n > 10
function ab_ring, n
  if n GT 10 then begin
     print, 'Ring addressing approximation fails for n > 10'
     return, -1
  endif

  adj = 0.1*n < 0.98
  a_ring = -n>round((n + adj)*sin((findgen(6*n)+2*n)/(6*n)*2*!pi))<n
  b_ring = -n>round((n + adj)*sin((findgen(6*n))/(6*n)*2*!pi))<n

  return, transpose([[a_ring], [b_ring]])
end

;; ----------------------------------------------------------------------------------------------------
;; Main Procedure
;; ----------------------------------------------------------------------------------------------------
pro mlfaddress, input, output, RING = RING, CLOCK = CLOCK, PHI = PHI

ndim = size(input, /dimensions)  ;; Determine whether input is spiral or 2D coordinate based on dimensions of input

nmax = 10
bundle_sizes = 1 + total(findgen(nmax+1)*6,/cum)  ;; Define Possible (total) bundle sizes up to nmax rings

;; -- For converting from spiral index to 2D coordinates --
if ndim NE 2 then begin
   wfull = where( (input-bundle_sizes) GT 0)  ;; Which bundle sizes are filled
   nfull = max(wfull)                         ;; The largest bundle size filled

   fring = input - bundle_sizes[nfull]  ;; The spiral distance on the outer most ring

   n = nfull+1 ;; The rank (number) of the outer most ring
   RING = n

   output = (ab_ring(n))[*, fring]  

   print, output

   ;; plot, a_ring, ps=-8, yr=(n+1)*[-1,1] 
   ;; oplot, b_ring, ps=-8, col=!red
   ;; wshow

;; -- For converting from 2D coordinates to spiral index --
endif else begin
   a = input[0]
   b = input[1]

   if a*b LT 0 then nn = max(abs([a,b])) else nn = total(abs([a,b]))  ;; Determine what ring we're on (nn)

   RING = nn

   ab_ring = ab_ring(nn)  ;; Get coordinates for this ring

   wmatch = where(ab_ring[0,*] EQ a AND ab_ring[1,*] EQ b, nmatch)  ;; Match to the list of coordinates

   clock = wmatch  ;; Clocking position of the fiber in the outer ring (0 = x-axis); counter-clockwise

   phi = float(clock)/(6*Ring) * 360  ;; Position angle, counter-clockwise, with respect to x-axis (a-axis)

   output = wmatch + bundle_sizes[nn-1]
   print, output
endelse

end
