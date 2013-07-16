pro manga_makeccdspec_part1

galflux=readfits('./inputspec.fits')
zsize=(size(galflux))[1]
inpwave=readfits('./inputwave.fits')

; Read in BOSS calibration file
  calibfile=strcompress(getenv('MANGACORE_DIR')+'/'+getenv('MANGACORE_VER')+'calib/calibmatrix.fits',/remove_all)
  calibmatrix=readfits(calibfile)
  wave=calibmatrix[0,*]
  nwave=(size(calibmatrix))[2]
; calibmatrix[0,*]=wavelengths
; calibmatrix[1,*]=throughput
; calibmatrix[2,*]=erg flux - photons conversion
; calibmatrix[3,*]=sky flux in input photons

; Put input values on proper wavelength grid
tempspec=gaussfold(inpwave,galflux,2.51)
tempspec2=interpol(tempspec,inpwave,wave)
  zerolow=where(wave lt inpwave[0])
  zerohigh=where(wave gt inpwave[zsize-1])
  if (zerolow[0] ne -1) then tempspec2[zerolow]=0.
  if (zerohigh[0] ne -1) then tempspec2[zerohigh]=0.
galflux=tempspec2

; Steven's code assumes input in units of 1e-17 erg/s/cm2/Ang
; Make everything copies of the same thing?
skyonly=calibmatrix[3,*]*calibmatrix[1,*]/calibmatrix[2,*]
totalflux=galflux+skyonly

nlines=1500

OutputBlock=fltarr(nwave,nlines+1)
OutputBlock[*,0]=wave

for i=1,nlines do begin
  OutputBlock[*,i]=totalflux[*]
  if (i mod 30 eq 0) then OutputBlock[*,i]=skyonly[*]
endfor

writefits,'block_manga.fits',OutputBlock

return
end
