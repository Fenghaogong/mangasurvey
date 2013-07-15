forward_function DRLCosmo_DA
forward_function DRLCosmo_DL

; Make MaNGA input files from CALIFA files
pro mkcalifa_ngc2916

cube=readfits('NGC2916.V1200.rscube.fits.gz')
wave=findgen(1701)*0.7+3650.

; Pixels are each 1'' in size
; Default was 1e-16 units, covert to 1e-17
;cube=10.*cube
writefits,'ngc2916crap.fits',cube*10.
stop
redshift=0.
ratio=100.
i=0
; Originally 24''=Re, z=0.013, 127 bundle radius is 15''
; NOT WELL MATCHED, place at proper distance for Re=15''
while (ratio gt 15./24.) do begin
  i=i+1
  redshift=0.013
  DLum0=DRLCosmo_DL(0.732,redshift,0.238,0.762,0.0); Angular diameter distance in Mpc
  DAng0=DRLCosmo_DA(0.732,redshift,0.238,0.762,0.0); Angular diameter distance in Mpc

  redshift=0.013+i*0.001
  DLum=DRLCosmo_DL(0.732,redshift,0.238,0.762,0.0); Angular diameter distance in Mpc
  DAng=DRLCosmo_DA(0.732,redshift,0.238,0.762,0.0); Angular diameter distance in Mpc

  ratio=DAng0/DAng
  print,ratio
;  break
endwhile

print,'redshift=',redshift
  DLum0=DRLCosmo_DL(0.732,0.013,0.238,0.762,0.0); Angular diameter distance in Mpc
  DAng0=DRLCosmo_DA(0.732,0.013,0.238,0.762,0.0); Angular diameter distance in Mpc
  DLum=DRLCosmo_DL(0.732,redshift,0.238,0.762,0.0); Angular diameter distance in Mpc
  DAng=DRLCosmo_DA(0.732,redshift,0.238,0.762,0.0); Angular diameter distance in Mpc

; Scale flux
cube=cube*10.*(DLum0/DLum)*(DLum0/DLum)

; Pixel size
angscale=1.0*DAng0/DAng

writefits,'ngc2916_hires_19input.fits',cube
writefits,'ngc2916_hires_inpwave.fits',wave
print,'Pixel size is ',angscale,' arcsec/pixel'



return
end
