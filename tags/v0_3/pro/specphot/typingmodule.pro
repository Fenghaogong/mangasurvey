function typingmodule, objflux,loglam,objivar,dispimg,iphoto,observation,plugmap,kindx=kindx,modflux=modflux

   spectroid=1
   nfiber = 560
   plateid = observation.plate
   maxmjd = observation.mjd

   npix = (size(objflux,/dimen))[0]
   ;----------
   ; For each star, find the best-fit model.

   !p.multi = [0,2,3]
   nphoto = n_elements(iphoto)
   modflux = 0 * objflux
   for ip=0L, nphoto-1 do begin
      thisfiber = iphoto[ip] + 1 + nfiber * (spectroid[0] - 1)
      splog, prelog='Fiber '+string(thisfiber,format='(I4)')

      plottitle = 'PLATE=' + string(plateid[0], format='(i4.4)') $
       + ' MJD=' + string(maxmjd, format='(i5.5)') $
       + ' Spectro-Photo Star' $
       + ' Fiber ' + strtrim(thisfiber,2)

      ; Find the best-fit model -- evaluated for each exposure [NPIX,NEXP]
      thismodel = spflux_bestmodel(loglam[*,*,ip], objflux[*,*,ip], $
       objivar[*,*,ip], dispimg[*,*,ip], kindx=thisindx, plottitle=plottitle)

      ; Also evaluate this model over a big wavelength range [3006,10960] Ang.
      tmploglam = 3.4780d0 + lindgen(5620) * 1.d-4
      tmpdispimg = 0 * tmploglam + 1.0 ; arbitrarily select this resolution
      tmpdispimg = interpol(dispimg[*,*,ip],loglam[*,*,ip],tmploglam)
      tmpflux = spflux_read_kurucz(tmploglam, tmpdispimg, $
       iselect=thisindx.imodel)

      ; The returned models are redshifted, but not fluxed or
      ; reddened.  Do that now...  we compare data vs. model reddened.
      extcurve1 = ext_odonnell(10.^loglam[*,*,ip], 3.1)
      thismodel = thismodel $
       * 10.^(-extcurve1 * 3.1 * plugmap[iphoto[ip]].sfd_ebv / 2.5)
      extcurve2 = ext_odonnell(10.^tmploglam, 3.1)
      tmpflux = tmpflux $
       * 10.^(-extcurve2 * 3.1 * plugmap[iphoto[ip]].sfd_ebv / 2.5)

      ; Now integrate the apparent magnitude for this spectrum,
      ; The units of FTHRU are such that m = -2.5*alog10(FTHRU) + (48.6-2.5*17)
      ; Note that these computed magnitudes, THISMAG, should be equivalent
      ; to THISINDX.MAG in the case of no reddening.
      wavevec = 10.d0^tmploglam
      flambda2fnu = wavevec^2 / 2.99792e18
      fthru = filter_thru(tmpflux * flambda2fnu, waveimg=wavevec, /toair)
      thismag = -2.5 * alog10(fthru) - (48.6-2.5*17)

      ; Compute SCALEFAC = (plugmap flux) / (uncalibrated flux)
      if (tag_exist(plugmap, 'CALIBFLUX')) then begin
         scalefac = plugmap[iphoto[ip]].calibflux[2] $
          / 10.^((22.5-thismag[2])/2.5)
         ; Reject this star if we don't know its flux.
         if (plugmap[iphoto[ip]].calibflux[2] LE 0) then begin
            splog, 'Warning: Rejecting std star in fiber = ', $
             iphoto[ip] + 1 + nfiber * (spectroid[0] - 1), $
             ' with unknown calibObj flux'
            qfinal[ip] = 0
         endif
      endif else begin
         splog, 'WARNING: No CALIBFLUX for zero-pointing the fluxes'
         scalefac = 10.^((thismag[2] - plugmap[iphoto[ip]].mag[2])/2.5)
      endelse
      thismodel = thismodel * scalefac

      modflux[*,*,ip] = thismodel
      if (ip EQ 0) then kindx = replicate( create_struct( $
       'PLATE', 0L, $
       'MJD', 0L, $
       'FIBERID', 0L, $
       'QGOOD', 0, $
       thisindx, $
       'MODELFLUX', fltarr(npix)), $
       nphoto)
      copy_struct_inx, thisindx, kindx, index_to=ip
      kindx[ip].plate = plateid[0]
      kindx[ip].mjd = maxmjd
      kindx[ip].fiberid = thisfiber
      splog, prelog=''
   endfor
   !p.multi = 0
   return,0

end
