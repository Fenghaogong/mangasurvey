; program manga_quicklook
;
; Interactive program for viewing spectra of specific
; fibers from a MaNGA exposure (SOS or full).
;
; Author:
;   David R. Law (Dunlap Institute; drlaw@di.utoronto.ca)
;   First written: 12/19/2012
;   Last modified: 01/21/2013
;
; REVISION HISTORY:
;   v1.0: 19-Dec-2012  D. Law
;       First written; spectrum display interaction based on
;       continuum-fitting code by George Becker.
;   v1.1: 28-Dec-2012 D. Law
;       Fixed mouse click location errors with some window sizes.
;       Auto-fits windows to screen real estate.
;       Keep same spectrum range for new spectra.
;   v1.2: 09-Jan-2013 D. Law
;       Auto-finds wset and slitmap files, handles reading spectra
;       with new slitmap files.  Handles full and SOS pipeline data.
;   v1.3: 16-Jan-2013 D. Law
;       Modified for new directory structure
;   v1.4: 21-Jan-2013 D. Law
;       Fixed FITS header bug

Pro win_initialize

common cf_state, stateS, stateBM, statePL
common cf_data, bmap, bhdr, smap, shdr, xsS, ysS, xsBM, ysBM, xsPL, ysPL, $
                b1sci,r1sci,b1wave,r1wave,flux,wave,binflux,binwave,index, $
                indextype,bmindex,cenra,cendec,filepath,plate,mjd,scalemm, $
                detector

device,get_screen_size=screen
;stateBM.plot_xsz = stateBM.plot_xsz  < (screen(0)-100)
;stateBM.plot_ysz = stateBM.plot_ysz  < (screen(1)-100)
scalemm=60./1000.

stateS = {              $
          plotwindow: 0L, $        ; Plot window number
          plot_xsz: long(screen(0)*2./3.), $       ; Screen x-size of plot window
          plot_ysz: long(screen(1)/3), $        ; Screen y-size of plot window
          ;plot_xsz: 800L, $       ; Screen x-size of plot window
          ;plot_ysz: 400L, $        ; Screen y-size of plot window
          xdispmin: 0d0, $         ; Min display x-value 
          xdispmax: 0d0, $         ; Max display x-value
          ydispmin: 0d0, $         ; Min display x-value 
          ydispmax: 0d0, $         ; Max display x-value
          default_xdispmin: 0d0, $ ; Default min display x-value 
          default_xdispmax: 0d0, $ ; Default max display x-value
          default_ydispmin: 0d0, $ ; Default min display x-value 
          default_ydispmax: 0d0, $ ; Default max display x-value
          binsz: 1L, $             ; Bin size for displayed spectrum
          histogram: 0L, $         ; Flag to plot in histogram mode
          minwave: 0L, $           ; Minimum wavelength for the current order
          maxwave: 0L, $           ; Maximum wavelength for the current order
          xcoord: 0., $            ; Data x-coordinate of mouse click
          ycoord: 0., $            ; Data y-coordinate of mouse click
          order: 0L, $             ; Index of current order
          norders: 0L, $           ; Number of orders
          npix: 0L, $              ; Number of pixels in each order
          buttonbase_xsz: 0L, $    ; Screen x-size of button base
          buttonbase_ysz: 0L, $    ; Screen y-size of button base
          base_id: 0L, $           ; WIDGET IDS
          subbase_id: 0L, $        ;
          buttonbase_id: 0L, $     ;
          plotwindow_id: 0L, $     ; 
          loadexp_button_id: 0L, $ ;
          bluespec_button_id: 0L, $    ;
          redspec_button_id: 0L, $     ;
          savefits_button_id: 0L, $ ;
          quit_button_id: 0L $     ;
        }

stateBM = {              $
          plotwindow: 1L, $        ; Plot window number
          plot_xsz: long(screen(0)/2), $       ; Screen x-size of plot window
          plot_ysz: long(screen(1)/2.), $        ; Screen y-size of plot window
          ;plot_xsz: 700L, $       ; Screen x-size of plot window
          ;plot_ysz: 700L, $        ; Screen y-size of plot window
          xdispmin: 0d0, $         ; Min display x-value 
          xdispmax: 0d0, $         ; Max display x-value
          ydispmin: 0d0, $         ; Min display x-value 
          ydispmax: 0d0, $         ; Max display x-value
          default_xdispmin: 0d0, $ ; Default min display x-value 
          default_xdispmax: 0d0, $ ; Default max display x-value
          default_ydispmin: 0d0, $ ; Default min display x-value 
          default_ydispmax: 0d0, $ ; Default max display x-value
          xcoord: 0., $            ; Data x-coordinate of mouse click
          ycoord: 0., $            ; Data y-coordinate of mouse click
          order: 0L, $             ; Index of current order
          norders: 0L, $           ; Number of orders
          npix: 0L, $              ; Number of pixels in each order
          buttonbase_xsz: 0L, $    ; Screen x-size of button base
          buttonbase_ysz: 0L, $    ; Screen y-size of button base
          base_id: 1L, $           ; WIDGET IDS
          subbase_id: 0L, $        ;
          buttonbase_id: 0L, $     ;
          plotwindow_id: 1L, $     ; 
          quit_button_id: 0L $     ;
        }

statePL = {              $
          plotwindow: 2L, $        ; Plot window number
          plot_xsz: long(screen(0)/2), $       ; Screen x-size of plot window
          plot_ysz: long(screen(1)/2.), $        ; Screen y-size of plot window
          ;plot_xsz: 700L, $       ; Screen x-size of plot window
          ;plot_ysz: 700L, $        ; Screen y-size of plot window
          xdispmin: 0d0, $         ; Min display x-value 
          xdispmax: 0d0, $         ; Max display x-value
          ydispmin: 0d0, $         ; Min display x-value 
          ydispmax: 0d0, $         ; Max display x-value
          default_xdispmin: 0d0, $ ; Default min display x-value 
          default_xdispmax: 0d0, $ ; Default max display x-value
          default_ydispmin: 0d0, $ ; Default min display x-value 
          default_ydispmax: 0d0, $ ; Default max display x-value
          xcoord: 0., $            ; Data x-coordinate of mouse click
          ycoord: 0., $            ; Data y-coordinate of mouse click
          order: 0L, $             ; Index of current order
          norders: 0L, $           ; Number of orders
          npix: 0L, $              ; Number of pixels in each order
          buttonbase_xsz: 0L, $    ; Screen x-size of button base
          buttonbase_ysz: 0L, $    ; Screen y-size of button base
          base_id: 2L, $           ; WIDGET IDS
          subbase_id: 0L, $        ;
          buttonbase_id: 0L, $     ;
          plotwindow_id: 2L, $     ; 
          quit_button_id: 0L $     ;
        }

end

;------------------------------------------------------------------------
; Widget events (main)
;------------------------------------------------------------------------

Pro spec_window_event,event

common cf_state
common cf_data

widget_control,event.id,get_uvalue=eventval

case eventval of 
    'plotwindow' : begin
                      if event.clicks then begin
                         ; Mouse event: act when mouse button is pressed
                         if (event.press gt 0) then mouse_event,event
                      endif else begin
                         ; Keyboard event: act when key is pressed
                         if (event.press gt 0) then keyboard_event,event
                      endelse
                   end
    'loadexp'       : begin
                      readBOSSfiles
                      readslitmap
                      ; Default spectrum to setup with
                      index=0
                      indextype='SINGLE'
                      detector='b1'
                      selectspectrum,detector
                      plotspec
                      plotpl
                   end
    'bluespec'       : begin
                      detector='b1'
                      selectspectrum,detector
                      plotspec
                   end
    'redspec'       : begin
                      detector='r1'
                      selectspectrum,detector
                      plotspec
                   end
    'savefits'       : begin
                      writespectrum
                      print,'Not implemented.'
;                      print,'Wrote to outputspec.fits and outputwave.fits'
                   end
    'quit'       : begin
                      widget_control,event.top,/destroy
                   end
    else         : print,'Unknown event'
endcase

end

Pro bm_window_event,event

common cf_state
common cf_data

widget_control,event.id,get_uvalue=eventval

case eventval of 
    'plotwindow' : begin
                      if event.clicks then begin
                         ; Mouse event: act when mouse button is pressed
                         if (event.press gt 0) then mouse_eventBM,event
                      endif
                   end
    else         : print,'Unknown event'
endcase

end

Pro pl_window_event,event

common cf_state
common cf_data

widget_control,event.id,get_uvalue=eventval

case eventval of 
    'plotwindow' : begin
                      if event.clicks then begin
                         ; Mouse event: act when mouse button is pressed
                         if (event.press gt 0) then mouse_eventPL,event
                      endif
                   end
    else         : print,'Unknown event'
endcase

end

;------------------------------------------------------------------------
; Mouse events
;------------------------------------------------------------------------

Pro mouse_event,event

common cf_state
common cf_data

; Restore scaling information
wset,stateS.plotwindow
!X.S = xsS
!Y.S = ysS

; Convert coordinates from device to data
coords = convert_coord(event.x,event.y,$
                       /device,/to_data)

stateS.xcoord = coords(0)
stateS.ycoord = coords(1)

case event.press of
   1: begin  ; Left
         ; Get new x-range
         print,'Select second x-bound'
         cursor,xselect,yselect,/down,/data
         stateS.xdispmin = min([xselect,stateS.xcoord])
         stateS.xdispmax = max([xselect,stateS.xcoord])
         print,'click1'
      end
   2: begin  ; Center
         ; Get new x/y-range
         print,'Select second bound'
         cursor,xselect,yselect,/down,/data
         stateS.xdispmin = min([xselect,stateS.xcoord])
         stateS.xdispmax = max([xselect,stateS.xcoord])
         stateS.ydispmin = min([yselect,stateS.ycoord])
         stateS.ydispmax = max([yselect,stateS.ycoord])
      end
   4: begin  ; Right
         ; Get new y-range
         print,'Select second y-bound'
         cursor,xselect,yselect,/down,/data
         stateS.ydispmin = min([yselect,stateS.ycoord])
         stateS.ydispmax = max([yselect,stateS.ycoord])
      end
   else: print,'Unknown mouse button'
endcase
;print,'stop1'
plotspec

end

; Figures out which fiber matches in a bundle and puts
; info in the global variable index
pro bm_whichfiber,xselect,yselect

common cf_data

xoff=xselect-bmap.xpmm/scalemm
yoff=yselect-bmap.ypmm/scalemm

roff=sqrt(xoff*xoff+yoff*yoff)
; Point to fiber it fell in, if that fiber is alive
bmindex=where((roff lt 1.0)and(bmap.gbu eq 1))

end

Pro mouse_eventBM,event

common cf_state
common cf_data

; Restore scaling information
wset,stateBM.plotwindow
!X.S = xsBM
!Y.S = ysBM

; Convert coordinates from device to data
coords = convert_coord(event.x,event.y,$
                       /device,/to_data)

xselect=coords(0)
yselect=coords(1)
;print,xselect,yselect

; Find which fiber was closest
bm_whichfiber,xselect,yselect

; If it found a valid fiber, highlight the fiber in green (color=140)
; and plot it in the spectrum window
if bmindex ne -1 then begin
  ;print,bmap[index].fnum
  plotbm
  wset,stateBM.plotwindow
  mldrawcirc,bmap[bmindex].xpmm/scalemm,bmap[bmindex].ypmm/scalemm,2.0,strcompress(string(bmap[bmindex].ise),/remove_all),color=140
  ; Figure out what fiber number on the slit it corresponds to
  stemp=smap[where(smap.ifuname eq indextype)]
  temp=where(stemp.fnum eq bmap[bmindex].ise)
  index=where(smap.fiberid eq stemp[temp].fiberid)

  ; And plot it
  selectspectrum,detector
  plotspec
endif

;print,'stop2'

end

; Figures out which fiber matches on a plate click and puts
; info in the global variable index, with type indextype
pro pl_whichfiber,xselect,yselect

common cf_data

xoff=(xselect-smap.ra)*cos(cendec*!PI/180.)
yoff=yselect-smap.dec

roff=sqrt(xoff*xoff+yoff*yoff)
index=where(roff eq min(roff))

indextype='SINGLE'

; If it found an IFU, just point to one fiber for now
if (size(index))[1] gt 1 then index=index[0]

; If it was an IFU (name starts with 'ma'), set index type to IFU name
if strmid(smap[index].ifuname,0,2) eq 'ma' then indextype=smap[index].ifuname

end

Pro mouse_eventPL,event

common cf_state
common cf_data

; Restore scaling information
wset,statePL.plotwindow
!X.S = xsPL
!Y.S = ysPL

; Convert coordinates from device to data
coords = convert_coord(event.x,event.y,$
                       /device,/to_data)

xselect=coords(0)
yselect=coords(1)
;print,xselect,yselect

; Find which fiber was closest
pl_whichfiber,xselect,yselect
;print,'index=',index
selptx=dblarr(1)
selpty=dblarr(1)
selptx[0]=smap[index].ra
selpty[0]=smap[index].dec

; If valid fiber, which size?  Assign overlap color appropriately
if (smap[index].fsize eq 2.) then pcol=255
if (smap[index].fsize eq 3.) then pcol=150
if (smap[index].fsize eq 5.) then pcol=210
if indextype ne 'SINGLE' then pcol=100

; If it found a valid fiber, highlight the fiber in appropriate color
if index ne -1 then begin
  plotpl
  wset,statePL.plotwindow
  oplot,selptx,selpty,psym=7,color=pcol,thick=4
endif

; If it was a single fiber, plot the spectrum
if indextype eq 'SINGLE' then begin
  selectspectrum,detector
  plotspec
endif

; If it was an IFU, plot the IFU to select fiber
if indextype ne 'SINGLE' then begin
  readbmap
  plotbm
endif

end

;------------------------------------------------------------------------
; Keyboard events
;------------------------------------------------------------------------

Pro keyboard_event,event

common cf_state
common cf_data

; Convert coordinates from device to data
coords = convert_coord(event.x,event.y,$
                       /device,/to_data)
stateS.xcoord = coords(0)
stateS.ycoord = coords(1)

; Determine which key was pressed
if (event.ch ne 0) then begin
   case event.ch of
      8:    key = 'backspace'
      9:    key = 'tab'
      13:   key = 'return'
      27:   key = 'escape'
      32:   key = 'space'
      127:  key = 'delete'
      else: key = string(event.ch)
   endcase
endif else begin
   case event.key of
      5:    key = 'leftarrow'
      6:    key = 'rightarrow'
      7:    key = 'uparrow'
      8:    key = 'downarrow'
      9:    key = 'pageup'
      10:   key = 'pagedown'
      11:   key = 'home'
      12:   key = 'end'
      else: key = 'unknown'
   endcase
endelse

case key of
   ;;;
   ;;; Display keys
   ;;;
   'h': begin   ; Toggle histogram mode
           if (stateS.histogram eq 0) then stateS.histogram = 1 else $
              stateS.histogram = 0
           plotspec
        end
   ']': begin   ; Step up in wavelength
           if (stateS.xdispmax lt stateS.maxwave) then begin
              dwave = stateS.xdispmax - stateS.xdispmin
              stateS.xdispmax = stateS.xdispmax + 0.95*dwave
              stateS.xdispmin = stateS.xdispmax - dwave
              plotspec
           endif
        end
   '[': begin   ; Step down in wavelength
           if (stateS.xdispmin gt stateS.minwave) then begin
              dwave = stateS.xdispmax - stateS.xdispmin
              stateS.xdispmin = stateS.xdispmin - 0.95*dwave
              stateS.xdispmax = stateS.xdispmin + dwave
              plotspec
           endif
        end
   'p': begin   ; Pan to the cursor wavelength position
           dwave = stateS.xdispmax - stateS.xdispmin
           stateS.xdispmin = stateS.xcoord - 0.5*dwave
           stateS.xdispmax = stateS.xcoord + 0.5*dwave
           plotspec
        end
   'i': begin   ; Zoom in at the cursor position (and rescale in y)
           dwave = stateS.xdispmax - stateS.xdispmin
           stateS.xdispmin = stateS.xcoord - 0.25*dwave
           stateS.xdispmax = stateS.xcoord + 0.25*dwave
           display_ls = where(wave(*,stateS.order) ge stateS.xdispmin and $
                              wave(*,stateS.order) le stateS.xdispmax,n_display)
           if (n_display gt 0) then begin
              disp_flux_range = max(flux(display_ls,stateS.order)) - $
                                min(flux(display_ls,stateS.order))
              stateS.ydispmin = 0 < (min(flux(display_ls,stateS.order)) - $
                                    0.05*disp_flux_range)
              stateS.ydispmax = max(flux(display_ls,stateS.order)) + $
                                    0.05*disp_flux_range
           endif
           plotspec
        end
   'o': begin   ; Zoom out at the cursor position (and rescale in y)
           dwave = stateS.xdispmax - stateS.xdispmin
           stateS.xdispmin = stateS.xcoord - 1.0*dwave
           stateS.xdispmax = stateS.xcoord + 1.0*dwave
           display_ls = where(wave(*,stateS.order) ge stateS.xdispmin and $
                              wave(*,stateS.order) le stateS.xdispmax,n_display)
           if (n_display gt 0) then begin
              disp_flux_range = max(flux(display_ls,stateS.order)) - $
                                min(flux(display_ls,stateS.order))
              stateS.ydispmin = 0 < (min(flux(display_ls,stateS.order)) - $
                                    0.05*disp_flux_range)
              stateS.ydispmax = max(flux(display_ls,stateS.order)) + $
                                    0.05*disp_flux_range
           endif
           plotspec
        end
   'w': begin   ; Reset plot range
           stateS.xdispmin = stateS.default_xdispmin
           stateS.xdispmax = stateS.default_xdispmax
              stateS.ydispmin = stateS.default_ydispmin
              stateS.ydispmax = stateS.default_ydispmax
           plotspec
        end

   '1': cf_bin_arrays,1
   '2': cf_bin_arrays,2
   '3': cf_bin_arrays,3
   '4': cf_bin_arrays,4
   '5': cf_bin_arrays,5
   '6': cf_bin_arrays,6
   '7': cf_bin_arrays,7
   '8': cf_bin_arrays,8
   '9': cf_bin_arrays,9
   '0': cf_bin_arrays,10
   'tab': print,'Tab key pressed.  Be careful of highlighted buttons!'
   ;;; 

   's': begin   ; Print current cursor wavelength
      print,'Cursor wavelength: ',stateS.xcoord
        end
   ;;;
   ;;; Other keys
   ;;;
   '?': begin
           print,'Display Commands'
           print,'     left click : Set x range'
           print,'    right click : Set y range'
           print,'   center click : Set x and y ranges'
           print,'              [ : Pan left'
           print,'              ] : Pan right'
           print,'              p : Pan to the cursor position'
           print,'              i : Zoom in'
           print,'              o : Zoom out'
           print,'              w : Display whole order'
           print,'            0-9 : Bin spectrum'
           print,'              h : Plot in histogram mode (toggle)'
           print,'              s : Print cursor wavelength'
        end
   else: print,'No command for '+key+' key.  Press ? for menu.'
endcase

end

;------------------------------------------------------------------------
; Resize window events
;------------------------------------------------------------------------

Pro resize_spec_window,event

common cf_state

; Resize plotting window
widget_control,stateS.plotwindow_id,$
               draw_xsize=((event.x - 6) > stateS.buttonbase_xsz),$
               draw_ysize=((event.y - stateS.buttonbase_ysz - 9) > 1)

; Redraw plot
plotspec

end

Pro resize_bm_window,event

common cf_state

; Resize plotting window
widget_control,stateBM.plotwindow_id,$
               draw_xsize=((event.x - 6) > stateBM.buttonbase_xsz),$
               draw_ysize=((event.y - stateBM.buttonbase_ysz - 9) > 1)

; Redraw plot
plotbm

end

Pro resize_pl_window,event

common cf_state

; Resize plotting window
widget_control,statePL.plotwindow_id,$
               draw_xsize=((event.x - 6) > statePL.buttonbase_xsz),$
               draw_ysize=((event.y - statePL.buttonbase_ysz - 9) > 1)

; Redraw plot
plotpl

end

;------------------------------------------------------------------------
; Initialize widgets
;------------------------------------------------------------------------

Pro create_spec_window,group=group

common cf_state

; Top level base (for resizing)

stateS.base_id = widget_base(title='Spec',tlb_size_events=1)

; Create base to actually hold all the widgets and direct events
; to the correct proceedure.

stateS.subbase_id = widget_base(stateS.base_id,/column,/base_align_right,$
                               event_pro='spec_window_event')

; Adjust plot windows to fit screen, if necessary
device,get_screen_size=screen
;stateS.plot_xsz = screen(0)
;stateS.plot_ysz = screen(1)/4
;stateS.plot_xsz = stateS.plot_xsz  < (screen(0)-100)
;stateS.plot_ysz = stateS.plot_ysz  < (screen(1)-100)

stateS.plotwindow_id = widget_draw(stateS.subbase_id,$
                                  scr_xsize=stateS.plot_xsz,$
                                  scr_ysize=stateS.plot_ysz,$
                                  uvalue='plotwindow',$
                                  keyboard_events=1,$
                                  /button_events)

stateS.buttonbase_id = widget_base(stateS.subbase_id,/row)

stateS.loadexp_button_id = widget_button(stateS.buttonbase_id,value='Load Exposure',$
                                     uvalue='loadexp')

stateS.bluespec_button_id = widget_button(stateS.buttonbase_id,value='Blue Spectrum',$
                                     uvalue='bluespec')

stateS.redspec_button_id = widget_button(stateS.buttonbase_id,value='Red Spectrum',$
                                     uvalue='redspec')

stateS.savefits_button_id = widget_button(stateS.buttonbase_id,value='Save FITS',$
                                     uvalue='savefits')

stateS.quit_button_id = widget_button(stateS.buttonbase_id,value='Quit Program',$
                                     uvalue='quit')

; Get size of button row
buttonbase_geom = widget_info(stateS.buttonbase_id,/geometry)
stateS.buttonbase_xsz = buttonbase_geom.xsize
stateS.buttonbase_ysz = buttonbase_geom.ysize

; Realize widgets
widget_control,stateS.base_id,/realize

; Determine window numbers
widget_control,stateS.plotwindow_id, get_value = tmp_value
stateS.plotwindow = tmp_value

end

Pro create_bm_window,group=group

common cf_state

device,get_screen_size=screen

; Top level base (for resizing)

stateBM.base_id = widget_base(title='BM',tlb_size_events=1,xoffset=screen[0]/2,yoffset=screen[1]/2)

; Create base to actually hold all the widgets and direct events
; to the correct proceedure.

stateBM.subbase_id = widget_base(stateBM.base_id,/column,/base_align_right,$
                               event_pro='bm_window_event')

; Adjust plot windows to fit screen, if necessary
stateBM.plot_xsz = stateBM.plot_xsz  < (screen(0)-100)
stateBM.plot_ysz = stateBM.plot_ysz  < (screen(1)-100)

stateBM.plotwindow_id = widget_draw(stateBM.subbase_id,$
                                  scr_xsize=stateBM.plot_xsz,$
                                  scr_ysize=stateBM.plot_ysz,$
                                  uvalue='plotwindow',$
                                  keyboard_events=1,$
                                  /button_events)

; Realize widgets
widget_control,stateBM.base_id,/realize

; Determine window numbers
widget_control,stateBM.plotwindow_id, get_value = tmp_value
stateBM.plotwindow = tmp_value

end

Pro create_pl_window,group=group

common cf_state

device,get_screen_size=screen

; Top level base (for resizing)

statePL.base_id = widget_base(title='PL',tlb_size_events=1,yoffset=screen[1]/2)

; Create base to actually hold all the widgets and direct events
; to the correct proceedure.

statePL.subbase_id = widget_base(statePL.base_id,/column,/base_align_right,$
                               event_pro='pl_window_event')

; Adjust plot windows to fit screen, if necessary
statePL.plot_xsz = statePL.plot_xsz  < (screen(0)-100)
statePL.plot_ysz = statePL.plot_ysz  < (screen(1)-100)

statePL.plotwindow_id = widget_draw(statePL.subbase_id,$
                                  scr_xsize=statePL.plot_xsz,$
                                  scr_ysize=statePL.plot_ysz,$
                                  uvalue='plotwindow',$
                                  keyboard_events=1,$
                                  /button_events)

; Realize widgets
widget_control,statePL.base_id,/realize

; Determine window numbers
widget_control,statePL.plotwindow_id, get_value = tmp_value
statePL.plotwindow = tmp_value

end

;------------------------------------------------------------------------
; Draw plot
;------------------------------------------------------------------------

Pro plotspec

common cf_state
common cf_data

wset,stateS.plotwindow

fit_thick = 2.0

; Define circle plot symbol
plotsym,0,thick=fit_thick,/fill

; Histogram mode?
if (stateS.histogram eq 1) then  data_psym = 10 else data_psym = 0

   ; Plot original arrays
   ; Zero line
   plot,[-1e9,1e9],[0,0],linestyle=1,$
                  xrange=[stateS.xdispmin,stateS.xdispmax],/xs,$
                  yrange=[stateS.ydispmin,stateS.ydispmax],/ys,$
                  title='Row '+strtrim(stateS.order,2)
   ; Flux
   oplot,binwave,binflux,psym=data_psym

; Label it
xylab_x=(stateS.xdispmax-stateS.xdispmin)*1./5.+stateS.xdispmin
xylab_y=(stateS.ydispmax-stateS.ydispmin)*4./5.+stateS.ydispmin
xylab_y2=(stateS.ydispmax-stateS.ydispmin)*3.8/5.+stateS.ydispmin
xylab_y3=(stateS.ydispmax-stateS.ydispmin)*3.6/5.+stateS.ydispmin
xyouts,xylab_x,xylab_y,strcompress('FIBERID '+string(smap[index].fiberid)),alignment=0.5,charsize=2
xyouts,xylab_x,xylab_y2,smap[index].ifuname,alignment=0.5,charsize=2
xyouts,xylab_x,xylab_y3,strcompress('FNUM '+string(smap[index].fnum)),alignment=0.5,charsize=2

; Save scaling information
xsS=!X.S
ysS=!Y.S

end

Pro plotbm

common cf_state
common cf_data

wset,stateBM.plotwindow
mlplotbm,bmap,bhdr

; Save scaling information
xsBM=!X.S
ysBM=!Y.S

end

Pro plotpl

common cf_state
common cf_data

wset,statePL.plotwindow
mlplotplate,smap,shdr

; Save scaling information
xsPL=!X.S
ysPL=!Y.S

end

;------------------------------------------------------------------------
; Bin arrays
;------------------------------------------------------------------------

Pro cf_bin_arrays,binsz

common cf_state
common cf_data

binsz = long(binsz)

stateS.binsz = binsz

n_bins = floor(1.*stateS.npix / binsz)

if (binsz eq 1) then begin
   binflux  = flux
   binwave  = wave
endif else begin
   binflux  = fltarr(n_bins,stateS.norders)
   binwave  = dblarr(n_bins,stateS.norders)
   i = 0L
   for i=0,n_bins-1 do binflux(i,*) = $
                           total(flux(i*binsz:(i+1)*binsz-1,*),1)
   for i=0,n_bins-1 do binwave(i,*) = $
                           total(wave(i*binsz:(i+1)*binsz-1,*),1)
   binflux  = binflux / binsz
   binwave  = binwave / binsz
endelse

plotspec

end

; Read in data from a file
pro readBOSSfiles

common cf_state
common cf_data

b1fluxfile = dialog_pickfile( title='Select b1 flux file', $
                                     filter=['sci-*b*.fits;spFrame-b1*'], $
                                     get_path=filepath, $
                                     /fix_filter, $
                                     /MUST_EXIST)

; Parse the filename to get information

; Set input flavor.
; NONE=unknown, SOS is SOS, FULL is full BOSS pipeline
; This sets how to read in ancillary info like wavelength solution
; and red frames
flavor='NONE'
temp=strpos(b1fluxfile,'sci-')
if temp ne -1 then flavor='SOS'
temp=strpos(b1fluxfile,'spFrame')
if temp ne -1 then flavor='FULL'

if flavor eq 'NONE' then begin
  print,'Input file flavor unrecognized- exit!'
  exit
endif

head=headfits(b1fluxfile)

; Process SOS files
if flavor eq 'SOS' then begin
  r1fluxfile=mlstrreplace(b1fluxfile,'-b1-','-r1-')

  b1wsetfile=file_search(filepath+'wset*b1*.fits')
  if (size(b1wsetfile))[1] ne 1 then begin
    ; If no wset file, error
    if (size(b1wsetfile))[0] eq 0 then begin
      print,'ERROR: No wset*b1.fits file in input data directory'
      exit
    endif
    ; If more than one wset file, use the first
    if (size(b1wsetfile))[0] ne 0 then b1wsetfile=b1wsetfile[0]  
  endif

  r1wsetfile=mlstrreplace(b1wsetfile,'-b1','-r1')

  ; Read all data from files
  b1sci=readfits(b1fluxfile)
  b1wset=mrdfits(b1wsetfile,1); Array of legendre coefficients
  r1sci=readfits(r1fluxfile)
  r1wset=mrdfits(r1wsetfile,1); Array of legendre coefficients

  ; Convert from legendre coefficient to CCD wavelength solution
  traceset2xy,b1wset,dummy,b1loglam; Convert from coefficients to log wavelengths
  b1wave=10.^b1loglam; Array of wavelengths for each fiber
  traceset2xy,r1wset,dummy,r1loglam; Convert from coefficients to log wavelengths
  r1wave=10.^r1loglam; Array of wavelengths for each fiber
endif

; Process full pipeline files
if flavor eq 'FULL' then begin
  r1fluxfile=mlstrreplace(b1fluxfile,'-b1-','-r1-')

  ; Read science data
  b1sci=mrdfits(b1fluxfile,9)
  r1sci=mrdfits(r1fluxfile,9)

  ; Read wavelength solution
  b1wset=mrdfits(b1fluxfile,3)
  r1wset=mrdfits(r1fluxfile,3)
  ; Convert from legendre coefficient to CCD wavelength solution
  traceset2xy,b1wset,dummy,b1loglam; Convert from coefficients to log wavelengths
  b1wave=10.^b1loglam; Array of wavelengths for each fiber
  traceset2xy,r1wset,dummy,r1loglam; Convert from coefficients to log wavelengths
  r1wave=10.^r1loglam; Array of wavelengths for each fiber
endif

mjd=strcompress(fxpar(head,'MJD'),/remove_all)
plate=strcompress(fxpar(head,'PLATEID'),/remove_all)

end

; Parses the current data file and puts the spectrum of the
; selected fiber into a vector
pro selectspectrum,ccd

common cf_state
common cf_data

if ccd eq 'b1' then begin
  flux=b1sci[*,index]
  wave=b1wave[*,index]
endif

if ccd eq 'r1' then begin
  flux=r1sci[*,index]
  wave=r1wave[*,index]
endif

; Determine number of pixels, ranges
sz = size(flux)
stateS.npix    = sz(1)
stateS.minwave = min(wave)
stateS.maxwave = max(wave)
wave_range = stateS.maxwave - stateS.minwave
flux_range = max(flux) - min(flux)
stateS.default_xdispmin = min(wave) - 0.02*wave_range
stateS.default_xdispmax = max(wave) + 0.02*wave_range
stateS.default_ydispmin = 0 < (min(flux) - 0.05*flux_range)
stateS.default_ydispmax = max(flux) + 0.05*flux_range

; Set plot limits to the default if no limits already in buffer
if stateS.xdispmin eq 0. then begin
  stateS.xdispmin = stateS.default_xdispmin
  stateS.xdispmax = stateS.default_xdispmax

  stateS.ydispmin = stateS.default_ydispmin
  stateS.ydispmax = stateS.default_ydispmax
endif

; Set the initial binned arrays to be the input arrays
binflux  = flux
binwave  = wave

end

; Write the current displayed spectrum to a file
; Currently hard-set to outputspec.fits and outputwave.fits
pro writespectrum

common cf_data

;writefits,'outputspec.fits',flux
;writefits,'outputwave.fits',wave

end

pro readbmap

common cf_data

bmfile=strcompress(getenv('MANGACORE_DIR')+'/'+getenv('MANGACORE_VER')+'/metrology/56280/v1.1/'+indextype+'_bmp_v1.1.par',/remove_all)
bmap=yanny_readone(bmfile,'BUNDLEMAP',hdr=bhdr)

end

pro readslitmap

common cf_data

slitmapfile=strcompress(getenv('MANGACORE_DIR')+'/'+getenv('MANGACORE_VER')+'/slitmaps/slitmap-'+string(plate)+'-'+string(mjd)+'.par',/remove_all)
print,slitmapfile
smap=yanny_readone(slitmapfile,'SMAP',hdr=shdr)

cenra=smap[0].cenra
cendec=smap[0].cendec

end

;------------------------------------------------------------------------
; Main program - see beginning of file for header
;------------------------------------------------------------------------

Pro manga_quicklook

common cf_state
common cf_data

; Initialize common blocks and settings 
win_initialize
print,stateS.npix

readBOSSfiles

readslitmap

; Default spectrum to setup with
index=0
indextype='SINGLE'
detector='b1'
selectspectrum,detector

; Use rainbow color indices
device,decomp=0
loadct,39

; Create a new window if one does not already exist
if not(xregistered('spec_window')) then create_spec_window
if not(xregistered('bm_window')) then create_bm_window
if not(xregistered('pl_window')) then create_pl_window

plotspec

; Plot the plate
plotpl

; Print help reminder
print,"Type '?' for menu."

; Register widgets- leader is stateS
xmanager,'spec_window',stateS.base_id,group_leader=stateS.base_id,$
         event_handler='resize_spec_window',/just_reg
xmanager,'bm_window',stateBM.base_id,group_leader=stateS.base_id,$
         event_handler='resize_bm_window',/just_reg
xmanager,'pl_window',statePL.base_id,group_leader=stateS.base_id,$
         event_handler='resize_pl_window',/just_reg

return
end

