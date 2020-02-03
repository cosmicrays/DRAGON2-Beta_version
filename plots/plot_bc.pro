PRO plot_bc, filename, phi, PS=ps
;+
; NAME:
;     PLOT_BC
; PURPOSE:
;     Plot B/C from DRAGON output against HEAO-3 data
;
; CALLING SEQUENCE:
;     plot_bc, filename, phi [, /PS]
;
; INPUTS:
;     filename = The input file
;     phi = modulation potential [GV]
;
; OPTIONAL INPUT:
;     /PS = The plot is saved in BC.eps file
;
; EXAMPLES:
;     plot_bc,'settings_spectrum.fits.gz',0.5
;
; MODIFICATION HISTORY:
;     C. Evoli Original version May 2013
;-
;                       Check the parameters.
  On_error, 2
  compile_opt idl2
  
  np = N_params()
  IF (np NE 2) THEN BEGIN
     print, "Syntax: plot_bc, filename, phi"
     RETURN
  ENDIF
  
; plot canvas 
  xt='Log Kinetic Energy [GeV/n]'
  yt='B/C'

  IF Keyword_Set(PS) THEN BEGIN
     set_plot,'ps'
     DEVICE,/ENCAPSULATED,FILENAME='BC.eps',xsize=16,ysize=16
  ENDIF ELSE BEGIN
     window,0
  ENDELSE
  
  plot,[-1,3],[0,0.5],xtitle=xt,ytitle=yt,/NoData,xstyle=1,ystyle=1,charsize=2.,xthick=4.5,$
       ythick=4.5,charthick=3.2,thick=4.5
  
; extract model
  fits_info,filename,N_ext=N,/silent
  
  h = headfits(filename, exten=0)
  E = FXPAR(h,'Ekmin') * FXPAR(h,'Ekin_fac')^dindgen(FXPAR(h,'dimE'))
  B10 = dblarr(FXPAR(h,'dimE'))
  B11 = dblarr(FXPAR(h,'dimE'))
  C12 = dblarr(FXPAR(h,'dimE'))
  C13 = dblarr(FXPAR(h,'dimE'))
  C14 = dblarr(FXPAR(h,'dimE'))
  
  FOR i=1,N DO BEGIN
     h = headfits(filename, exten=i)
     IF FXPAR(h,'A') EQ 10 AND FXPAR(h,'Z_') EQ 5 then B10 = B10 + readfits(filename, exten=i, /silent)
     IF FXPAR(h,'A') EQ 11 AND FXPAR(h,'Z_') EQ 5 then B11 = B11 + readfits(filename, exten=i, /silent)
     IF FXPAR(h,'A') EQ 12 AND FXPAR(h,'Z_') EQ 6 then C12 = C12 + readfits(filename, exten=i, /silent)
     IF FXPAR(h,'A') EQ 13 AND FXPAR(h,'Z_') EQ 6 then C13 = C13 + readfits(filename, exten=i, /silent)
     IF FXPAR(h,'A') EQ 14 AND FXPAR(h,'Z_') EQ 6 then C14 = C14 + readfits(filename, exten=i, /silent)    
  endfor

  oplot,alog10(E),(B10+B11)/(C12+C13+C14),linestyle=2, thick=3;, color=cgColor("blue")

  print, (B10+B11)/(C12+C13+C14)
  
; modulate model
  IF phi GT 0. THEN BEGIN
     mp = .939d
     F = 5./10.*phi & B10_mod = ((E+mp)^2-mp^2)/((E+mp+F)^2-mp^2)*10.^CSPLINE(alog10(E),alog10(B10),alog10(E+F))
     F = 5./11.*phi & B11_mod = ((E+mp)^2-mp^2)/((E+mp+F)^2-mp^2)*10.^CSPLINE(alog10(E),alog10(B11),alog10(E+F))
     F = 6./12.*phi & C12_mod = ((E+mp)^2-mp^2)/((E+mp+F)^2-mp^2)*10.^CSPLINE(alog10(E),alog10(C12),alog10(E+F))
     F = 6./13.*phi & C13_mod = ((E+mp)^2-mp^2)/((E+mp+F)^2-mp^2)*10.^CSPLINE(alog10(E),alog10(C13),alog10(E+F))
     F = 6./14.*phi & C14_mod = ((E+mp)^2-mp^2)/((E+mp+F)^2-mp^2)*10.^CSPLINE(alog10(E),alog10(C14),alog10(E+F))

     oplot,alog10(E),(B10_mod+B11_mod)/(C12_mod+C13_mod+C14_mod),linestyle=0, thick=3;, color=cgColor("blue")
     ;oplot,alog10(E),(B10+B11)/(C12+C13+C14),linestyle=1, thick=3;, color=cgColor("blue")

     print,     (B10_mod+B11_mod)/(C12_mod+C13_mod+C14_mod)
  ENDIF

; plot PAMELA data from CRDB ( http://lpsc.in2p3.fr/cosmic-rays-db/ )
  data = READ_ASCII('HEAO3_BC.txt',data_start=2) & data = data.field01
  oploterror,alog10(data[0,*]),data[3,*],data[8,*],psym=3,/LOBAR
  oploterror,alog10(data[0,*]),data[3,*],data[9,*],psym=3,/HIBAR

  IF Keyword_Set(PS) THEN BEGIN
     DEVICE,/CLOSE
     set_plot,'X'
  ENDIF
  
  return 
END

