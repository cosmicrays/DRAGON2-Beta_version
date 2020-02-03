PRO plot_protons, filename, phi, PS=ps
;+
; NAME:
;     PLOT_PROTONS
; PURPOSE:
;     Plot Proton spectrum from DRAGON output against PAMELA data
;
; CALLING SEQUENCE:
;     plot_protons, filename, phi [, /PS]
;
; INPUTS:
;     filename = The input file
;     phi = modulation potential [GV]
;
; OPTIONAL INPUT:
;     /PS = The plot is saved in protons.eps file
;
; EXAMPLES:
;     plot_protons,'example_spectrum.fits.gz',0.5
;
; MODIFICATION HISTORY:
;     C. Evoli Original version May 2013
;-
; Check the parameters.
  On_error, 2
  compile_opt idl2
  
  np = N_params()
  IF (np NE 2) THEN BEGIN
     print, "Syntax: plot_protons, filename, phi"
     RETURN
  ENDIF
  
; plot canvas 
  xt='Log Kinetic Energy [GeV]'
  yt='E!E2!NJ!Dp!N [GeV m!E-2!N s!E-1!N sr!E-1!N]'
  
  IF Keyword_Set(PS) THEN BEGIN
     set_plot,'ps'
     DEVICE,/ENCAPSULATED,FILENAME='protons.eps',xsize=16,ysize=16
  ENDIF ELSE BEGIN
     window,0
  ENDELSE
  
  plot,[-1,3],[1d0,1d4],xtitle=xt,ytitle=yt,/NoData,charsize=2,/ylog
  
; extract model
  fits_info,filename,N_ext=N,/silent
  
  h = headfits(filename, exten=0)
  E = FXPAR(h,'Ekmin') * FXPAR(h,'Ekin_fac')^dindgen(FXPAR(h,'dimE'))
  J = dblarr(FXPAR(h,'dimE'))
  
  FOR i=1,N DO BEGIN
     h = headfits(filename, exten=i)
     IF FXPAR(h,'A') EQ 1 AND FXPAR(h,'Z_') EQ 1 then J = J + readfits(filename, exten=i, /silent)
  endfor
  
  oplot,alog10(E),E^2.0*J
  
; modulate model
  IF phi GT 0. THEN BEGIN
     mp = .939d
     F = 1.d0*phi & J_mod = ((E+mp)^2-mp^2)/((E+mp+F)^2-mp^2)*10.^CSPLINE(alog10(E),alog10(J),alog10(E+F))
     oplot,alog10(E),E^2.0*J_mod,linestyle=2
  ENDIF

; plot PAMELA data from CRDB ( http://lpsc.in2p3.fr/cosmic-rays-db/ )
  data = READ_ASCII('PAMELA_protons.txt',data_start=2) & data = data.field01
  oploterror,alog10(data[0,*]),data[3,*],data[8,*],psym=3,/LOBAR
  oploterror,alog10(data[0,*]),data[3,*],data[9,*],psym=3,/HIBAR

  IF Keyword_Set(PS) THEN BEGIN
     DEVICE,/CLOSE
     set_plot,'X'
  ENDIF
  
  return 
END
