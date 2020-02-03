PRO plot_leptons, filename, phi, PS=ps
;+
; NAME:
;     PLOT_LEPTONS
; PURPOSE:
;     Plot lepton spectra from DRAGON output against PAMELA data
;
; CALLING SEQUENCE:
;     plot_lepton, filename, phi [, /PS]
;
; INPUTS:
;     filename = The input file
;     phi = modulation potential [GV]
;
; OPTIONAL INPUT:
;     /PS = The plot is saved in lepton.eps file
;
; EXAMPLES:
;     plot_leptons,'example_spectrum.fits.gz',0.5
;
; MODIFICATION HISTORY:
;     C. Evoli Original version May 2013
;-
; Check the parameters.
  On_error, 2
  compile_opt idl2
  
  np = N_params()
  IF (np NE 2) THEN BEGIN
     print, "Syntax: plot_leptons, filename, phi"
     RETURN
  ENDIF
  
; plot canvas 
  xt='Log Kinetic Energy [GeV]'
  yt='e!E+!N/(e!E+!N+e!E-!N)'
  
  IF Keyword_Set(PS) THEN BEGIN
     set_plot,'ps'
     DEVICE,/ENCAPSULATED,FILENAME='leptons.eps',xsize=16,ysize=16
  ENDIF ELSE BEGIN
     window,0
  ENDELSE
  
  plot,[-1,3],[1d-2,1d0],xtitle=xt,ytitle=yt,/NoData,charsize=2,/ylog
  
; extract model
  fits_info,filename,N_ext=N,/silent
  
  h = headfits(filename, exten=0)
  E = FXPAR(h,'Ekmin') * FXPAR(h,'Ekin_fac')^dindgen(FXPAR(h,'dimE'))
  J_ele = dblarr(FXPAR(h,'dimE')) & J_pos = J_ele & J_extra = J_ele
  
  FOR i=1,N DO BEGIN
     h = headfits(filename, exten=i)
     IF FXPAR(h,'A') EQ 0 AND FXPAR(h,'Z_') EQ  1 AND FXPAR(h,'EXTRA') EQ 0 then J_pos = J_pos + readfits(filename, exten=i, /silent)
     IF FXPAR(h,'A') EQ 0 AND FXPAR(h,'Z_') EQ -1 AND FXPAR(h,'EXTRA') EQ 0 then J_ele = J_ele + readfits(filename, exten=i, /silent)
  
     IF FXPAR(h,'extra') EQ 1 THEN J_extra = J_extra + readfits(filename, exten=i, /silent)
endfor
  
  J_ele = J_ele + J_extra
  J_pos = J_pos + J_extra

  oplot,alog10(E),J_pos/(J_pos+J_ele)

; modulate model
  IF phi GT 0. THEN BEGIN
     mp = .511d-3
     
     F_ele = 1.d0*phi & J_ele_mod = ((E+mp)^2-mp^2)/((E+mp+F_ele)^2-mp^2)*10.^CSPLINE(alog10(E),alog10(J_ele),alog10(E+F_ele))
     F_pos = 1.d0*phi & J_pos_mod = ((E+mp)^2-mp^2)/((E+mp+F_pos)^2-mp^2)*10.^CSPLINE(alog10(E),alog10(J_pos),alog10(E+F_pos))
     
     oplot,alog10(E),J_pos_mod/(J_pos_mod+J_ele_mod),linestyle=2
  ENDIF

; plot PAMELA data from CRDB ( http://lpsc.in2p3.fr/cosmic-rays-db/ )
  data = READ_ASCII('PAMELA_posratio.txt',data_start=2) & data = data.field01
  oploterror,alog10(data[0,*]),data[3,*],data[8,*],psym=3,/LOBAR
  oploterror,alog10(data[0,*]),data[3,*],data[9,*],psym=3,/HIBAR

  IF Keyword_Set(PS) THEN BEGIN
     DEVICE,/CLOSE
     set_plot,'X'
  ENDIF
  
  return 
END
