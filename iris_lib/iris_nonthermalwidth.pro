; NAME:
;
; iris_nonthermalwidth
;
; PURPOSE:
;
; Compute thermal and nonthermal widths in the unit of km/s
;
; CATEGORY:
;
; Data analysis.
;
; CALLING SEQUENCE:
;
; Wnt_v=iris_nonthermalwidth(element,ion,wavelength,Wobs_v,instr_fwhm,ti_use=ti_use,ti_max=ti_max,Wt_v=Wt_v)
;
; INPUTS:
;
; element - element of the line, e.g., 'Si'
; 
; ion - ionization state, e.g., 'IV'
;
; wavelength - wavelength of the line
; 
; Wobs_v - observed (1/e) line width in the unit of km/s, Wobs_v can be an array of any dimension
; 
; instr_fwhm - the instrumental FWHM (in Angstroms)
;       
; KEYWORD PARAMETERS:
; 
; ti_use - the ion temperature, if don't want to use the default one in Harry Warren's file of eis_width2velocity.dat
;
; OUTPUTS:
;
; The nonthermal velocity is the default output.
;
; OPTIONAL OUTPUTS:
;
; ti_max - The peak temperature in the ionization balance.
; 
; Wt_v - thermal width in the unit of km/s
;
; COMMON BLOCKS:
;
; w2v_com: Stores a list of thermal velocities for each ion.
;
; RESTRICTIONS:
;
; Velocities can be computed only for those ions covered by the
; ionization balance calculation (see Harry Warren's file of eis_width2velocity.dat). 
; Also note that for many singly ionized ions the temperature
; of formation is assumed to be 10,000 K and this may be an
; overstatement.
;
; HISTORY: Written by Hui Tian at CfA, April 8, 2013. The routinely was modified fro IRIS from Harry Warren's EIS routine eis_width2velocity.pro
;

function iris_nonthermalwidth,element,ion,wavelength,Wobs_v,instr_fwhm,$
      ti_use=ti_use,ti_max=ti_max,Wt_v=Wt_v

element = strupcase(element)
ion     = strupcase(ion)

cc = 2.9979E+5 ;; speed of light in km/sec

;;read Dr. Harry Warren's file of temperature and thermal velocity eis_width2velocity.dat

common w2v_com,init,data,elements,ions,t_max,v_therm

if n_elements(init) eq 0 then init = 1

if init then begin
  file = find_with_def('eis_width2velocity.pro',!path)
  file = str_replace(file,'.pro','.dat')
  data     = rd_tfile(file,/auto)
  elements = strupcase(reform(data(0,*)))
  ions     = strupcase(reform(data(1,*)))
  t_max    = float(reform(data(2,*)))
  v_therm  = float(reform(data(3,*)))
  init     = 0
endif

match = where(element eq elements and ion eq ions,count)


if (count gt 0) then begin

  ;; thermal (1/e) width in the unit of km/s
  if not(keyword_set(ti_use)) then begin
    ti_max = t_max(match(0))
    Wt_v    = v_therm(match(0))
  endif else begin
    mass   = eis_element2mass(element,/kg)
    Wt_v    = sqrt(2*(10.0^ti_use)*1.3806488E-23/mass)/1.0E+3
    ti_max = ti_use
  endelse


  ;; instrumental (1/e) width in the unit of km/s
  instr_v=instr_fwhm/sqrt(alog(2.))/2./wavelength*cc

  ;; nonthermal FWHM in Angstroms
  Wnt_v_sq = (Wobs_v^2 - instr_v^2 - Wt_v^2)
  Wnt_v = Wobs_v
  sub_good=where(Wnt_v_sq ge 0)
  if sub_good[0] ne -1 then Wnt_v[sub_good] = sqrt(Wnt_v_sq[sub_good]) else message,"Observed width smaller than thermal width",/informational
  sub_bad=where(Wnt_v_sq lt 0)
  if sub_bad[0] ne -1 then Wnt_v[sub_bad] = -1
  
    
  return,Wnt_v

endif else message," Ion not found: "+element+" "+ion,/informational


end
