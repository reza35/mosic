
;+
; NAME:
;
; eis_width2velocity
;
; PURPOSE:
;
; Compute thermal and nonthermal velocities.
;
; CATEGORY:
;
; Data analysis.
;
; CALLING SEQUENCE:
;
; to compute nonthermal velocities:
; v=eis_width2velocity(element,ion,wavelength,width)
; v=eis_width2velocity(element,ion,wavelength,width,instr_fwhm=instr_fwhm)
;
; INPUTS:
;
; The element, ionization state, and wavelength must be passed. If a
; nonthermal velocity is desired the observed line FWHM (in Angstroms)
; must also be passed. The keyword instr_fwhm is used to specify the
; instrumental FWHM (in Angstroms). If only the thermal velocity and
; corresponding FWHM is desired the only the element, ionization
; state, and wavelength need to be passed. 
;
; Note that FWHM = 2.355*Standard Deviation of Gaussian Distribution.
;
; OPTIONAL INPUTS:
;
; instr_fwhm: The instrumental FWHM in Angstroms.
;       
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;
; The nonthermal velocity is the default output.
;
; OPTIONAL OUTPUTS:
;
; ti_max: The peak temperature in the ionization balance.
;
; COMMON BLOCKS:
;
; w2v_com: Stores a list of thermal velocities for each ion.
;
; SIDE EFFECTS:
;
; ?
;
; RESTRICTIONS:
;
; Velocities can be computed only for those ions covered by the
; ionization balance calculation (see the eis_width2velocity.dat
; file). Also note that for many singly ionized ions the temperature
; of formation is assumed to be 10,000 K and this may be an
; overstatement.
;
; PROCEDURE:
;
; FWHM^2 = (instr_fwhm)^2 + 4*ln(2)*(wavelength/c)^2(v_t^2+v_nt^2)
;
; EXAMPLE:
;
; Some example thermal widths:
;
; v_t = eis_width2velocity('C','II',1335,/thermal,ti_max=t,therm_fwhm=fwhm)
;
;   ion      wave A   temp K fwhm A
;  C   II   1335.00  2.2e+04  0.041
;  C   IV   1548.00  1.1e+05  0.105
;  N    V   1238.00  1.8e+05  0.101
;  O    I   1355.00  1.0e+04  0.024
;  O   IV   1401.00  1.7e+05  0.102
;  Si  II   1304.00  1.4e+04  0.021
;  Si III   1206.00  3.0e+04  0.028
;  Si  IV   1393.00  6.3e+04  0.047
;  S    I   1472.00  1.0e+04  0.019
;
; Compare with Feldman, Doschek, and Patterson, ApJ, 209:270-281, 1976
;
; An example nonthermal velocity:
;
; v_nt = eis_width2velocity('S','I',1900.29,0.04,ti_max=t)
; v_nt = 3.0 km/sec
; t    = 10,000 K
; 
; Compare with Mariska, Feldman, and Doschek, ApJ, 226:698-705, 1978
;
; MODIFICATION HISTORY:
;
;   HPW 01-SEP-1996: Orginal version.
;
;   HPW 23-JAN-1997: Modified to use only FWHM instead of standard
;                    deviation of Gaussian fit.
;
;   HPW 15-MAY-2007: Renamed eis_width2velocity.
;
;-

function eis_width2velocity,input_element,input_ion,wavelength,obs_fwhm,$
                            instr_fwhm=instr_fwhm,$
                            therm_fwhm=therm_fwhm,$
                            nt_fwhm=nt_fwhm,$
                            thermal=thermal,ti_max=ti_max,$
                            ti_use=ti_use

element = strupcase(input_element)
ion     = strupcase(input_ion)

cc = 2.9979E+5 ;; speed of light in km/sec

;; ----------------------------------------------------------------------
;; READ LIST OF THERMAL VELOCITES BASED 
;; ----------------------------------------------------------------------

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

;; ----------------------------------------------------------------------
;; COMPUTE WIDTHS AND VELOCITIES
;; ----------------------------------------------------------------------

if (count gt 0) then begin

  ;; thermal velocity and temperature of formation
  if not(keyword_set(ti_use)) then begin
    ti_max = t_max(match(0))
    v_t    = v_therm(match(0))
  endif else begin
    mass   = eis_element2mass(element)
    v_t    = sqrt(2*(10.0^ti_use)*1.6022E-12/11604.0/mass)/1.0E+5
    ti_max = ti_use
    print,v_t
  endelse

  ;; instrumental FWHM in Angstroms
  if keyword_set(instr_fwhm) then dl_i = instr_fwhm else dl_i = 0

  ;; thermal FWHM in Angstroms
  dl_t = sqrt(4*alog(2))*wavelength*v_t/cc

  ;; return thermal velocity not nonthermal velocity
  therm_fwhm = dl_t
  if keyword_set(thermal) then return,v_t

  ;; observed FWHM in Angstroms
  dl_o = obs_fwhm

  ;; nonthermal FWHM in Angstroms
  dl_nt_2 = (dl_o^2 - dl_i^2 - dl_t^2)

  ;; return nonthermal velocity
  if dl_nt_2 gt 0 then begin
    v_nt    = cc*(sqrt(dl_nt_2)/wavelength)/sqrt(4*alog(2))
    nt_fwhm = sqrt(dl_nt_2)
  endif else begin
    v_nt    = -1.0
    nt_fwhm = -1.0
    message,"Observed width smaller than thermal width",/informational
  endelse
  return,v_nt

endif else message," Ion not found: "+element+" "+ion,/informational

return,-1.0
end
