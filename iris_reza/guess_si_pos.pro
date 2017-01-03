pro guess_si_pos, inpath, inp, dispersion, SiIV_range1, SiIV_range2

;+
;===============================================================  
; procedure : guess_si_pos.pro
;  
; purpose : to find line profiles in the Si IV spectra and create a wavelength scale
;  
; 2015/12/30 :  created
; 2016/04/19 :  disp=5 pm was added.
;  
; R.Rezaei @ IAC                         e-mail:  rrezaei@iac.es      
;===============================================================
;

  
on_error, 2

speed_of_light = 299792.458     ; km/s in vacuum

dim = get_screen_size(RESOLUTION=resolution)

u = max(inp[SiIV_range1[0]:SiIV_range1[1]], pos)
master_si1 = pos
u = max(inp[SiIV_range2[0]:SiIV_range2[1]], pos)
master_si2 = pos

n1 = SiIV_range1[1] - SiIV_range1[0] + 1.
n2 = SiIV_range2[1] - SiIV_range2[0] + 1.

if (dispersion lt 1.3) then si_w_accept = 50.   ; maximum acceptable line width
if (dispersion gt 2.5) then si_w_accept = 25.

del_si = ((140300.0d * 50.0 / speed_of_light)/ dispersion) ; = max delta lambda in px  equivalent to 50 km/s

inp -= (min(inp[SiIV_range2[0]:SiIV_range2[1]]) -1)

save, filename=inpath+'si_lines_position.sav', master_si1, master_si2, n1, n2, del_si, inp, si_w_accept

end  
