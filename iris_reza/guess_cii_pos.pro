pro guess_cii_pos, inpath, inp, dispersion

;+
;===============================================================  
; procedure : guess_cii_pos.pro
;  
; purpose : to find line profiles in the C II spectra and create a wavelength scale
;  
; 2015/12/30 :  created
; 2016/02/18 :  disp=5 pm was added.
;  
; R.Rezaei @ IAC                         e-mail:  rrezaei@iac.es      
;===============================================================
;

  
on_error, 2

speed_of_light = 299792.458     ; km/s in vacuum

dim = get_screen_size(RESOLUTION=resolution)

ca_ilo = fltarr(2)
ca_ihi = fltarr(2)

if (dispersion gt 2.4) then begin
   shit = max(inp, pmax)
   line2 = pmax
   lpff, inp[line2-10:line2+10], posmax
   line2 = line2-10. + posmax
   
   shit = max(inp[0:line2 -20], pmax)
   line1 = pmax
   lpff, inp[line1-10:line1+10], posmax
   line1 = line1-10. + posmax
endif else begin
   shit = max(inp, pmax)
   line2 = pmax
   lpff, inp[line2-20:line2+20], posmax
   line2 = line2-20. + posmax
   
   shit = max(inp[0:line2 -40], pmax)
   line1 = pmax
   lpff, inp[line1-20:line1+20], posmax
   line1 = line1-20. + posmax
endelse   

if (dispersion gt 4.) then begin
   shit = max(inp, pmax)
   line2 = pmax
   lpff, inp[line2-5:line2+5], posmax
   line2 = line2 - 5. + posmax
   
   shit = max(inp[0:line2 -10], pmax)
   line1 = pmax
   lpff, inp[line1-5:line1+5], posmax
   line1 = line1-5. + posmax
endif


if 0 then begin
satisfied = 0.

print, 'if marked positions are NOT ok, press N or n'
antw = get_kbrd()
if (antw eq 'N')or(antw eq 'n') then begin
   print
   print, 'first click on the peak and then a red-wing position acceptable as maximum velocity: first line'
   read_plot_click, 2, vv, vv, pos_si
   master_c1 = pos_si[0]
   lpff, vv[master_c1-9: master_c1+9], ppos
   master_c1 = ppos + master_c1-9
   del_c1 = ((133500.0d * 20.0 / speed_of_light)/ dispersion) < abs(pos_si[0] - pos_si[1]) 
   ca_ilo[0] = (master_c1 - 130.) > 0.
   ca_ihi[0] = (master_c1 + 60.) < (n_elements(inp) - 20.)
   print

   print, 'first click on the peak and then a red-wing position acceptable as maximum velocity: second line'
   read_plot_click, 2, vv, vv, pos_si
   master_c2 = pos_si[0]
   lpff, vv[master_c2-9: master_c2+9], ppos
   master_c2 = ppos + master_c2-9
   del_c2 =  ((133600.0d * 20.0 / speed_of_light)/ dispersion) < abs(pos_si[0] - pos_si[1]) 
   ca_ilo[1] = (master_c2 - 50.) > 0.
   ca_ihi[1] = (master_c2 + 130.) < (n_elements(inp) - 20.)
endif
endif

master_c1 = line1
del_c1 = ((133500.0d * 20.0 / speed_of_light)/ dispersion)
ca_ilo[0] = (master_c1 - 130.) > 0.
ca_ihi[0] = (master_c1 + 60.) < (line2 - 20.)

master_c2 = line2
del_c2 = ((133500.0d * 20.0 / speed_of_light)/ dispersion)
ca_ilo[1] = (master_c2 - 130.) >  ca_ihi[0]
ca_ihi[1] = (master_c2 + 60.) < (n_elements(inp) - 20.)

ca_ilo = round(ca_ilo)
ca_ihi = round(ca_ihi)

save, filename=inpath+'cii_lines_position.sav', master_c1, master_c2, del_c1, del_c2, ca_ilo, ca_ihi

end  
