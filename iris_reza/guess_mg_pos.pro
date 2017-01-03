pro guess_mg_pos, inpath, inp, dispersion, do_h_line=do_h_line

;+
;===============================================================  
; procedure : guess_mg_pos.pro
;  
; purpose : to find line profiles in the Mg II k/h spectra and create a
; wavelength scale
;  
; 2015/12/30 :  created
; 2016/01/14 :  disp=5 pm was added.
;  
; R.Rezaei @ IAC                         e-mail:  rrezaei@iac.es      
;===============================================================
;

  
on_error, 2

speed_of_light = 299792.458     ; km/s in vacuum

if (n_elements(do_h_line) eq 0) then do_h_line = 0 else do_h_line = 1
dim = get_screen_size(RESOLUTION=resolution)

if (dispersion lt 3.) then begin
     out = [104, 111, 117, 65, 94, 125, 148, 229, 237, 324, 334]  ; dummy variable
     if (do_h_line eq 1) then out = [out, 386, 392, 398, 342, 378, 406, 424]
     window, xs=dim[0], ys=600
     plot, inp,/xst
     loadct, 40,/silent
     for i=0, 2   do oplot, [out[i],out[i]],[0, 9000], color=245, linestyle=1
     for i=3, 6   do oplot, [out[i],out[i]],[0, 2000], color=75, linestyle=2
     for i=7, 10  do oplot, [out[i],out[i]],[0, 2000], color=215
     for i=11, 13  do oplot, [out[i],out[i]],[0, 9000], color=245, linestyle=1
     for i=14, 17 do oplot, [out[i],out[i]],[0, 2000], color=75, linestyle=2
     loadct,0,/silent
     mg_xaxis = (findgen(n_elements(inp))- out[1])*dispersion*1.0d-3 + 279.5528d  ; x-axis of the Mg channel in nm
  endif else begin
     ;print, 'TBD.......'
     out = [107, 111, 114, 85, 102, 119, 130, 161, 167, 215, 223]
     if (do_h_line eq 1) then out = [out, 248, 251, 254, 227, 243, 259, 274]
     window, xs=dim[0], ys=600
     plot, inp,/xst
     loadct, 40,/silent
     for i=0, 2   do oplot, [out[i],out[i]],[0, 9000], color=245, linestyle=1
     for i=3, 6   do oplot, [out[i],out[i]],[0, 2000], color=75, linestyle=2
     for i=7, 10  do oplot, [out[i],out[i]],[0, 2000], color=215
     for i=11, 13  do oplot, [out[i],out[i]],[0, 9000], color=245, linestyle=1
     for i=14, 17 do oplot, [out[i],out[i]],[0, 2000], color=75, linestyle=2
     loadct,0,/silent
     mg_xaxis = (findgen(n_elements(inp))- out[1])*dispersion*1.0d-3 + 279.5528d  ; x-axis of the Mg channel in nm
  endelse   

print, 'if marked positions are NOT ok, press N or n'
antw = get_kbrd()
if (antw ne 'N')and(antw ne 'n') then begin
   master_k2v = out[0]
   master_k3 = out[1]
   master_k2r = out[2]
   master_k1v1 = out[3]
   master_k1v2 = out[4]
   master_k1r1 = out[5]
   master_k1r2 = out[6]
   master_k1v = (master_k1v1 + master_k1v2)/2.
   master_k1r = (master_k1r1 + master_k1r2)/2.

   ca_ilo = intarr(2)
   ca_ihi = intarr(2)
   ca_ilo[0] = out[7]   &  ca_ihi[0] = out[8]
   ca_ilo[1] = out[9]   &  ca_ihi[1] = out[10]

   if (do_h_line eq 1) then begin
      master_h2v = out[11]
      master_h3 = out[12]
      master_h2r = out[13]
      master_h1v1 = out[14]
      master_h1v2 = out[15]
      master_h1r1 = out[16]
      master_h1r2 = out[17]
      master_h1v = (master_h1v1 + master_h1v2)/2.
      master_h1r = (master_h1r1 + master_h1r2)/2.
   endif
endif else begin

vv = inp   
print, 'click Carefully on the Violet peak, K3, and Red peak, respectively.'
read_plot_click, 3, vv, vv, out
master_k2v = out[0]
master_k3 = out[1]
master_k2r = out[2]
print
print, 'click Carefully on two points for range of the Violet and Red minima of K line, respectively.'
read_plot_click, 4, vv, vv, out
master_k1v1 = out[0]
master_k1v2 = out[1]
master_k1r1 = out[2]
master_k1r2 = out[3]
master_k1v = (master_k1v1 + master_k1v2)/2.
master_k1r = (master_k1r1 + master_k1r2)/2.

if (do_h_line eq 1) then begin
  print, 'click Carefully on the Violet peak, H3, and Red peak, respectively.'
  read_plot_click, 3, vv, vv, out
  master_h2v = out[0]
  master_h3 = out[1]
  master_h2r = out[2]

  print, 'click Carefully on two points for range of the Violet and Red minima of H line, respectively.'
  read_plot_click, 4, vv, vv, out
  master_h1v1 = out[0]
  master_h1v2 = out[1]
  master_h1r1 = out[2]
  master_h1r2 = out[3]
  master_h1v = (master_h1v1 + master_h1v2)/2.
  master_h1r = (master_h1r1 + master_h1r2)/2.
  print
endif

print, 'click Carefully on the left/right borders of two photospheric blend.'
read_plot_click, 4, vv, vv, out
ca_ilo = intarr(2)
ca_ihi = intarr(2)
ca_ilo[0] = out[0]   &  ca_ihi[0] = out[1]
ca_ilo[1] = out[2]   &  ca_ihi[1] = out[3]
endelse

del_mg = ((279500.0d * 20.0 / speed_of_light)/ dispersion) ; velocities within 20 km/s

if (do_h_line eq 1) then begin
save, filename=inpath+'mg_lines_position.sav', del_mg, ca_ilo, ca_ihi, inp, dispersion, mg_xaxis, $
       master_k1r, master_k1r2, master_k1r1, master_k1v2, master_k1v1, master_k2r, master_k3, master_k2v, $
       master_h1r, master_h1r2, master_h1r1, master_h1v2, master_h1v1, master_h2r, master_h3, master_h2v
endif else begin
save, filename=inpath+'mg_lines_position.sav', del_mg, ca_ilo, ca_ihi, inp, dispersion, mg_xaxis, $
       master_k1r, master_k1r2, master_k1r1, master_k1v2, master_k1v1, master_k2r, master_k3, master_k2v
endelse

end  
