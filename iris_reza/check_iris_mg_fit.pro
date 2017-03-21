pro check_iris_mg_fit, inpath

;+
;===============================================================
; procedure : check_iris_mg_fit.pro
;  
; purpose : to check multi-Gaussian fit to Mg II h/k lines
;
; inpath :  the directory containing the save file 
;
; May 26, 2016 : created
;  
; R.Rezaei @ IAC                         e-mail:  rrezaei@iac.es      
;===============================================================
;-

device, decomposed=0
  
restore, inpath + 'iris_mg_analyzed_scan_0.sav';, /v
restore, inpath + 'iris_cont_analyzed_scan_0.sav';, /v

loadct,0,/silent
window, 14, title='reduced chi-square statistics, Mg II fits'
s = 12.
x= reform(k_spar[*, *, 4]) & q=where(x ne 0.) & plot_histogram, x(q), xrange=[0, s]
loadct, 40,/silent
x= reform(k_dpar[*, *, 7]) & q=where(x ne 0.) & plot_histogram, x(q), xrange=[0, s], oplott=70
x= reform(k_tpar[*, *, 8]) & q=where(x ne 0.) & plot_histogram, x(q), xrange=[0, s], oplott=250
loadct, 0,/silent
 
k_fit_gauss = k_fit_gauss[*,*,0:47,*]
tay=size(k_fit_gauss)
window, 0, xs=tay[1]*3, ys=tay[2], tit='Continuum, Amplitude'
window, 2, xs=700, ys=550, tit='profile'

vsig = k2v_b * mg_wing  &  vsig = shift(vsig, [0,1])
osig = k3_b * mg_wing  &  osig = shift(osig, [0,2])
tvscl, mg_wing

finished = 0
while (finished eq 0) do begin
  cursor,x , y, /device, /change

; calculate position 
px = x mod (tay[1]) & py = y mod (tay[2])

loadct, 0,/silent
wset, 0
tvscl, reform(em_lines[*,*,1]) ;cont ;mg_wing
tvscl, vsig, tay[1], 0
tvscl, osig, tay[1]*2, 0
;-------------------------------
;---   draw cursor position
;-------------------------------
loadct,40,/silent
for i=0,2 do plots,[px,px]+tay[1]*i,[0,tay[2]],/device, color=250
for i=0,0 do plots,[0,tay[1]*3],[py,py]+tay[2]*i,/device, color=250
loadct,0,/silent

;x = px 
;y = py

xyouts, 2, 2, num2string(fix(px))+','+num2string(fix(py)), /device, color=0,charsize=2

; quit: right mouse button
if !mouse.button eq 4 then finished = 1

;----------------------------------
;--  left button: inversion result
;----------------------------------
starter = 0
if !mouse.button eq 1  then starter = 1
if !mouse.button eq 2  then starter = 2
wait, .1
x2 = (findgen(50) - 25.)* dispersion_nuv

if (starter eq 1) then begin
     !p.font=0
     print, px, py, k_spar[px,py,4], k_dpar[px,py,7], k_tpar[px,py,8], k_tpar[px, py,2], k_tpar[px, py,5]
     wset, 2
     plot, k_fit_gauss[px, py, *, 0],/xsty;,/ylog, yr=[0.1, 1d4]
     ee = 0.29 * sqrt(k_fit_gauss[px, py, *, 0])
     errplot, k_fit_gauss[px, py, *, 0] - ee, k_fit_gauss[px, py, *, 0] + ee,color=220
     loadct,40,/silent
     oplot, k_fit_gauss[px, py, *, 2], color=175,thick=2.5, linestyle=2
     oplot, k_fit_gauss[px, py, *, 3], color=245,thick=2.5
     loadct,0,/silent
     print, '------------------------------------------------------------------------'
     wset, 0
  endif

  if (starter eq 2) then begin
     !p.font=-1
     set_plot, 'ps'
     device, filename='~/iris_mg_prof_'+strtrim(px,1)+'_'+strtrim(py,1)+'.eps', /encapsulated, /color
     x = (findgen(tay[3]) - tay[3]/2.)* dispersion_nuv
     print, px, py, k_spar[px,py,4], k_dpar[px,py,7], k_tpar[px,py,8], k_tpar[px, py,2], k_tpar[px, py,5]
     plot, x, k_fit_gauss[px, py, *, 0],/xsty,xthick=2,ythick=2,thick=2,chars=1.3,xtit='!7Dk!3 [pm]', $
           ytit='intensity [DN]',charthick=2, yr=[0, max(k_fit_gauss[px,py,*, 0])*1.1]
     ee = ir_error(reform(k_fit_gauss[px, py, *, 0]), /nuv, /dark)
     errplot, x, k_fit_gauss[px, py, *, 0] - ee, k_fit_gauss[px, py, *, 0] + ee, thick=2
     loadct,40,/silent  
     oplot, x, k_fit_gauss[px, py, *, 3], color=245,thick=2.5
     loadct,0,/silent
     set_plot, 'x'
     !p.font=0
  endif


  
endwhile


stop

end
