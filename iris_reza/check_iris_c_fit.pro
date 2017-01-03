pro check_iris_c_fit, inpath

;+
;===============================================================
; procedure : check_iris_c_fit.pro
;  
; purpose : to check multi-Gaussian fit to the C II lines
;
; inpath :  the directory containing the save file 
;
; May 26, 2016 : created
;  
; R.Rezaei @ IAC                         e-mail:  rrezaei@iac.es      
;===============================================================
;-

device, decomposed=0
  
restore, inpath + 'iris_cii_analyzed_scan_0.sav', /v

loadct,0,/silent
window, 14, title='reduced chi-square statistics, C II fits'
s = 9
x= reform(cii_s[*, *, 4])  & q=where(x ne 0.)  &  plot_histogram, x(q), xrange=[0, s]
loadct,40,/silent
x= reform(cii_d[*, *, 5])  & q=where(x ne 0.)  &  plot_histogram, x(q), xrange=[0, s], oplott=70
x= reform(cii_q[*, *, 10]) & q=where(x ne 0.)  &  plot_histogram, x(q), xrange=[0, s], oplott=250
loadct,0,/silent

tay=size(cii_fit_gauss)
window, 0, xs=tay[1]*3, ys=tay[2], tit='Continuum, Amplitude'
window, 2, xs=700, ys=550, tit='profile'

vsig = alog10(reform(cii_s[*,*,2])>5.) 
osig = reform(cii_s[*,*,3])
tvscl, reform(cii_s[*,*,0]) >.1 

finished = 0
while (finished eq 0) do begin
  cursor,x , y, /device, /change

; calculate position 
px = x mod (tay[1]) & py = y mod (tay[2])

loadct, 0,/silent
wset, 0
tvscl, reform(cii_s[*,*,0]) < 10. > 0.
tvscl, vsig, tay[1], 0
tvscl, osig, tay[1]*2, 0
;-------------------------------
;---   draw cursor position
;-------------------------------
loadct,40,/silent
for i=0,2 do plots,[px,px]+tay[1]*i,[0,tay[2]],/device, color=250
for i=0,0 do plots,[0,tay[1]*3],[py,py]+tay[2]*i,/device, color=250
loadct,0,/silent

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

  if (starter eq 1) then begin
     print, px, py, cii_s[px,py,4], cii_d[px,py,5], cii_q[px,py,10], cii_p[px,py,13]
     wset, 2
     plot, cii_fit_gauss[px, py, *, 0],/xsty
     ee = 0.5 * sqrt(cii_fit_gauss[px, py, *, 0])
     errplot, cii_fit_gauss[px, py, *, 0] - ee, cii_fit_gauss[px, py, *, 0] + ee,color=220
     loadct,40,/silent
     oplot, cii_fit_gauss[px, py, *, 1], color=175,thick=2.5, linestyle=2
     oplot, cii_fit_gauss[px, py, *, 3], color=245,thick=2.5
     loadct,0,/silent
     print, '------------------------------------------------------------------------'
     wset, 0
  endif
  if (starter eq 2) then begin
     !p.font=-1
     set_plot, 'ps'
     device, filename='~/iris_cii_prof_'+strtrim(px,1)+'_'+strtrim(py,1)+'.eps', /encapsulated, /color
     x = (findgen(tay[3]) - master_c1)* dispersion_fuv * 1.0e-3 + 133.4532

     print, px, py, cii_s[px,py,4], cii_d[px,py,5], cii_q[px,py,10], cii_p[px,py,13]
     plot, x, cii_fit_gauss[px, py, *, 0],/xsty,xthick=2,ythick=2,thick=2,chars=1.3,xtit='!7k!3 [pm]', $
           ytit='intensity [DN]',charthick=2;,/ylog, yr=[1.0, max(si2_fit_gauss[px, py, *, 0])*1.2]
     ee = ir_error(reform(cii_fit_gauss[px, py, *, 0]), /fuv, /dark)
     errplot, x, cii_fit_gauss[px, py, *, 0] - ee, cii_fit_gauss[px, py, *, 0] + ee, thick=2
     loadct,40,/silent  
     oplot, x, cii_fit_gauss[px, py, *, 2], color=245,thick=5.5
     loadct,0,/silent
     set_plot, 'x'
     !p.font=0
  endif
  
endwhile


stop

end
