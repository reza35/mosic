pro check_iris_si_fit, inpath

;+
;===============================================================
; procedure : check_iris_si_fit.pro
;  
; purpose : to check multi-line fit to the Si IV 140.277 + O IV 1401 + ...
;
; inpath :  the directory containing the save file  
;
; Dec 15, 2014 : created
;  
; R.Rezaei @ IAC                         e-mail:  rrezaei@iac.es      
;===============================================================
;

device, decomposed=0
  
restore, inpath + 'iris_si_analyzed_scan_0.sav', /v

window, 14, title='reduced chi-square statistics, Si IV + O IV fits'

loadct,0,/silent
s = 5
x= reform(si_iv2[*, *, 4]) & q=where(x ne 0.) & plot_histogram, x(q), xrange=[0, s]
loadct,40,/silent
x= reform(si_iv2_v[*, *, 7]) & q=where(x ne 0.) & plot_histogram, x(q), xrange=[0, s],  oplott=70
x= reform(si_iv2_q[*, *, 7]) & q=where(x ne 0.) & plot_histogram, x(q), xrange=[0, s],  oplott=170
x= reform(si_iv2_z[*, *, 10]) & q=where(x ne 0.) & plot_histogram, x(q), xrange=[0, s], oplott=200
x= reform(si_iv2_f[*, *, 13]) & q=where(x ne 0.) & plot_histogram, x(q), xrange=[0, s], oplott=250
loadct,0,/silent


tay = size(si2_fit_gauss)
window, 0, xs=tay[1]*3, ys=tay[2], tit='Continuum, Si IV amplitude, chi-sqaure'
window, 2, xs=700, ys=550, tit='profile'

vsig = alog10(reform(si_iv2_q[*,*,1])>1.)<3.  &  vsig = shift(vsig, [0, 5])
osig = reform(si_iv2_q[*,*,7]) < 6.   &   osig = shift(osig, [0, 5])

tvscl, cont

finished = 0
while (finished eq 0) do begin
  cursor,x , y, /device, /change

; calculate position 
px = x mod (tay[1]) & py = y mod (tay[2])

loadct, 0,/silent
wset, 0
tvscl, cont
tvscl, vsig, tay[1], 0
tvscl, osig, tay[1]*2, 0
;-------------------------------
;---   draw cursor position
;-------------------------------
loadct,40,/silent
for i=0,2 do plots,[px,px]+tay[1]*i,[0,tay[2]],/device, color=250
for i=0,0 do plots,[0,tay[1]*3],[py,py]+tay[2]*i,/device, color=250
loadct,0,/silent

px = px 
py = py - 4

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
xx = (findgen(tay[3]))* 1.296e-3 + 139.90 - 0.02 ;dispersion_fuv

  if (starter eq 1) then begin
     print, px, py, si_iv2[px,py,4], si_iv2_v[px,py,7], si_iv2_q[px,py,7], si_iv2_z[px,py,12]
     wset, 2
     plot, xx, si2_fit_gauss[px, py, *, 0],/xsty;,/ylog, yr=[-1., 1d4]
     ee =0.5 * sqrt(si2_fit_gauss[px, py, *, 0])
     errplot, xx, si2_fit_gauss[px, py, *, 0] - ee, si2_fit_gauss[px, py, *, 0] + ee,color=220
     loadct,40,/silent
     oplot, xx, si2_fit_gauss[px, py, *, 3], color=175,thick=2.5, linestyle=2
     oplot, xx, si2_fit_gauss[px, py, *, 5], color=245,thick=2.5
     loadct,0,/silent
     print, '---------------------------------'
     wset, 0
  endif

  if (starter eq 2) then begin
     device, decomposed=0
     !p.font=-1
     set_plot, 'ps'
     s=3.
     device, filename='~/iris_si_prof_'+strtrim(px,1)+'_'+strtrim(py,1)+'.eps', /encapsulated, /color
     
     print, px, py, si_iv2[px,py,4], si_iv2_q[px,py,7], si_iv2_z[px,py,10], si_iv2_f[px,py,13]
     plot, xx, si2_fit_gauss[px, py, *, 0]+s,/xsty,xthick=2,ythick=2,thick=2,chars=1.3,xtit='!7k!3 [nm]', $
           ytit='intensity [DN]',charthick=2,/ylog, yr=[1.0,1d4]
     ee = ir_error(reform(si2_fit_gauss[px, py, *, 0]), /fuv, /dark)
     errplot, xx, si2_fit_gauss[px, py, *, 0]+s - ee, si2_fit_gauss[px, py, *, 0]+s + ee, thick=1, color=120
     oplot, xx, si2_fit_gauss[px, py, *, 0]+s, thick=3
     loadct,40,/silent  
     oplot, xx, si2_fit_gauss[px, py, *, 5]+s, color=245,thick=5.5
     ;oplot, [0,0]+ 140.151, [1,1000]
     ;oplot, [0,0]+ 139.9957, [1, 1000]
     loadct,0,/silent
     set_plot, 'x'
     !p.font=0
  endif


  
endwhile


stop

end
