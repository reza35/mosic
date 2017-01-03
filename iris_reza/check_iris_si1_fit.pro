pro check_iris_si1_fit, inpath

;+
;===============================================================
; procedure : check_iris_si1_fit.pro
;  
; purpose : to check multi-Gaussian fit to Si IV 139.37 nm
;
; inpath :  the directory containing the save file  
;
; May 27, 2016 : created
;  
; R.Rezaei @ IAC                         e-mail:  rrezaei@iac.es      
;===============================================================
;-

device, decomposed=0
  
restore, inpath + 'iris_si_analyzed_scan_0.sav', /v

window, 14, title='reduced chi-square statistics, Si IV 139.3 nm fit'

loadct,0,/silent
s = 10
x= reform(si_iv2[*, *, 4]) & q=where(x ne 0.) & plot_histogram, x(q), xrange=[0, s]
loadct,40,/silent
x= reform(si_iv2_v[*, *, 7]) & q=where(x ne 0.) & plot_histogram, x(q), xrange=[0, s],  oplott=70
x= reform(si_iv2_q[*, *, 7]) & q=where(x ne 0.) & plot_histogram, x(q), xrange=[0, s],  oplott=170
x= reform(si_iv2_z[*, *, 10]) & q=where(x ne 0.) & plot_histogram, x(q), xrange=[0, s], oplott=200
x= reform(si_iv2_f[*, *, 13]) & q=where(x ne 0.) & plot_histogram, x(q), xrange=[0, s], oplott=250
loadct,0,/silent


tay=size(cont)
window, 0, xs=tay[1]*3, ys=tay[2], tit='Continuum, Si IV amplitude, O IV amplitude'
window, 2, xs=700, ys=550, tit='profile'

vsig = alog10(reform(si_iv1[*,*,4])>1.)
osig = alog10(reform(si_iv1_v[*,*,7])>1.)
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
py = py

xyouts, 2, 2, num2string(fix(px))+','+num2string(fix(py)), /device, color=0,charsize=2

; quit: right mouse button
if !mouse.button eq 4 then finished = 1

; middle mouse button: switch between (changing)
if !mouse.button eq 2 then begin
   stop
endif

;----------------------------------
;--  left button: inversion result
;----------------------------------
starter = 0
if !mouse.button eq 1  then starter = 1
wait, .1

  if (starter eq 1) then begin
     print, px, py
     wset, 2
     plot, si1_fit_gauss[px, py, *, 0],/xsty
     ee = 0.5 * sqrt(si1_fit_gauss[px, py, *, 0])
     errplot, si1_fit_gauss[px, py, *, 0] - ee, si1_fit_gauss[px, py, *, 0] + ee, color=220
     loadct,40,/silent
     oplot, si1_fit_gauss[px, py, *, 1], color=175,thick=2.5, linestyle=2
     oplot, si1_fit_gauss[px, py, *, 2], color=245,thick=2.5
     loadct,0,/silent
     print, '---------------------------------'
     wset, 0
  endif
endwhile


stop

end
