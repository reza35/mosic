pro check_iris_o_fit, inpath

;+
;===============================================================
; procedure : check_iris_o_fit.pro
;  
; purpose : to check multi-line fit to O I and C I lines
;
; inpath :  the directory containing the save file 
;
; May 26, 2016 : created
;  
; R.Rezaei @ IAC                         e-mail:  rrezaei@iac.es      
;===============================================================
;-

device, decomposed=0
  
restore, inpath + 'iris_oi_analyzed_scan_0.sav', /v

window, 14, title='reduced chi-square statistics, O I + C I fits'
s=15.
x= reform(occ[*, *, 4]) & q=where(x ne 0.) & plot_histogram, x(q), xrange=[0, s]
loadct, 40,/silent
x= reform(occ_q[*, *, 6]) & q=where(x ne 0.) & plot_histogram, x(q), xrange=[0, s], oplott=245
loadct,0, /silent

tay = size(occ_fit)
window, 0, xs=tay[1]*3, ys=tay[2], tit='Continuum, chi-square, chi-square'
window, 2, xs=700, ys=550, tit='profile'

vsig = alog10(reform(occ[*,*,2])>2. < 100.)
;vsig = alog10(reform(occ_q[*,*,3])>1.)
osig = alog10(reform(occ_q[*,*,3])>1. < 100.)
tvscl, reform(occ[*,*,0])

finished = 0
while (finished eq 0) do begin
  cursor,x , y, /device, /change

; calculate position 
px = x mod (tay[1]) & py = y mod (tay[2])

loadct, 0,/silent
wset, 0
tvscl, reform(occ[*,*,0])
tvscl, vsig, tay[1], 0
tvscl, osig, tay[1]*2, 0
;-------------------------------
;---   draw cursor position
;-------------------------------
loadct,40,/silent
for i=0,2 do plots,[px,px]+tay[1]*i,[0,tay[2]],/device, color=250
for i=0,0 do plots,[0,tay[1]*3],[py,py]+tay[2]*i,/device, color=250

xyouts, 2, 2, num2string(fix(px))+','+num2string(fix(py)), /device, color=250,charsize=2
loadct,0,/silent

; quit: right mouse button
if !mouse.button eq 4 then finished = 1

;----------------------------------
;--  left button: inversion result
;----------------------------------
starter = 0
if !mouse.button eq 1  then starter = 1
if !mouse.button eq 2  then starter = 2
wait, .1
xx = (findgen(tay[3]))* 1.296e-3 + 135.37 - 0.007 ;dispersion_fuv

  if (starter eq 1) then begin
     print, px, py, occ[px,py,4], occ_q[px,py,6]
     wset, 2
     plot, occ_fit[px, py, *, 0],/xsty
     ee = 0.5 * sqrt(occ_fit[px, py, *, 0])
     errplot, occ_fit[px, py, *, 0] - ee, occ_fit[px, py, *, 0] + ee,color=220
     loadct,40,/silent
     oplot, occ_fit[px, py, *, 1], color=175,thick=2.5, linestyle=2
     oplot, occ_fit[px, py, *, 2], color=245,thick=2.5
     loadct,0,/silent
     print, '------------------------------------------------------------------------'
     wset, 0
  endif
  if (starter eq 2) then begin
     !p.font=-1
     set_plot, 'ps'
     device, filename='~/iris_o_prof_'+strtrim(px,1)+'_'+strtrim(py,1)+'.eps', /encapsulated, /color
     plot, xx, occ_fit[px, py, *, 0],/xsty,xthick=2,ythick=2,thick=2,chars=1.3,xtit='!7k!3 [nm]', $
           ytit='intensity [DN]',charthick=2, yr=[min(occ_fit[px, py, *, 0]) > (-5.0), max(occ_fit[px, py, *, 0])*1.05], /ystyle
;
;     /xsty,xthick=2,ythick=2,thick=2,chars=1.3,xtit='pixel', ytit='intensity [DN]',charthick=2
     ee = 0.5 * sqrt(occ_fit[px, py, *, 0])
     errplot, xx, occ_fit[px, py, *, 0] - ee, occ_fit[px, py, *, 0] + ee, thick=2
     loadct,40,/silent
     ;oplot, xx, occ_fit[px, py, *, 1], color=80,thick=4, linestyle=0
     oplot, xx, occ_fit[px, py, *, 2], color=245,thick=5
     loadct,0,/silent     
     set_plot, 'x'
     !p.font=0
     print, px, py, occ[px,py,4], occ_q[px,py,6]
  endif
endwhile


stop

end
