pro get_systematic_velocity, curve, output

;+
;===============================================================
; procedure : get_systematic_velocity.pro
; 
; purpose : to calculate a polynomial fit for orbital velocity modulations
;           interactively select an order for the temporal variation.
;           for large maps with nx=400, usually an order ~10 works.
;           
; Jun 17, 2015: original
;
; R.Rezaei @ IAC                         e-mail:  rrezaei@iac.es      
;===============================================================
;-

on_error, 2
satisfy = 'N'
device, decomposed=0

curve = median(curve, 3)
bbc = where(curve ne 0.)
q=where(curve eq 0., count) & if (count ge 1) then curve(q)= median(curve[bbc])

tay = size(curve)
window,32,/free

while ((satisfy ne 'y') AND (satisfy ne 'Y') ) do begin
  loadct,0, /silent
  plot, curve, /yno, charsize=1.3, psym=-1

  print, 'please enter the poly_fit order: 1=linear, ..., or negative to stop'
  read, order, PROMPT='Enter order: '
  xaux = findgen(tay[1])
  if (order gt 0) then res = poly_fit(xaux, curve, order, yfit=yfit)
  if (order eq 0) then yfit = fltarr(tay[1])+mean(curve)
  if (order lt 0) then stop
  loadct, 40, /silent
  oplot, xaux, yfit, color=245, thick=2
  print, 'if it is ok, pres y'
  satisfy = get_kbrd()
endwhile

wset, 32
loadct, 0, /silent
plot, curve
loadct, 40, /silent
oplot, yfit, color=245, thick=2
loadct, 0
wait, 2.
wdelete, 32
output = yfit - mean(yfit)

end  
