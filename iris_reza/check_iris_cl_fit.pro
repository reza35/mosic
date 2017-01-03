pro check_iris_cl_fit, infile

;+
;===============================================================
; procedure : check_iris_cl_fit.pro
;  
; purpose : to check the single-Gaussian fit to Cl I line
;
; infile :  the save file containing results of the Gaussian fitting  
;
; May 26, 2016 : created
;  
; R.Rezaei @ IAC                         e-mail:  rrezaei@iac.es      
;===============================================================
;-

device, decomposed=0
  
restore, infile

tay=size(cl[*,*,0])
window, 0, xs=tay[1]*3, ys=tay[2], tit='Continuum, Amplitude'
window, 2, xs=700, ys=550, tit='profile'

vsig = alog10(reform(cl[*,*,0])>1.)
osig = alog10(reform(cl[*,*,1])>1.)
tvscl, reform(occ[*,*,1])

finished = 0
while (finished eq 0) do begin
  cursor,x , y, /device, /change

; calculate position 
px = x mod (tay[1]) & py = y mod (tay[2])

loadct, 0,/silent
wset, 0
tvscl, reform(cl[*,*,0])
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
     print, px, py, cl[px,py,4]
     wset, 2
     plot, cl_fit_gauss[px, py, *, 0],/xsty
     ee = 0.5 * sqrt(cl_fit_gauss[px, py, *, 0])
     errplot, cl_fit_gauss[px, py, *, 0] - ee, cl_fit_gauss[px, py, *, 0] + ee,color=220
     loadct,40,/silent
     oplot, cl_fit_gauss[px, py, *, 1], color=245,thick=2.5
     loadct,0,/silent
     print, '------------------------------------------------------------------------'
     wset, 0
  endif
endwhile


stop

end
