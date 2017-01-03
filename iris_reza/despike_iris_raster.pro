function despike_iris_raster, inp, x1, x2

;+
;===============================================================
; function :  despike_iris_raster.pro
;  
; purpose:  to despike IRIS FUV spectra
; It uses a very simplistic despike algorithm:
; if distance between a point in the original and median filtered
; 2D raster is larger than 40 DN, then take it as cosmic ray.   
;
; inp : 2D FUV raster spectra
; x1,x2 : location of spectral lines (to be masked) 
;  
; 12/11/2014 : original
; 28/12/2015 : bug fix
;  
; R.Rezaei @ IAC                         e-mail:  rrezaei@iac.es      
;===============================================================
;-


np = n_elements(inp)

if ((x1 - 20.) lt 0.) then x1 = 20
if ((x2+20) ge (np-1)) then x2 = np-21

out = inp
u = median(median(inp, 15),3)
check = abs(u - inp)

bad =  where(check gt 40., count)

if (count ge 1) then begin 
     out(bad) = u(bad)
endif
out[x1-15:x1+15,*] = inp[x1-15:x1+15,*]
out[x2-15:x2+15,*] = inp[x2-15:x2+15,*]

return, out

end
