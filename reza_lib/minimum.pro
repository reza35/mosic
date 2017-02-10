pro minimum, xin, arr, xout

;+
;===============================================================
; procedure: minimum.pro
;
; purpose: fit a parabola to a peak
;
; Aug 18, 2005 : created
; Nov 09, 2011 : updated 
;
; R.Rezaei @ KIS                         e-mail:  rrezaei@iac.es      
;===============================================================
;-       

s = size(arr)
nb = s[1]
xout = fltarr(9)

nc = 7
xf = findgen(nc) + xin - 3
if (max(abs(arr)) le 0.02)or(nb lt 8) then xf = findgen(5) + xin - 2.0d
yf = arr[xf]

  res = poly_fit(xf, yf, 2, /double)
  yres = res[0] + res[1]*xf + res[2] * xf^2
  ;oplot, xf, yres;, linestyle=1
 
  xout[0:2] = res
  xout[3] = -res[1]/(2*res[2])
  xout[4] = min(yres)
  xout[5] = min(xf)
  xout[6] = res[0] + res[1]*xout[5] + res[2] * xout[5]^2
  xout[7] = max(xf)
  xout[8] = res[0] + res[1]*xout[7] + res[2] * xout[7]^2

if (total(~finite(xout)) gt 0.) then begin

arr = smooth(arr, 5)
s = size(arr)
nb = s[1]

nc = 7
xf = findgen(nc) + xin - 3
if (max(abs(arr)) le 0.02) then xf = findgen(5) + xin - 2

yf = arr[xf]
;print,xin, min(xf), max(xf)

  res = poly_fit(xf, yf, 2, /double)
  yres = res[0] + res[1]*xf + res[2] * xf^2
  ;oplot, xf, yres;, linestyle=1
 
  xout[0:2] = res
  xout[3] = -res[1]/(2*res[2])
  xout[4] = min(yres)
  xout[5] = min(xf)
  xout[6] = res[0] + res[1]*xout[5] + res[2] * xout[5]^2
  xout[7] = max(xf)
  xout[8] = res[0] + res[1]*xout[7] + res[2] * xout[7]^2

endif


end
