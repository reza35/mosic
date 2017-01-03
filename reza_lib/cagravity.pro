function grav, a1, a2, a3, arr

g1 = total(abs(arr(a1:a2)),1) - total(abs(arr(a2:a3)),1)
return, g1

end


pro cagravity, a1, a2, arr, gravit, amplitude

;+
;===============================================================
; procedure: cagravity.pro
;
; purpose: calculate center-of-gravity wavelength for 
;          an emission profile
;
; May 01, 2007 : created
;
; R.Rezaei @ KIS                         e-mail:  rrezaei@iac.es      
;===============================================================
;-       

  on_error, 2
  
  num = n_elements(arr)
  if (num lt 6) then begin
     print, 'array is too short.'
     return
  endif   
  k1 = a2-a1+1
  s = findgen(k1)+a1
  r = fltarr(k1)
  for i=0, k1-1 do r[i] = grav(a1, i+a1, a2, arr)
  x = min(abs(r), m1)

  x = findgen(5) + m1- 2.
  y = r(x)

  m1 += float(a1)
  x += float(a1)

  res = poly_fit(x, y, 1, /double)
  erg = where(finite(res,/nan), count)

  if (count eq 0)and(min(abs(res[1])) gt 1.0d-3) then begin
    gravit = - res[0]/res[1] 
  endif else begin
    res = poly_fit(x, smooth(y,3,/edge_truncate), 1, /double)
    erg = where(finite(res,/nan), count)
    if (count eq 0)and(min(abs(res[1])) gt 1.0d-3) then begin
       gravit = - res[0]/res[1] 
    endif else begin
       gravit = float(num)/2.0
    endelse   
  endelse   
    
  y = arr(x)
  res = poly_fit(x, y, 2, /double)
  amplitude = res[0] + res[1]*gravit + res[2]*(gravit)^2
  
end
