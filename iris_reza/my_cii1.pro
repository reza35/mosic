function my_cii1, x, y, e, fit0, fit1, range0, range1, dlambda, good, double = double

;+
;===============================================================
; function : my_cii1.prp
; purpose : performs a diuble-Gaussian fit with 5 free parameters to
; the C II lines at 135 nm range
;
; x/y/e : x/y of the profile and its error
; fit0/fit1 : initial guess values  
; range0/range1 : the range for parameters in fit0/fit1
; dlambda : scaling the range, usually 1.0d  
;  
; f = x1*exp((x-x2)/x3)^2 + x4*exp((x-x5)/x3)^2
;  
; Feb 01, 2016 : created 
;  
; R.Rezaei @ IAC                         e-mail:  rrezaei@iac.es      
;===============================================================
;-  

  xlimited=[1,1]
  dlambda= 1.0d

  fit0 = fit0 * 1.0d
  fit1 = fit1 * 1.0d
  range0 = range0 * 1.0d
  range1 = range1 * 1.0d
 
  ; First Component Peak Intensity
  x1={limited:xlimited, limits:fit0[0] + range0[0]*[-dlambda, dlambda]*fit0[0], value:fit0[0], mpmaxstep:2.0d1} 
  ; First Component Line Center Position
  x2={limited:xlimited, limits:fit0[1] + range0[1]*[-dlambda, dlambda]*fit0[1], value:fit0[1], mpmaxstep:1.0d-1}
  ; First Component Gaussian Width
  x3={limited:xlimited, limits:fit0[2] + range0[2]*[-dlambda, dlambda]*fit0[2], value:fit0[2], mpmaxstep:1.9d} ;

  ; Second Component Peak Intensity
  x4={limited:xlimited, limits:fit1[0] + range1[0]*[-dlambda, dlambda]*fit1[0], value:fit1[0], mpmaxstep:2.0d1} 
  ; Second Component  line Center Position
  x5={limited:xlimited, limits:fit1[1] + range1[1]*[-dlambda, dlambda]*fit1[1], value:fit1[1], mpmaxstep:1.0d-1} 

  parinfo=[x1, x2, x3, x4, x5]
  param = parinfo.value
  ;print, reform(param)

  res = mpfitfun('double_cii', x[good], y[good], e[good], param, parinfo=parinfo, /quiet, nprint=0, errmsg=errmsg, $
    maxiter = 1000, dof = dof, bestnorm = bestnorm, yfit = yfit, double = double,status = status, perror=perr)

  if (n_elements(perr) eq 0)and(status ge 0.) then perr = res - res

  result = {i1:res[0],p1:res[1],w1:res[2], i2:res[3], p2:res[4], fit:yfit, status:status,sigma:perr}
  return, result
end

