function my_cii1n, x, y, e, fit0, range0, dlambda, good, double = double

;+
;===============================================================
; function : my_cii1n.prp
; purpose : performs a double-Gaussian fit with 4 free parameters to
; the C II lines at 135 nm range
;
; x/y/e : x/y of the profile and its error
; fit0 : initial guess values  
; range0 : the range for parameters in fit0/fit1
; dlambda : scaling the range, usually 1.0d  
;  
; f = x1*0.82 *exp((x-x2)/x3)^2 + x1*exp((x-x4)/x3)^2 
;  
; Feb 01, 2016 : created 
;  
; R.Rezaei @ IAC                         e-mail:  rrezaei@iac.es      
;===============================================================
;-  

  xlimited=[1,1]
  dlambda= 1.0d

  fit0 = fit0 * 1.0d
  range0 = range0 * 1.0d
 
  ; First Component Peak Intensity
  x1={limited:xlimited, limits:fit0[0] + range0[0]*[-dlambda, dlambda]*fit0[0], value:fit0[0], mpmaxstep:2.0d} 
  ; First Component Line Center Position
  x2={limited:xlimited, limits:fit0[1] + range0[1]*[-dlambda, dlambda]*fit0[1], value:fit0[1], mpmaxstep:1.0d-1}
  ; First Component Gaussian Width
  x3={limited:xlimited, limits:fit0[2] + range0[2]*[-dlambda, dlambda]*fit0[2], value:fit0[2], mpmaxstep:1.9d} ;

  ; Second Component  line Center Position
  x4={limited:xlimited, limits:fit0[3] + range0[3]*[-dlambda, dlambda]*fit0[3], value:fit0[3], mpmaxstep:1.0d-1} 

  parinfo=[x1, x2, x3, x4]
  param = parinfo.value
  ;print, reform(param)

  res = mpfitfun('double_ciin', x[good], y[good], e[good], param, parinfo=parinfo, /quiet, nprint=0, errmsg=errmsg, $
    maxiter = 1000, dof = dof, bestnorm = bestnorm, yfit = yfit, double = double,status = status, perror=perr)

  if (n_elements(perr) eq 0)and(status ge 0.) then perr = res - res

  result = {i1:res[0],p1:res[1],w1:res[2], p2:res[3], fit:yfit, status:status,sigma:perr}
  return, result
end

