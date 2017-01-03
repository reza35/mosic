function my_dgf_mg, x, y, e, fit0, fit1, range0, range1, dlambda, good, double = double

;+
;===============================================================
; function : my_dgf.pro

; purpose : performs a double-Gaussian fit to the input line profile
; using 6 free parameters.
; It is very similar in approch to my_dgf.pro, and the main difference
; is to have no parameter for the continuum
;
; x/y/e : x/y of the profile and its error
; fit0/fit1 : initial guess values  
; range0/range1 : the range for parameters in fit0/fit1
; dlambda : scaling the range, usually 1.0d  
;  
; f = x0*exp((x-x1)/x2)^2 + x3*exp((x-x4)/x5)^2
;  
; Dec 19, 2014 : created 
; Apr 26, 2016 : improved documentation
; May 11, 2016 : continuum intensity can be negative
; May 13, 2016 : range of parameters updated  
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
  
  peak = fit0[1] > 9.0d

  ; First Component Peak Intensity
  x1={limited:xlimited, limits:fit0[0] + range0[0]*[-dlambda, dlambda]*peak, value:fit0[0], mpmaxstep:2000.0d}
  ; First Component Line Center Position
  x2={limited:xlimited, limits:fit0[1] + range0[1]*[-dlambda, dlambda]*fit0[1], value:fit0[1], mpmaxstep:3.0d}
  ; First Component Gaussian Width
  x3={limited:xlimited, limits:fit0[2] + range0[2]*[-dlambda, dlambda]*fit0[2], value:fit0[2], mpmaxstep:0.2d}

  peak = fit1[0] > 9.
  ; Second Component Peak Intensity
  x4={limited:xlimited, limits:fit1[0] + range1[0]*[-dlambda, dlambda]*fit1[0], value:fit1[0], mpmaxstep:1.0d2}
  ; Second Component Line Center Position
  x5={limited:xlimited, limits:fit1[1] + range1[1]*[-dlambda, dlambda]*fit1[1], value:fit1[1], mpmaxstep:3.0d}
  ; Second Component Gaussian Width
  x6={limited:xlimited, limits:fit1[2] + range1[2]*[-dlambda, dlambda]*fit1[2], value:fit1[2], mpmaxstep:0.2d} 

  parinfo=[x1, x2, x3, x4, x5, x6]
  param = parinfo.value
  ;print
  ;print, parinfo.value
  ;print, '++++++++++++++++++++++++++++++++++++++++++++++'
  ;print, parinfo.limits
  ;print, '----------------------------------------------'
  res = mpfitfun('doublegauss_mg', x[good], y[good], e[good], param, parinfo=parinfo, nprint=1, errmsg=errmsg, $
    maxiter = 2000, dof = dof, bestnorm = bestnorm, yfit = yfit, double = double,status = status, /quiet, perror=perr)

  if (n_elements(perr) eq 0)and(status ge 0.) then perr = res - res

  result = {i1:res[0], p1:res[1], w1:res[2], i2:res[3], p2:res[4], w2:res[5], fit:yfit, status:status,sigma:perr}
  return, result
end

