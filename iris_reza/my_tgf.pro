function my_tgf, x, y, e, fit0, fit1, range0, range1, dlambda, good, double = double

;+
;===============================================================
; function : my_tgf.pro
;
; purpose : performs a penta-Gaussian fit with 7 free parameters to
; the Si IV 1403, O IV 1401/1399, and two other O IV and S IV lines.
; The idea is to have a controlled fit so parameters cannot return out
; of raneg values.
;  
; x/y/e : x/y of the profile and its error
; fit0/fit1 : initial guess values  
; range0/range1 : the range for parameters in fit0/fit1
; dlambda : scaling the range, usually 1.0d  
;  
; f = sum of five gaussian functions
;
; Dec 19, 2014 : created 
; Apr 26, 2016 : improved documentation
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

  ; x can be FIXED or LIMITED
  
  ; First Component Peak Intensity
  x1={limited:xlimited, limits:fit0[0] + range0[0]*[-dlambda, dlambda]*fit0[0], value:fit0[0], mpmaxstep:2.0d3} 
  ; First Component Line Center Position
  x2={limited:xlimited, limits: fit0[1] +range0[1]*[-dlambda, dlambda]*fit0[1], value:fit0[1], mpmaxstep:5.0d}
  ; First Component Gaussian Width
  x3={limited:xlimited, limits:fit0[2] + range0[2]*[-dlambda, dlambda]*fit0[2], value:fit0[2], mpmaxstep:5.0d-1} ;

  ; Second Component Peak Intensity
  x4={limited:xlimited, limits:[0.0d, fit1[0] + range1[0]*dlambda*fit1[0]], value:fit1[0],mpmaxstep:1.0d} 
  ; 3-rd Component Peak Intensity
  x5={limited:xlimited, limits:[0.0d, fit1[1] + range1[1]*dlambda*fit1[1]], value:fit1[1],mpmaxstep:1.0d} 
  ; 4-th Component Peak Intensity
  x6={limited:xlimited, limits:[0.0d, fit1[2] + range1[2]*dlambda*fit1[2]], value:fit1[2],mpmaxstep:1.0d} 
  ; 5-th Component Gaussian Width
  x7={limited:xlimited, limits:[0.0d, fit1[3] + range1[3]*dlambda*fit1[3]], value:fit1[3],mpmaxstep:1.0d} ;

  parinfo=[x1,x2,x3,x4, x5, x6, x7]
  param = parinfo.value
  ;print, parinfo.limits
  
  res = mpfitfun('pentagauss', x[good], y[good], e[good], param, parinfo=parinfo,$
    maxiter = 2000, dof = dof, bestnorm = bestnorm, yfit = yfit, double = double,status = status,  /quiet, perror=perr)

  if (n_elements(perr) eq 0)and(status ge 0.) then perr = res - res

result = {i1:res[0],p1:res[1],w1:res[2], i2:res[3], i3:res[4], i4:res[5], i5:res[6], fit:yfit, status:status,sigma:perr}
;print, res
  return, result
end

