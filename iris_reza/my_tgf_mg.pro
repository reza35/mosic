function my_tgf_mg, x, y, e, fit0, fit1, range0, range1, dlambda, good, double = double

;+
;===============================================================
; function : my_tgf_mg.pro
;
; purpose : performs a triple-Gaussian fit with 8 free parameters to
; the Mg II h or k lines.
;
; x/y/e : x/y of the profile and its error
; fit0/fit1 : initial guess values  
; range0/range1 : the range for parameters in fit0/fit1
; dlambda : scaling the range, usually 1.0d  
;  
; f = sum of three gaussian functions
;  
; Apr 26, 2015 : created 
; May 16, 2016 : slight modification 
;  
; R.Rezaei @ IAC                         e-mail:  rrezaei@iac.es      
;===============================================================
;-  

  xlimited=[1,1]
  
  fit0 = fit0 * 1.0d
  fit1 = fit1 * 1.0d
  range0 = range0 * 1.0d
  range1 = range1 * 1.0d

  ; limited is Boolian and comes with limits > mpfit parameters.
   
  x1={limited:xlimited, limits:fit0[0] + range0[0]*[-1.,1.]*fit0[0], value:fit0[0]} ; First Component Peak Intensity
  x2={limited:xlimited, limits:fit0[1] + range0[1]*[-1.,1.]*fit0[1], value:fit0[1]} ; First Component Line Center Position
  x3={limited:xlimited, limits:fit0[2] + range0[2]*[-1.,1.]*fit0[2], value:fit0[2]} ; First Component Gaussian Width
  x4={limited:xlimited, limits:fit1[0] + range1[0]*[-1.,1.]*fit1[0], value:fit1[0]} ; H2v Peak Intensity
  x5={limited:xlimited, limits:fit1[1] + range1[1]*[-1.,1.]*fit1[1], value:fit1[1]} ; H2v Line Center Position
  x6={limited:xlimited, limits:fit1[2] + range1[2]*[-1.,1.]*fit1[2], value:fit1[2]} ; H2v Gaussian Width
  x7={limited:xlimited, limits:fit1[3] + range1[3]*[-1.,1.]*fit1[3], value:fit1[3]} ; H2r Peak Intensity
  x8={limited:xlimited, limits:fit1[4] + range1[4]*[-1.,1.]*fit1[4], value:fit1[4]} ; H2r Line Center Position
  x9={limited:xlimited, limits:fit1[5] + range1[5]*[-1.,1.]*fit1[5], value:fit1[5]};; H2r Gaussian Width

  
  parinfo=[x1, x2, x3, x4, x5, x6, x7, x8, x9]
  param = parinfo.value
  ;print, parinfo.limits
  
  res = mpfitfun('tgf_mg', x[good], y[good], e[good], param, parinfo=parinfo,$
    maxiter = 2000, dof = dof, bestnorm = bestnorm, yfit = yfit, double = double,status = status,  /quiet, perror=perr)

  if (n_elements(perr) eq 0)and(status ge 0.) then perr = res - res

  result = {i1:res[0],p1:res[1],w1:res[2], i2:res[3], p2:res[4],w2:res[5],i3:res[6],p3:res[7], w3:res[8],$
            fit:yfit, status:status,sigma:perr}
  return, result
end

