function lvg, x, y, e, fit0, fit1, range, dlambda, good, double = double

;+
;===============================================================
; function : lvg.pro
;
; purpose : performs a penta-Gaussian fit with ONE free parameters to
; the Si IV 1403, O IV 1401/1399, and two other O IV and S IV lines.
; It actually fine-tunes the line positions w.r.t. Si IV which needs
; fine-tunning in case of large velocity gradients.  
;  
; x/y/e : x/y of the profile and its error
; fit0/fit1 : initial guess values  
; range0/range1 : the range for parameters in fit0/fit1
; dlambda : scaling the range, usually 1.0d  
;  
; f = sum of five gaussian functions
;
; June 10, 2016 : created 
;  
; R.Rezaei @ IAC                         e-mail:  rrezaei@iac.es      
;===============================================================
;-  
  common shared_pars, mp
  xlimited=[1,1]
  dlambda *= 1.0d

  fit0 = fit0 * 1.0d
  fit1 = fit1 * 1.0d
  range = range * 1.0d

  ; x can be FIXED or LIMITED
  ;print, dlambda
  ; First Component Peak Intensity
  x1={fixed, value:fit0[0]} 
  ; First Component Line Center Position
  x2={fixed, value:fit0[1]}
  ; First Component Gaussian Width
  x3={fixed, value:fit0[2]} ;

  ; Second Component Peak Intensity
  x4={fixed, value:fit1[0]} 
  ; 3-rd Component Peak Intensity
  x5={fixed, value:fit1[1]} 
  ; 4-th Component Peak Intensity
  x6={fixed, value:fit1[2]} 
  ; 5-th Component Gaussian Width
  x7={fixed, value:fit1[3]} ;
  x8={limited:xlimited, limits:[-dlambda, dlambda]*fit1[4], value:0.2d, mpmaxstep:0.5d} ;
  
  mp = [reform(fit0), reform(fit1[0:3])]

  parinfo=[x8]
  param = parinfo.value
  ;print, parinfo.limits
  ;stop
  res = mpfitfun('pentagauss_lvg', x[good], y[good], e[good], param, parinfo=parinfo,$
    maxiter = 2000, dof = dof, bestnorm = bestnorm, yfit = yfit, double = double,status = status,  /quiet, perror=perr)

  if (n_elements(perr) eq 0)and(status ge 0.) then perr = res - res

  result = {p2:res[0], fit:yfit, status:status,sigma:perr}

  return, result
end

