function my_sgf, x, y, e, fit0, range0, dlambda, good, double = double

;+
;===============================================================
; function : my_sgf.pro
;
; purpose : performs a single-Gaussian fit to the input line profile
; using 4 free parameters.
; The continuum level is calculated in advance and is only fine-tuned
; in a small range.
;
; x/y/e : x/y of the profile and its error
; fit0 : initial guess values  
; range0 : the range for parameters in fit0
; dlambda : scaling the range, usually 1.0d  
;
; The common block tells the program to increase the stepsize in the line
; width, which helps fitting the Mg II line profiles.  
;  
; f = x0 + x1*exp((x-x2)/x3)^2
;  
; Dec 19, 2014 : created 
; Apr 26, 2016 : improved documentation
; May 11, 2016 : continuum intensity can be negative  
; May 17, 2016 : common block to share the line tag   
;  
; R.Rezaei @ IAC                         e-mail:  rrezaei@iac.es      
;===============================================================
;-  
common share_line, die_linie                       

  xlimited=[1,1]
  if (n_elements(dlambda) eq 0) then dlambda= 1.0d

  fit0 = fit0 * 1.0d
  range0 = range0 * 1.0d
  
  peak = fit0[1] > 6.0d
  level = abs(fit0[0]) > 5.0d
  
  ; Peak Intensity
  x1={limited:xlimited, limits:fit0[1] + range0[1]*[-dlambda, dlambda]*peak, value:fit0[1], mpmaxstep:2.0d3} 
  ; Line Center Position
  x2={limited:xlimited, limits:fit0[2] + range0[2]*[-dlambda, dlambda]*fit0[2], value:fit0[2], mpmaxstep:5.0d}

  if (die_linie eq 'Mg') then begin
    ; Background Intensity (can be negative)
    x0={limited:xlimited, limits:fit0[0] + range0[0]*[-dlambda, dlambda]*level, value:fit0[0], mpmaxstep:5.0d}   
    ; Gaussian Width
    x3={limited:xlimited, limits:fit0[3] + range0[3]*[-dlambda, dlambda]*fit0[3], value:fit0[3], mpmaxstep:2.0d} 
  endif else begin
   ; Background Intensity (can be negative)
    x0={limited:xlimited, limits:fit0[0] + range0[0]*[-dlambda, dlambda]*abs(fit0[0]), value:fit0[0], mpmaxstep:1.0d}
   ; Gaussian Width
    x3={limited:xlimited, limits:fit0[3] + range0[3]*[-dlambda, dlambda]*fit0[3], value:fit0[3], mpmaxstep:5.0d-1}
  endelse   
  
  parinfo=[x0,x1,x2,x3]
  param = parinfo.value
  ;print, parinfo.limits
  
  res = mpfitfun('rgauss', x[good], y[good], e[good], param, parinfo=parinfo, /quiet, nprint=5, errmsg=errmsg, $ ;
    maxiter = 2000, dof = dof, bestnorm = bestnorm, double = double,status = status, perror=perr)

;  if (n_elements(perr) eq 0)and(status ge 0.) then perr = res - res
  if (n_elements(perr) eq 0)and(status ge 0.) then perr = res - res
  ;print, errmsg
  result = {b:res[0], i1:res[1], p1:res[2], w1:res[3], status:status, sigma:perr}
  return, result
end

