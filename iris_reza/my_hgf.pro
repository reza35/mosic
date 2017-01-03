function my_hgf, x, y, e, fit0, fit1, range0, range1, dlambda, good, double = double

;+
;===============================================================
; function : my_hgf.pro
;  
; purpose : performs a hexa-Gaussian fit with 13 free parameters to
; the Si IV 1403, O IV 1401/1399, and two other O IV and S IV lines.
; It is very similar in approch to my_pgf.pro, and the main difference
; is to have a double Gaussian fit to Si IV line.
;
; x/y/e : x/y of the profile and its error
; fit0/fit1 : initial guess values  
; range0/range1 : the range for parameters in fit0/fit1
; dlambda : scaling the range, usually 1.0d  
;  
; f = sum of six gaussian functions
;  
; May 03, 2016 : created 
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
  x1={limited:xlimited, limits:fit0[0] + range0[0]*[-dlambda, dlambda]*fit0[0], value:fit0[0]}
  ; First Component Line Center Position
  x2={limited:xlimited, limits:fit0[1] + range0[1]*[-dlambda, dlambda]*fit0[1], value:fit0[1]} 
  ; First Component Gaussian Width
  x3={limited:xlimited, limits:fit0[2] + range0[2]*[-dlambda, dlambda]*fit0[2], value:fit0[2]}
  ; Second Component Peak Intensity
  x4={limited:xlimited, limits:fit0[3] + range0[3]*[-dlambda, dlambda]*fit0[3], value:fit0[3]}
  ; Second Component Line Center Position
  x5={limited:xlimited, limits:fit0[4] + range0[4]*[-dlambda, dlambda]*fit0[4], value:fit0[4]} 
  ; Second Component Gaussian Width
  x6={limited:xlimited, limits:fit0[5] + range0[5]*[-dlambda, dlambda]*fit0[5], value:fit0[5]}

  ; Third Component Peak Intensity
  x7={limited:xlimited, limits:fit1[0] + range1[0]*[-dlambda, dlambda]*fit1[0], value:fit1[0]} 
  ; Third Component Line Center Position
  x8={limited:xlimited, limits:fit1[1] + range1[1]*[-dlambda, dlambda]*fit1[1], value:fit1[1]}
  ; Third Component Gaussian Width
  x9={limited:xlimited, limits:fit1[2] + range1[2]*[-dlambda, dlambda]*fit1[2], value:fit1[2]}

  ; Fourth component Peak Intensity
  x10={limited:xlimited, limits:fit1[3] + range1[3]*[-dlambda, dlambda]*fit1[3], value:fit1[3]} 
  ; Fifth component Peak Intensity
  x11={limited:xlimited, limits:fit1[4] + range1[4]*[-dlambda, dlambda]*fit1[4], value:fit1[4]} 
  ; Sixth component Peak Intensity
  x12={limited:xlimited, limits:fit1[5] + range1[5]*[-dlambda, dlambda]*fit1[5], value:fit1[5]} 
  ; Seventh component Peak Intensity (for Fe II and Ni II lines)
  x13={limited:xlimited, limits:fit1[6] + range1[6]*[-dlambda, dlambda]*fit1[6], value:fit1[6]} 

  parinfo=[x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13]
  param = parinfo.value

  res = mpfitfun('hexagauss', x[good], y[good], e[good], param, parinfo=parinfo,$
    maxiter = 2000, dof = dof, bestnorm = bestnorm, yfit = yfit, double = double,status = status,  /quiet, perror=perr)

  if (n_elements(perr) eq 0)and(status ge 0.) then perr = res - res

  result = {i1:res[0],p1:res[1],w1:res[2], i2:res[3], p2:res[4], w2:res[5], $
            i3:res[6], p3:res[7], w3:res[8], $
            i4:res[9], i5:res[10], i6:res[11], i7:res[12],$
            fit:yfit, status:status,sigma:perr}

  return, result
end

