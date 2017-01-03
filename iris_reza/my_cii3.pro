function my_cii3, x, y, e, fit0, fit1, range0, range1, dlambda, good, double = double

;+
;===============================================================
; function : my_cii3.prp
;
; purpose : performs a penta-Gaussian fit with 13 free parameters to
; both C II lines at 135 nm range and teh Ni II 133.52 nm
;
; x/y/e : x/y of the profile and its error
; fit0/fit1 : initial guess values  
; range0/range1 : the range for parameters in fit0/fit1
; dlambda : scaling the range, usually 1.0d  
;  
; f = sum of five Gaussians
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
  x1={limited:xlimited, limits:fit0[0] + range0[0]*[-dlambda, dlambda]*fit0[0], value:fit0[0], mpmaxstep:2.0d} 
  ; First Component Line Center Position
  x2={limited:xlimited, limits:fit0[1] + range0[1]*[-dlambda, dlambda]*fit0[1], value:fit0[1], mpmaxstep:1.0d-1}
  ; First Component Gaussian Width
  x3={limited:xlimited, limits:fit0[2] + range0[2]*[-dlambda, dlambda]*fit0[2], value:fit0[2], mpmaxstep:1.9d} ;
  ; Second Component Peak Intensity
  x4={limited:xlimited, limits:fit0[3] + range0[3]*[-dlambda, dlambda]*fit0[3], value:fit0[3], mpmaxstep:2.0d} 
  ; Second Component  line Center Position
  x5={limited:xlimited, limits:fit0[4] + range0[4]*[-dlambda, dlambda]*fit0[4], value:fit0[4], mpmaxstep:1.0d-1} 
  ; Third Component Peak Intensity
  x6={limited:xlimited, limits:fit0[5] + range0[5]*[-dlambda, dlambda]*fit0[5], value:fit0[5], mpmaxstep:2.0d} 
  ; third Component Line Center Position
  x7={limited:xlimited, limits:fit0[6] + range0[6]*[-dlambda, dlambda]*fit0[6], value:fit0[6], mpmaxstep:1.0d-1}
  ; Third Component Gaussian Width
  x8={limited:xlimited, limits:fit0[7] + range0[7]*[-dlambda, dlambda]*fit0[7], value:fit0[7], mpmaxstep:1.9d} ;
  ; Fourth Component Peak Intensity
  x9={limited:xlimited, limits:fit0[8] + range0[8]*[-dlambda, dlambda]*fit0[8], value:fit0[8], mpmaxstep:2.0d} 
  ; Fourth Component  line Center Position
  x10={limited:xlimited,limits:fit0[9] + range0[9]*[-dlambda, dlambda]*fit0[9], value:fit0[9], mpmaxstep:1.0d-1} 

  ; Ni II Peak intensity
  x11={limited:xlimited, limits:fit1[0] + range1[0]*[-dlambda, dlambda]*fit1[0], value:fit1[0], mpmaxstep:1.0d-1} 
  ; Ni II line center position
  x12={limited:xlimited, limits:fit1[1] + range1[1]*[-dlambda, dlambda]*fit1[1], value:fit1[1], mpmaxstep:5.0d-1} 
  ; Ni II line width
  x13={limited:xlimited, limits:fit1[2] + range1[2]*[-dlambda, dlambda]*fit1[2], value:fit1[2], mpmaxstep:1.0d-1}

  
  parinfo=[x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13]
  param = parinfo.value


  res = mpfitfun('penta_cii', x[good], y[good], e[good], param, parinfo=parinfo, /quiet, nprint=0, errmsg=errmsg, $
    maxiter = 1000, dof = dof, bestnorm = bestnorm, yfit = yfit, double = double,status = status, perror=perr)

  if (n_elements(perr) eq 0)and(status ge 0.) then perr = res - res

  result = {i1:res[0],p1:res[1],w1:res[2], i2:res[3], p2:res[4], i3:res[5],p3:res[6],w2:res[7], i4:res[8], p4:res[9], $
            i5:res[10], p5:res[11], w3:res[12], $  ; Ni II parameters
            fit:yfit, status:status, sigma:perr}
  return, result
end

