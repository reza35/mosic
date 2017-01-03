;;             NAME : doublegauss
; 
;          PURPOSE : Double Gaussian function. Called by dgf_1lp.pro.
;
; CALLING SEQUENCE : dg =  doublegauss(xarr,p)
; 
;           INPUTS : xarr - wavelength vector
;                    
;                    p - double gaussian parameters, [background, 1st component peak intensity/centroid/width, 2nd component peak intensity/centroid/width]
;
;          HISTORY : Written by Hui Tian at CfA, April 4, 2013

function doublegauss, xarr, p
 
  model = p[0] + p[1]*exp(-((xarr - p[2])/p[3])^2) + p[4]*exp(-((xarr - p[5])/p[6])^2)

  return, model
end