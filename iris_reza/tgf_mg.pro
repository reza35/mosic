function tgf_mg, xarr, p
 
;+
;===============================================================
; function : tgf_mg.pro
; 
; purpose : triple  Gaussian function.
; it is used to fit three gaussian to any Mg II h/k profile
;
; inputs :
;  xarr : wavelength vector                    
;  p : triple gaussian parameters
;
; Dec 11, 2014:  original
; Dec 25, 2015:  depending if the dispersion is 1.28 or 2.56, the
;                wavelength offset should be carefully selected !!
;
; R.Rezaei @ IAC                         e-mail:  rrezaei@iac.es      
;===============================================================
;

  model  = p[0]*exp(-((xarr - p[1])/p[2])^2) ; the first Guassian fit
  model += p[3]*exp(-((xarr - p[4])/p[5])^2) ; the second Guassian fit
  model += p[6]*exp(-((xarr - p[7])/p[8])^2) ; the third Guassian fit

  return, model
end
