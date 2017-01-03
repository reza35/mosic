function double_cii, xarr, p

;+
;===============================================================
; function : double_cii.pro
; 
; purpose : performs a double Gaussian fit.
; it is used to fit two C II lines at 133 nm range
; using only five degrees of freedom (width of the two C II lines are the same).
;
; the input profile is background subtracte.
;
; calling sequence : erg = double_cii(xarr,p)
; 
; inputs : xarr - wavelength vector
;                    
;  p - parameters, [1st component peak intensity/centroid/width, 
;      2nd component peak intensity/width, and the peak intensity for the 3rd component]
;
; Jun 17, 2015: original
; May 07, 2016: common block added for spectral dispersion
;
; R.Rezaei @ IAC                         e-mail:  rrezaei@iac.es      
;===============================================================
;-
 

  model = p[3]*exp(-((xarr - p[4])/p[2])^2) ; single Guassian fit to C II 133.57
  model += p[0] * exp(-((xarr - p[1])/p[2])^2)  ; C II 133.45 
 
  return, model
end
