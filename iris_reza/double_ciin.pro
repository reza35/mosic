function double_ciin, xarr, p

;+  
;===============================================================
; function : double_ciin.pro
; 
; purpose : performs a double Gaussian fit.
; it is used to fit two C II lines at 133 nm range
; using only four degrees of freedom (width of the two C II lines are
; the same, amplitude ratio is locked).
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
; May 10, 2016: original
;
; R.Rezaei @ IAC                         e-mail:  rrezaei@iac.es      
;===============================================================
;-

  model = p[0]*0.82*exp(-((xarr - p[3])/p[2])^2) ; single Guassian fit to C II 133.57
  model += p[0] * exp(-((xarr - p[1])/p[2])^2)  ; C II 133.45 
 
  return, model
end
