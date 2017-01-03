function quad_cii, xarr, p

;+
;===============================================================
; function : quad_cii.pro
; 
; purpose : performs a double Gaussian fit for each line of the .
;  two C II lines at 133 nm range
; using only ten degrees of freedom (width of the two C II lines are the same).
;
; the input profile is background subtracte.
;
; inpts : xarr - wavelength vector
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

  model = p[0]*exp(-((xarr - p[1])/p[2])^2)    ; C II 133.57, 1st
  model += p[3] * exp(-((xarr - p[4])/p[2])^2) ; C II 133.45, 1st
  
  model += p[5]*exp(-((xarr - p[6])/p[7])^2)   ; C II 133.57, 2nd
  model += p[8] * exp(-((xarr - p[9])/p[7])^2) ; C II 133.45, 2nd
 
  return, model
end
