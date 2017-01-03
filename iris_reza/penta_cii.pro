function penta_cii, xarr, p

;+  
;===============================================================
; function : penta_cii.pro
; 
; purpose : performs a double Gaussian fit for each line of the 
;  two C II lines at 133 nm range and a single Gaussian fot to Ni II
;  133.52 nm line
; using only 13 degrees of freedom (width of the two C II lines are the same).
;
; the input profile is background subtracted.
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
 
  model = p[0]*exp(-((xarr - p[1])/p[2])^2)    ; C II 133.57, 1st
  model += p[3] * exp(-((xarr - p[4])/p[2])^2) ; C II 133.45, 1st
  
  model += p[5]*exp(-((xarr - p[6])/p[7])^2)   ; C II 133.57, 2nd
  model += p[8] * exp(-((xarr - p[9])/p[7])^2) ; C II 133.45, 2nd

  model += p[10]*exp(-((xarr - p[11])/p[12])^2) ; Ni II 133.520
  
  return, model
end
