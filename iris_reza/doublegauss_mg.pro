function doublegauss_mg, xarr, p

;+
;===============================================================
; function : doublegauss_mg.pro
; 
; purpose : performs a double Gaussian fit to teh Mg II line core
;
; Jan 07, 2016: original
; May 07, 2016: common block added for spectral dispersion
;
; R.Rezaei @ IAC                         e-mail:  rrezaei@iac.es      
;===============================================================
;-

  model = p[0]*exp(-((xarr - p[1])/p[2])^2) + p[3]*exp(-((xarr - p[4])/p[5])^2)
  return, model
end
