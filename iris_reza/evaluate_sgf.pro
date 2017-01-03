function evaluate_sgf, init, erg

;+
;===============================================================
; function : evaluate_sgf.prp
;
; purpose : check the output of a single Gaussian fit and signals if
; the fit has failed.
;  
; May 13, 2016 : created 
;  
; R.Rezaei @ IAC                         e-mail:  rrezaei@iac.es      
;===============================================================
;-  

bad = 0.  
line_width = erg.w1
line_intensity = erg.i1

sigma = erg.sigma

if (line_intensity lt 0.0) then bad = 1.

if (line_width lt 0.2) then bad = 1.

if (max(sigma) gt 100.)or(sigma[1] gt 10.)or(sigma[3] gt 1.) then bad = 1.

return, bad
end  
