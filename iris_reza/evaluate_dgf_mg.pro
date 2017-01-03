function evaluate_dgf_mg, init, erg

;+
;===============================================================
; function : evaluate_dgf_mg.prp
;
; purpose : check the output of a double Gaussian fit and signal if
; the fit has failed.
;  
; May 13, 2016 : created 
;  
; R.Rezaei @ IAC                         e-mail:  rrezaei@iac.es      
;===============================================================
;-  

bad = 0.  
line_width1 = erg.w1
line_width2 = erg.w2

sigma = erg.sigma

if (line_width1 lt 0.2)or(line_width2 lt 0.2) then bad = 1.

if (max(sigma) gt 100.)or(sigma[0] gt 10.)or(sigma[2] gt 1.) then bad = 1.

if (sigma[3] gt 10.)or(sigma[5] gt 1.) then bad = 1.

return, bad
end  
