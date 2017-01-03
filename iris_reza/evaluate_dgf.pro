function evaluate_dgf, init, erg

;+
;===============================================================
; function : evaluate_dgf.prp
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

if (max(sigma) gt 100.)or(sigma[1] gt 10.)or(sigma[3] gt 1.) then bad = 1.

if (sigma[4] gt 10.)or(sigma[6] gt 1.) then bad = 1.

return, bad
end  
