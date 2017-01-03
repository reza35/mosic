function triple_occ, xarr, p

;+
;===============================================================
; function : triple_occ.pro
; 
; purpose : performs a triple Gaussian fit.
; it is used to fit three lines : O I 1356, C I 1354.288, & C I 1355.84
; using only six degrees of freedom (width of the two C I lines are the same).
;
; the input profile is background subtracte.
;
; inputs :
;  xarr : wavelength vector            
;  p : parameters, [1st component peak intensity/centroid/width, 
;      2nd component peak intensity/width, and the peak intensity for the 3rd component]
;
; Jun 17, 2015: original
; May 07, 2016: common block added for spectral dispersion
; May 31, 2016: dispersion of 5 mA/pixel was added.
;
; R.Rezaei @ IAC                         e-mail:  rrezaei@iac.es      
;===============================================================
;-  

 
common share_disp, z, u

; when dispersion = 5.1
  if (z gt 5.)and(z lt 5.3) then begin
     lamb_c1 = p[1] + 4.80d     ; in pixel unit  >> strong C line in red wing of O I line
     lamb_c2 = p[1] - 25.15      ; in pixel unit  >> weak  C line in blue wing of O I line
  endif
; when dispersion = 2.56
  if (z gt 2.4)and(z lt 2.6) then begin
     lamb_c1 = p[1] + 9.60d     ; in pixel unit  >> strong C line in red wing of O I line
     lamb_c2 = p[1] - 50.30      ; in pixel unit  >> weak  C line in blue wing of O I line
  endif
  if (z gt 1.2)and(z lt 1.3) then begin
; when dispersion = 1.28  
     lamb_c1 = p[1] + 19.20d    ; in pixel unit  >> strong C line in red wing of O I line
     lamb_c2 = p[1] - 100.60    ; in pixel unit  >> weak   C line in blue wing of O I line
  endif
  

  model = p[0]*exp(-((xarr - p[1])/p[2])^2) ; single Guassian fit to O I line at 1356 A
  model += p[3] * exp(-((xarr - lamb_c1)/p[4])^2)  ; C I 
  model += p[5] * exp(-((xarr - lamb_c2)/p[4])^2)  ; C I 
 
  return, model
end
