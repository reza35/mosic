function pentagauss, xarr, p

;+
;===============================================================
; function : pentagauss.pro
; 
; purpose : penta-Gaussian function.
; it is used to fit five lines in the range of Si I 1403 A, using only
; seven degrees of freedom like a double gaussian fit
; We assume the line width of S IV is equal to the one of Si IV.
;
; inputs : xarr - wavelength vector
;                    
;  p - parameters, [1st component peak intensity/centroid/width, 
;      peak intensity of the second, 3rd, 4th and 5th components]
;      width of O IV = 1.22 x width of Si IV
;
; Dec 11, 2014:  original
; Dec 25, 2015:  depending if the dispersion is 1.28 or 2.56, the
;                wavelength offset should be carefully selected !!
; Jun 04, 2016:  free dispersion parameter
;
;
; R.Rezaei @ IAC                         e-mail:  rrezaei@iac.es      
;===============================================================
;

common share_disp, z, u

df = z / 1.298

  lamb_s = p[1] + 158.760d /df ; in pixel unit  >> weak   S IV line at 1404.79 A
  lamb_o3 = p[1] + 160.56d /df ; in pixel unit  >> weak   O IV line at 1404.82 A
  lamb_o2 = p[1] - 235.60d /df ; in pixel unit  >> weak   O IV line at 1399.97 A
  lamb_o1 = p[1] - 126.32d /df ; in pixel unit  >> strong O IV line at 1401.16 A

  sig_o1 = p[2] * 1.22d      ; sigma of strong O IV line at 1401.16 A


  model = p[0]*exp(-((xarr - p[1])/p[2])^2) ; single Guassian fit to Si IV 1403

  model += p[3]*exp(-((xarr - lamb_o1)/sig_o1)^2)  ; O IV 1401.16 

  model += p[4]*exp(-((xarr - lamb_o2)/sig_o1)^2)  ; O IV 1399.97

  model += p[5]*exp(-((xarr - lamb_s)/p[2])^2)     ; S IV 1404.79

  model += p[6]*exp(-((xarr - lamb_o3)/sig_o1)^2)  ; O IV 1404.82

  
  return, model
end
