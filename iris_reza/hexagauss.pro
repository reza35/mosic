function hexagauss, xarr, p

;+
;===============================================================
; function : hexagauss.pro
; 
; purpose : hexa-Gaussian function.
; it is used to fit five lines in the range of Si I 1403 A, using 
; twelve degrees of freedom.
;
; We assume the line width of S IV is equal to the one of Si IV.
;
; calling sequence : 
; erg = my_hgf(x_prof, y_prof, error, fit0, fit1, range0, range1, dlambda, bbc, /double)
; 
; inputs : xarr - wavelength vector
;                    
;  p - parameters, [1st component peak intensity/centroid/width,
;      2nd component peak intensity/centroid/width,
;      3rd component peak intensity/centroid/width,
;      peak intensity of the 3rd, 4th and 5th components]
;      width of O IV is a free parameter
;
;  sigma of small lines are computed based on sigma of O IV 1401.16
;
; Dec 11, 2014 : original
; May 03, 2016 : improved
; May 07, 2016 : common block added for spectral dispersion
; Jun 04, 2016 : free dispersion parameter
; Jun 10, 2016 : several weak lines added with fix parameters
;
; R.Rezaei @ IAC                         e-mail:  rrezaei@iac.es      
;===============================================================
;-



common share_disp, z, u
common share_delta, dlam

df = z / 1.298


  lamb_x  = p[1] + dlam - 316.0d /df   ;  unknown  line at 1398.78 A
  lamb_ni  = p[1] + dlam - 298.0d /df  ;  weak  Ni II line at 1399.02 A
  lamb_o2 = p[1] + dlam - 235.60d /df  ;  the O IV line at 1399.78 A
  lamb_fe1 = p[1] + dlam - 222.60d /df ;  weak  Fe II line at 1399.96 A
  lamb_s1  = p[1] + dlam - 102.32d /df ;  weak  S I line at 1401.51 A
  lamb_s  = p[1] + dlam + 158.760d /df ;  weak  S IV line at 1404.81 A
  lamb_o3 = p[1] + dlam + 160.56d /df  ;  weak  O IV line at 1404.81 A
  lamb_fe2 = p[1] + dlam + 208.0d /df  ;  weak  Fe II line at 1405.61 A
  lamb_o4  = p[1] + dlam + 214.0d /df  ;  weak  O III line at 1406.45 A
  lamb_s2  = p[1] + dlam + 220.0d /df  ;  weak  S IV line at 1406.09 A
  lamb_o5  = p[1] + dlam + 244.0d /df  ;  weak  O IV line at 1407.38 A
  lamb_o6  = p[1] + dlam + 258.0d /df  ;  weak  O III line at 1407.6 A

  sig_s1 = p[8] * 0.25d      ; sigma of weak S I line at 1401.51 A
  sig_s2 = p[8] * 0.28d      ; sigma of weak S IV line at 1406.04 A
  sig_fe = p[8] * 0.17d      ; sigma of weak Fe II line at 1405.61 A
  sig_ni = p[8] * 0.33d      ; sigma of weak Ni II line at 1399.03 A
  sig_o3 = p[8] * 0.60d      ; sigma of O III line at 1405.80 A
  sig_x  = p[8] * 0.24d      ; sigma of unknown line at 1398.87 A

  model = p[0]*exp(-((xarr - p[1])/p[2])^2) + p[3]*exp(-((xarr - p[4])/p[5])^2) ; double Guassian fit to Si IV 1402.77 A

  model += p[6]*exp(-((xarr - p[7])/p[8])^2) ; single Guassian fit to O IV 1401.15 

  model += p[9]*exp(-((xarr - lamb_o2)/p[8])^2)  ; O IV 1399.97

  model += p[10]*exp(-((xarr - lamb_s)/p[2])^2)  ; S IV 1404.79

  model += p[11]*exp(-((xarr - lamb_o3)/p[8])^2) ; O IV 1404.82

  
  model += 0.45 * p[12] * exp(-((xarr - lamb_x)/sig_x)^2)   ;  1398.78
  model += 0.45 * p[12] * exp(-((xarr - lamb_ni)/sig_ni)^2)  ;  Ni II 1399.03
  model += 0.27 * p[12] * exp(-((xarr - lamb_fe1)/sig_fe)^2) ;  Fe II 1399.97
  
  model += p[12]*exp(-((xarr - lamb_s1)/sig_s1)^2) ; S I 1401.26

  model += 0.12 * p[12] * exp(-((xarr - lamb_fe2)/sig_fe)^2) ;  Fe II 1405.61
  model += 0.10 * p[6]  * exp(-((xarr - lamb_o4)/sig_o3)^2)  ;  O III 1405.80
  model += 0.20 * p[6]  * exp(-((xarr - lamb_s2)/sig_s1)^2)  ;  S IV 1406.04

  model += 0.05 * p[6] * exp(-((xarr - lamb_o5)/p[8])^2)   ;  O IV 1407.39
  model += 0.25 * p[6] * exp(-((xarr - lamb_o6)/sig_o3)^2) ;  O III 1407.70

  return, model
end
