function pentagauss2, xarr, p

;+
;===============================================================
; function : pentagauss2.pro
; 
; purpose : penta-Gaussian function.
; it is used to fit five lines in the range of Si I 1403 A, using 
; nine degrees of freedom.
; We assume the line width of S IV is equal to the one of Si IV.
;
; calling sequence : 
; erg = my_pgf(x_prof, y_prof, error, fit0, fit1, range0, range1, dlambda, bbc, /double)
; 
; inputs : xarr - wavelength vector
;                    
;  p - parameters, [1st component peak intensity/centroid/width,
;      2nd component peak intensity/centroid/width, 
;      peak intensity of the 3rd, 4th and 5th components]
;      width of O IV is a free parameter
;
;  sigma of small lines are computed based on sigma of stronger lines
;
; Dec 11, 2014:  original
; Dec 25, 2015:  depending if the dispersion is 1.28 or 2.56, the
;                wavelength offset should be carefully selected !!
; Apr 26, 2016:  improved
; May 04, 2016:  free dispersion parameter
; Jun 10, 2016:  several small lines added  with fix parameters
;
; R.Rezaei @ IAC                         e-mail:  rrezaei@iac.es      
;===============================================================
;-   


common share_disp, z, u
common share_delta, dlam

df = z / 1.298

  lamb_s  = p[1] + dlam + 158.760d / df     ; in pixel unit  >> weak   S IV line at 1404.79 A
  lamb_o3 = p[1] + dlam + 160.56d / df ; in pixel unit  >> weak   O IV line at 1404.82 A
  lamb_o2 = p[1] + dlam - 235.60d / df; in pixel unit  >> weak   O IV line at 1399.97 A

  lamb_x  = p[1] + dlam - 316.0d /df ; in pixel unit  >> unknown (?) line at 1398.78 A
  lamb_ni  = p[1] + dlam - 298.0d /df ; in pixel unit  >> weak  Ni II line at 1399.03 A
  lamb_fe1 = p[1] + dlam - 222.60d /df ; in pixel unit  >> weak  Fe II line at 1399.97 A
  lamb_s1  = p[1] + dlam - 102.32d /df ; in pixel unit  >> weak S I line at 1401.51 A
  lamb_fe2 = p[1] + dlam + 208.0d /df ; in pixel unit  >> weak  Fe II line at 1405.61 A
  lamb_o4  = p[1] + dlam + 214.0d /df ; in pixel unit  >> weak  O III line at 1405.80 A
  lamb_s2  = p[1] + dlam + 220.0d /df ; in pixel unit  >> weak  S IV line at 1406.04 A

  lamb_o5  = p[1] + dlam + 244.0d /df ; in pixel unit  >> weak  O IV line at 1407.39 A
  lamb_o6  = p[1] + dlam + 258.0d /df ; in pixel unit  >> weak  O III line at 1407.70 A

  sig_s1 = p[2] * 0.30d      ; sigma of S I line at 1401.51 A
  sig_o3 = p[5] * 0.90d      ; sigma of O III line at 1405.80 A
  sig_s2 = p[2] * 0.35d      ; sigma of S IV line at 1406.04 A
  sig_fe = p[2] * 0.21d      ; sigma of Fe II line at 1405.61 A
  sig_ni = p[2] * 0.31d      ; sigma of Ni II line at 1399.03 A


  model = p[0]*exp(-((xarr - p[1])/p[2])^2)  ; single Guassian fit to Si IV 1403

  model += p[3]*exp(-((xarr - p[4])/p[5])^2) ; single Guassian fit to O IV 1401.16 

  model += p[6]*exp(-((xarr - lamb_o2)/p[5])^2)  ; O IV 1399.97

  model += p[7]*exp(-((xarr - lamb_s)/p[2])^2)  ; S IV 1404.79

  model += p[8]*exp(-((xarr - lamb_o3)/p[5])^2)  ; O IV 1404.82


  model += 0.55 * p[9] * exp(-((xarr - lamb_x)/sig_ni)^2)   ;  1398.78
  model += 0.45 * p[9] * exp(-((xarr - lamb_ni)/sig_ni)^2)  ;  Ni II 1399.03
  model += 0.27 * p[9] * exp(-((xarr - lamb_fe1)/sig_fe)^2) ;  Fe II 1399.97

  model += p[9] * exp(-((xarr - lamb_s1)/sig_s1)^2) ; S I 1401.26 

  model += 0.12 * p[9] * exp(-((xarr - lamb_fe2)/sig_fe)^2) ;  Fe II 1405.61
  model += 0.10 * p[3] * exp(-((xarr - lamb_o4)/sig_o3)^2)  ;  O III 1405.80
  model += 0.20 * p[3] * exp(-((xarr - lamb_s2)/sig_s1)^2)  ;  S IV 1406.04

  model += 0.05 * p[3] * exp(-((xarr - lamb_o5)/p[5])^2)   ;  O IV 1407.39
  model += 0.25 * p[3] * exp(-((xarr - lamb_o6)/sig_o3)^2) ;  O III 1407.70

  return, model
end
