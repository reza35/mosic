function pentagauss_lvg, xarr, p

;+
;===============================================================
; function : pentagauss_lvg.pro
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
; sigma of small lines are computed based on sigma of stronger lines
;
; Jun 10, 2016:  original
;
; R.Rezaei @ IAC                         e-mail:  rrezaei@iac.es      
;===============================================================
;-   


common share_disp, z, u
common shared_pars, mp
;print, '-----', p[0]
  df = z / 1.298

  lamb_s  = mp[1] + p[0] + 158.760d /df ; in pixel unit  >> weak  S IV line at 1404.79 A
  lamb_o3 = mp[1] + p[0] + 160.56d /df ; in pixel unit  >> weak   O IV line at 1404.82 A
  lamb_o2 = mp[1] + p[0] - 235.60d /df ; in pixel unit  >> weak   O IV line at 1399.77 A
  lamb_o1 = mp[1] + p[0] - 126.32d /df ; in pixel unit  >> strong O IV line at 1401.16 A

  lamb_s1 = mp[1] + p[0] - 102.32d /df ; in pixel unit  >> weak S I line at 1401.51 A
  lamb_fe2 = mp[1] + p[0] + 208.0d /df ; in pixel unit  >> weak  Fe II line at 1405.61 A
  lamb_o4 = mp[1] + p[0] + 214.0d /df ; in pixel unit  >> weak  O III line at 1405.80 A
  lamb_s2 = mp[1] + p[0] + 220.0d /df ; in pixel unit  >> weak  S IV line at 1406.04 A
  lamb_fe1 = mp[1] + p[0] - 222.60d /df ; in pixel unit  >> weak  Fe II line at 1399.97 A
  lamb_ni = mp[1] + p[0] - 298.4d /df ; in pixel unit  >> weak  Ni II line at 1399.03 A

  
  sig_o1 = mp[2] * 1.12d      ; sigma of weak O IV line at 1401.16 A
  sig_s1 = mp[2] * 0.3d      ; sigma of weak S I line at 1401.51 A
  sig_s2 = mp[2] * 0.35d      ; sigma of weak S IV line at 1406.04 A
  sig_fe = mp[2] * 0.21d      ; sigma of weak Fe II line at 1405.61 A
  sig_ni = mp[2] * 0.31d      ; sigma of weak Ni II line at 1399.03 A

  model = mp[0]*exp(-((xarr - mp[1])/mp[2])^2) ; single Guassian fit to Si IV 1403

  model += mp[3]*exp(-((xarr - lamb_o1)/sig_o1)^2)  ; O IV 1401.16 

  model += mp[4]*exp(-((xarr - lamb_o2)/sig_o1)^2)  ; O IV 1399.97

  model += mp[5]*exp(-((xarr - lamb_s)/mp[2])^2)     ; S IV 1404.79

  model += mp[6]*exp(-((xarr - lamb_o3)/sig_o1)^2)  ; O IV 1404.82

;  model += 0.44 * mp[3]*exp(-((xarr - lamb_s1)/sig_s1)^2)  ; S I 1401.26 
;  model += 0.05 * mp[3]*exp(-((xarr - lamb_fe2)/sig_fe)^2) ;  Fe II 1405.61
;  model += 0.20 * mp[3]*exp(-((xarr - lamb_o4)/sig_fe)^2) ;  O III 1405.80
;  model += 0.37 * mp[3]*exp(-((xarr - lamb_s2)/sig_s1)^2) ;  S IV 1406.04
;  model += 0.12 * mp[3]*exp(-((xarr - lamb_fe1)/sig_fe)^2) ;  Fe II 1399.97
;  model += 0.29 * mp[3]*exp(-((xarr - lamb_ni)/sig_ni)^2) ;  Ni II 1399.03
 
  return, model
end
