pro random_dg_cii_fit, n_iter, px, py, ee, err_ave, bbc, m, sp, ergd, par, fit

;+
;===============================================================
; procedure : random_dg_cii_fit.pro
;
; purpose : run a double Gaussian fit with random initialization
;
; n_iter : number of iterations
; px, py, ee, err_ave : input vectors for the fit
; bbc : selection of points in the input vector which are actually fitted  
; m : a variable between [0, 1].
;     m scales the range of the input parameters.  
;
; par, fit :  output variables
;
; Dec 19, 2016 : created
;  
; R.Rezaei @ IAC                         e-mail:  rrezaei@iac.es      
;===============================================================
;-  

  
  f0 = fltarr(3, n_iter) 
  f1 = fltarr(2, n_iter) 
  h = m * 0.5d
  
  f0[0,*] = ((randomu(seed, n_iter) * m - h) + 1.0) * sp[2]  
  f0[1,*] = ((randomu(seed, n_iter) * m - h) + 1.0) * sp[1] 
  f0[2,*] = ((randomu(seed, n_iter) * m - h) + 1.0) * sp[3]
     
  f1[0,*] = ((randomu(seed, n_iter) * m - h) + 1.0) * sp[2] * 0.82  
  f1[1,*] = ((randomu(seed, n_iter) * m - h) + 1.0) * ergd.p2 

  counter = 0
  chisq2 = 1d6
  while (counter lt n_iter-1)and(chisq2 gt 1.) do begin
    fit0 = f0[*, counter]
    fit1 = f1[*, counter]
    
    range0 = [0.1, 0.1, 0.1]
    range1 = [0.1, 0.1]
    dlambda = 1.0d
    
    pyn  = py
    ergd = my_cii1(px[bbc], pyn[bbc], ee[bbc], fit0, fit1, range0, range1, dlambda[0], bbc, /double)
    fitdg_new = double_cii(px, [ergd.i1, ergd.p1, ergd.w1, ergd.i2, ergd.p2])
    chisq2_new = (1.0d/(n_elements(bbc) - 5.0d)) * total(((pyn[bbc] - fitdg_new[bbc])/err_ave[bbc])^2)

    if (chisq2_new lt chisq2) then begin
        chisq2 = chisq2_new
        fit = fitdg_new
        par = [ergd.p1, ergd.i1, ergd.w1, ergd.p2, ergd.i2, chisq2, reform(ergd.sigma)]
    endif   
    counter += 1
  endwhile

end  
