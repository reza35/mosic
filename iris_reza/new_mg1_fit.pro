pro new_mg1_fit, n_iter, px, py, ppy, ee, err_ave, bbc, m, ergs, par, fit

;+
;===============================================================
; procedure : new_mg1_fit.pro
;
; purpose : run the MG II single Gaussian fit with random initialization
;
; n_iter : number of iterations
; px, py, ee, err_ave : input vectors for the fit
; bbc : selection of points in the input vector which are actually fitted  
; m : a variable between [0, 1].
;     m scales the range of the input parameters.  
;
; par, fit :  output variables
;
; Jun 19, 2016 : created
;  
; R.Rezaei @ IAC                         e-mail:  rrezaei@iac.es      
;===============================================================
;-  

  
  f0 = fltarr(4, n_iter) 
  h = m * 0.5d
  
  f0[0,*] = ((randomu(seed, n_iter) * m - h) + 1.0) * ergs.b  
  f0[1,*] = ((randomu(seed, n_iter) * m - h) + 1.0) * ergs.i1  
  f0[2,*] = ((randomu(seed, n_iter) * m - h) + 1.0) * ergs.p1 
  f0[3,*] = ((randomu(seed, n_iter) * m - h) + 1.0) * ergs.w1
     

  counter = 0
  chisq1 = 1d6
  while (counter lt n_iter-1)and(chisq1 gt 5.) do begin
    fit0 = f0[*, counter]
    range0 = [0.2, 0.8, 0.3, 0.8]
    dlambda = 1.0d
    pyn  = py
    ergs_new = my_sgf(px[bbc], pyn[bbc], ee[bbc], fit0, range0, dlambda[0], bbc, /double)
    fitsg_new = rgauss(px, [ergs_new.b, ergs_new.i1, ergs_new.p1,  ergs_new.w1])
    chisq1_new = (1.0d/(n_elements(bbc) - 4.0d)) * total(((pyn[bbc] - fitsg_new[bbc])/err_ave[bbc])^2)

    if (chisq1_new lt chisq1) then begin
        chisq1 = chisq1_new
        fit = fitsg_new
        par = [ergs.b, ergs.p1, ergs.i1, ergs.w1, chisq1, reform(ergs.sigma)]
    endif   
    counter += 1
  endwhile

end  
