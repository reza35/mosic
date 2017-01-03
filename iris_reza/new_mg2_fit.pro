pro new_mg2_fit, n_iter, px, py, ppy, ee, err_ave, bbc, m, ergd, par, fit

;+
;===============================================================
; procedure : new_mg2_fit.pro
;
; purpose : run the MG II double Gaussian fit with random initialization
;
; n_iter : number of iterations
; px, py, ee, err_ave : input vectors for the fit
; bbc : selection of points in the input vector which are actually fitted  
; m : a scaling factor between [0, 1].
;     It scales the range of the input parameters.  
;
; par, fit :  output variables
;
; Jun 19, 2016 : created
;  
; R.Rezaei @ IAC                         e-mail:  rrezaei@iac.es      
;===============================================================
;   

  f0 = fltarr(3, n_iter)  ; h2v
  f1 = fltarr(3, n_iter)  ; h2r
  h = m * 0.5d
  
  f0[0,*] = ((randomu(seed, n_iter) * m - h) + 1.0) * ergd.i1  
  f0[1,*] = ((randomu(seed, n_iter) * m - h) + 1.0) * ergd.p1 
  f0[2,*] = ((randomu(seed, n_iter) * m - h) + 1.0) * ergd.w1
  
  f1[0,*] = ((randomu(seed, n_iter) * m - h) + 1.0) * ergd.i2
  f1[1,*] = ((randomu(seed, n_iter) * m - h) + 1.0) * ergd.p2
  f1[2,*] = ((randomu(seed, n_iter) * m - h) + 1.0) * ergd.w2
   

  counter = 0
  chisq2 = 1d6
  while (counter lt n_iter-1)and(chisq2 gt 3.) do begin
    fit0 = f0[*, counter]
    fit1 = f1[*, counter]
    range0 = [0.8, 0.3, 0.8]
    range1 = [0.8, 0.3, 0.8]
    dlambda = 1.0d
    pyn  = py
    ergd = my_dgf_mg(px[bbc], pyn[bbc], ee[bbc], fit0, fit1, range0, range1, dlambda[0], bbc, /double)
    fitdg_new = ergd.i1*exp(-((px - ergd.p1)/ergd.w1)^2) + ergd.i2*exp(-((px - ergd.p2)/ergd.w2)^2)
    chisq2_new = (1.0d/(n_elements(bbc) - 6.0d)) * total(((pyn[bbc] - fitdg_new[bbc])/err_ave[bbc])^2)

    if (chisq2_new lt chisq2) then begin
        chisq2 = chisq2_new
        fit = fitdg_new
        par = [ergd.p1, ergd.i1, ergd.w1, ergd.p2, ergd.i2, ergd.w2, chisq2, reform(ergd.sigma)]
    endif   
    counter += 1
  endwhile

end  
