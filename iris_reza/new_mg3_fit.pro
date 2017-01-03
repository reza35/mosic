pro new_mg3_fit, n_iter, px, py, ppy, ee, err_ave, bbc, m, ergt, par, fit

;+
;===============================================================
; procedure : new_mg3_fit.pro
;
; purpose : run the MG II triple Gaussian fit with random initialization
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
;   

  
  f0 = fltarr(3, n_iter)  ; h3
  f1 = fltarr(6, n_iter)  ; h2v and h2r
  h = m * 0.5d
  
  f0[0,*] = ((randomu(seed, n_iter) * m - h) + 1.0) * ergt.i1
  f0[1,*] = ((randomu(seed, n_iter) * m - h) + 1.0) * ergt.p1
  f0[2,*] = ((randomu(seed, n_iter) * m - h) + 1.0) * ergt.w1

  f1[0,*] = ((randomu(seed, n_iter) * m - h) + 1.0) * ergt.i2
  f1[1,*] = ((randomu(seed, n_iter) * m - h) + 1.0) * ergt.p2
  f1[2,*] = ((randomu(seed, n_iter) * m - h) + 1.0) * ergt.w2
  f1[3,*] = ((randomu(seed, n_iter) * m - h) + 1.0) * ergt.i3
  f1[4,*] = ((randomu(seed, n_iter) * m - h) + 1.0) * ergt.p3
  f1[5,*] = ((randomu(seed, n_iter) * m - h) + 1.0) * ergt.w3
   
  counter = 0
  chisq4 = 1d6
  while (counter lt n_iter-1)and(chisq4 gt 1.0) do begin
    fit0 = f0[*, counter]
    fit1 = f1[*, counter]
    range0 = [0.8, 0.3, 0.8]
    range1 = [0.8, 0.3, 0.4, 0.8, 0.2, 0.4]
    dlambda = 1.0d
    pyn  = py
    ergt = my_tgf_mg(px[bbc], pyn[bbc], ee[bbc], fit0, fit1, range0, range1, dlambda[0], bbc, /double)
    fittg_new = tgf_mg(px, [ergt.i1, ergt.p1, ergt.w1, ergt.i2, ergt.p2, ergt.w2, ergt.i3, ergt.p3, ergt.w3])
    chisq4_new = (1.0d/(n_elements(bbc) - 9.0d)) * total(((pyn[bbc] - fittg_new[bbc])/err_ave[bbc])^2)

    if (chisq4_new lt chisq4) then begin
        chisq4 = chisq4_new
        fit = fittg_new
        par = [ergt.p1, ergt.i1, ergt.w1, ergt.p2, ergt.i2, ergt.w2, ergt.p3, ergt.i3, ergt.w3, chisq4, reform(ergt.sigma)]
    endif   
    counter += 1
  endwhile

end  
