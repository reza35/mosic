function analyze_occ, iline, ca_ilo, ca_ihi, pos_line, disper, plt
;+
;===============================================================
; function : analyze_occ
; purpose : calculates parameters of a "Triple Gaussian fit" to
; C I 1354.288, O I 1355.598, and C I 1355.84    (blended)
;
; iline: a line profile
;
;
; - We use first a Single Gaussian fit. Then use its results as
;   initial guess for teh triple Gaussian fit
;
; list of output parameters:
;
; spar[0] :   baseline = continuum -------gaussian  fit
; spar[1] :   line-core position ------- gaussian fit
; spar[2] :   line-core amplitude ------- gaussian fit
; spar[3] :   line-core width ------- gaussian fit  > multiply by
;                                     2. sqrt(2. * alog(2.)) =  FWHM
; spar[4] :   chi^2 -------gaussian  fit
; spar[5]:   error in baseline ------- gaussian fit
; spar[6] :   error in line position -------gaussian  fit
; spar[7] :   error in amplitude ------- gaussian fit
; spar[8] :   error in line width ------- gaussian fit
;
;
; qpar[0] :   baseline = continuum -------gaussian  fit
; qpar[1] :   line-core position ------- gaussian fit O I 1355.598
; qpar[2] :   line-core amplitude ------- gaussian fit O I1355.598
; qpar[3] :   line-core width ------- gaussian fit  O I  > multiply by
;                                     2. sqrt(2. * alog(2.)) =  FWHM
; qpar[4] :   line-core amplitude ------- gaussian fit Cl I 1355.84
; qpar[5] :   line-core width ------- gaussian fit  Cl I 1355.84
; qpar[6] :   line-core amplitude ------- gaussian fit Cl I 1354.288
; qpar[7] :   line-core width ------- gaussian fit Cl I 1354.288
; qpar[8] :   chi^2 -------gaussian  fit
; qpar[9:14] :   1-sigma error of the fit

; band[0:2]   : H3, H2v, H2r, each 9 pixel wide
;
;
; Dec 08, 2014 : output for single Gaussian fit
; Dec 09, 2014 : program for double Gaussian fit using MPFIT
;                single Gaussian fit was completely re-written
;
; Jun 17, 2015 : applied first to Cl I + O I + Cl I data at 135 nm
; May 12, 2016 : chi-square statistic improved
; May 24, 2016 : new loop to double check failed fits / CR profiles
; Oct 03, 2016 : minor bug fix for spectral range
; Nov 21, 2016 : random initialization of the single Gaussian fits
;
; R.Rezaei @ IAC                         e-mail:  rrezaei@iac.es      
;===============================================================
;-

;on_error, 2
common share_ergm, av_fit

left =  0
right =  n_elements(iline)-1
if n_elements(winindex) eq 0 then windex  = 20 else windex = winindex
if n_elements(av_fit) eq 0 then av_fit = [1., 4., pos_line,  2.]

prof = iline

band1 = 1.0d
spar1 = fltarr(9)
dpar1 = fltarr(15)
qpar1 = fltarr(14)
;------------------------------------------------
;-- line core position
;------------------------------------------------
np = n_elements(prof[left:right])
corepos = pos_line - left
core_pos = corepos
posmax = pos_line
case disper of
   1.2980: dfac = 1.
   2.5960: dfac = 2.
   5.1920: dfac = 4.
endcase

fitsg = fltarr(np)
fitqg = fltarr(np)
chisq4 = 0.
erg = {spar1:fltarr(9), dpar1:fltarr(15), qpar1:fltarr(14), sprf:fltarr(np), $
       sfit:fltarr(np), dfit:fltarr(np), qfit:fltarr(np)}
;----------------------------------------------------------
;-- define the continuum wavelength range
;----------------------------------------------------------
py = reform(prof[left:right])
px = findgen(np)
bbc = where(py gt (-5.0), bppc) ;& if (bppc lt np) then py(where(py lt -10.)) = 0. 
if (bppc le 0) then return, erg

qk = [px[0:(posmax-120/dfac)>0], px[(posmax-70/dfac):(posmax-20.0/dfac)], px[(posmax+36/dfac)<(np-1):*]]

kont = good_mean(py(qk)) < 40.0 * dfac  ; kont will be affected by CRs
noise = stddev(py(qk))
o_max = 190.                      ; maximum amplitude allowed for the O I line (first guess)
;----------------------------------------------------------
;-- remove CRs in the continuum windows
;----------------------------------------------------------
u = py(qk)
q = where(u gt 20., count)
if (count gt 1) then u(q)= kont
py(qk) = u

;----------------------------------------------------------
;-- remove continuum slope, if present
;----------------------------------------------------------
res = poly_fit(px(qk), py(qk), 2, /double)
base_level = res[0] + res[1]*px + res[2] * px^2
py = py - base_level + kont
rms = stddev(py(qk))

sd = where(py le (-10.), count)
if (count gt 0) then py(sd) = kont

pye = py
bbc = where(py gt (-5.0), bppc) 
if (bppc le 20) then return, erg

ee = 0.5 * sqrt(abs(py)) > 1.0 ;ir_error(py, /fuv, /dark) < py ;0.5 * sqrt(abs(py)) > 0.5
err_ave = ee < (abs(py)*0.9)   ;make sure that the error value is larger than 0 and smaller than the data value
err_ave = err_ave > 1.0

sd = where(py gt 500., count)
if (count gt 0) then begin
   ee(sd) = 100.
   err_ave(sd) = 100.
endif
tmax = (max(py[posmax-10:(posmax+10)<(np-1)]) * 1.0d) > 1. < o_max

pyn = median(py,5)
c1 = max(pyn[(((posmax-120)>0)/dfac):(posmax-70/dfac)])
c2 = max(pyn[(posmax+10/dfac):(posmax+40/dfac)<(np-1)])
;--------------------------------------------
;-- initial guess, single Gaussain fit
;--------------------------------------------
fit0 = [kont, tmax,  posmax , 4.5/dfac + randomn(seed)*.1] * 1.0d

ppy = py

range0=[0.1, 0.2, 0.04, 0.4] ; 0.05 0.25 ;0.9 
dlambda = 1.0d
ee[0:((av_fit[2]-20/dfac)>0)]= 1.0d2
ee[((av_fit[2]+20/dfac)<(np-1)):*]= 1.0d2

ergs = my_sgf(px[bbc], ppy[bbc], ee[bbc], fit0, range0, dlambda[0], bbc, /double)
fitsg = rgauss(px, [ergs.b, ergs.i1, ergs.p1,  ergs.w1])
chisq1 = (1.0d/(n_elements(ppy[bbc]) - 4.0d)) * total(((ppy[bbc] - fitsg[bbc])/err_ave[bbc])^2)

is_bad = evaluate_sgf(fit0, ergs)
if (abs(ergs.p1 - posmax) gt 8.) then is_bad = 1.  ; we know that O I cannot have large Doppler shift

;------------------
if (ergs.i1 gt 1d2)or(ergs.w1 gt 9.0/dfac)or(is_bad eq 1.)or(n_elements(fitsg) lt np)or(chisq1 gt 10.0) then begin ; the first check
;print, 'aaaaaaaaaaaaaaaaa'
   fit0 = [kont, tmax, posmax, av_fit[3]+randomn(seed)*0.1]
     range0=[.02,  0.05, 0.04,  0.2] 
     ergs = my_sgf(px[bbc], ppy[bbc], ee[bbc], fit0, range0, dlambda[0], bbc, /double)
     fitsg = rgauss(px, [ergs.b, ergs.i1, ergs.p1,  ergs.w1])
     chisq1 = (1.0d/(n_elements(ppy[bbc]) - 4.0d)) * total(((ppy[bbc] - fitsg[bbc])/err_ave[bbc])^2)
     is_bad = evaluate_sgf(fit0, ergs)
endif   
;------------------
if (ergs.i1 gt 1d2)or(ergs.w1 gt 9.0/dfac)or(is_bad eq 1.)or(n_elements(fitsg) lt np)or(chisq1 gt 10.) then begin ; the second check
   ;print, '--------'
     fit0[3] = randomu(seed)*4.3/dfac + 1.
     glx = max(ppy[(posmax-8):(posmax+8)<(np-1)], p_ii) & p_ii += (posmax-8)
     fit0[1:2] = [glx<100., p_ii+1.]
     range0 = [0.4, 0.2, 0.04, 0.2] 
     ergs = my_sgf(px[bbc], ppy[bbc], ee[bbc], fit0, range0, dlambda[0], bbc, /double)
     fitsg = rgauss(px, [ergs.b, ergs.i1, ergs.p1,  ergs.w1])
     chisq1 = (1.0d/(n_elements(ppy[bbc]) - 4.0d)) * total(((ppy[bbc] - fitsg[bbc])/err_ave[bbc])^2)
     is_bad = evaluate_sgf(fit0, ergs)
endif  
;------------------
if (is_bad eq 1.)or(n_elements(fitsg) lt np)or(chisq1 gt 50.) then begin ; the third check
   ;print, ',,,,,,,,,,,,,,,,,,,,,,,', chisq1
     if (max(ppy) gt 100.) then ppy(where(ppy gt 100)) = 0.
     ppy = median(ppy, 7)
     fit0[3] = randomu(seed)*5.6/dfac + 1.
     glx = max(ppy[(posmax-8):(posmax+8)<(np-1)], p_ii) & p_ii += (posmax-8)
     fit0[1:2] = [glx < 100., p_ii+1.]
     range0 = [0.1, 0.1, 0.04, 0.1] 
     ergs = my_sgf(px[bbc], ppy[bbc], ee[bbc], fit0, range0, dlambda[0], bbc, /double)
     fitsg = rgauss(px, [ergs.b, ergs.i1, ergs.p1,  ergs.w1])
     chisq1 = (1.0d/(n_elements(ppy[bbc]) - 4.0d)) * total(((ppy[bbc] - fitsg[bbc])/err_ave[bbc])^2)
endif  
;------------------
if (ergs.i1 gt o_max)and(ergs.w1 lt 3.) then begin ; the last check
     fit0[3] = (randomu(seed)*5.9/dfac + 1.) > av_fit[3]
     fit0[1:2] = [20.0 + randomn(seed)*4., posmax]
     range0 = [0.01, 0.05, 0.02, 0.3] 
     ergs = my_sgf(px[bbc], ppy[bbc], ee[bbc], fit0, range0, dlambda[0], bbc, /double)
     fitsg = rgauss(px, [ergs.b, ergs.i1, ergs.p1,  ergs.w1])
     chisq1 = (1.0d/(n_elements(ppy[bbc]) - 4.0d)) * total(((ppy[bbc] - fitsg[bbc])/err_ave[bbc])^2)
endif  
;------------------
spar1=[ergs.b, ergs.p1, ergs.i1,  ergs.w1, chisq1, reform(ergs.sigma)]

;------------------------------------------------------------------------------
;  if the fit is not satisfactory, then use random initialization 
;------------------------------------------------------------------------------
if (chisq1 gt 1.) then begin
    random_sg_fit, 30, px, py, ppy, ee, err_ave, bbc, 0.6, ergs, spar_new, fitsg_new
    if (spar_new[4] lt chisq1) then begin
         chisq1 = spar_new[4]
         fitsg = fitsg_new
         spar1 = spar_new
    endif
endif   

if (abs(ergs.p1 - pos_line) gt 10.) then ergs.p1 = pos_line
line_core  = ergs.p1 > 5. < (np-6) 

sd = where(py gt 500., count)
if (count gt 0) then begin
   ee(sd) = 100.
   err_ave(sd) = 100.
endif
err_ave = err_ave > 0.5

if (max(fitsg)/rms gt 2.) then begin ;   6
;-------------------------------------------------------------------
;--  a Triple Guassian fit to the O I line profile + two C I lines 
;-------------------------------------------------------------------
  o_max = 80.0
   
  if (ergs.w1 lt 1.0d)or(ergs.w1 gt 8.0d) then ergs.w1 = 4.0d/dfac + randomn(seed)*.2
  if (ergs.i1 gt tmax*1.1)or(ergs.i1 lt tmax*0.9) then ergs.i1 = tmax
  if (abs(ergs.p1 - posmax) gt 8.) then ergs.p1 = posmax

  fit0 = [ergs.i1>1. < o_max, ergs.p1, ergs.w1 > .2]
  fit1 = [max(fitsg)*0.5d> 1. < 19., ergs.w1*1.15d + randomn(seed)*.05, max(fitsg)*0.25d> 0.5 <17.]

  range0=[0.1, 0.05, 0.05] 
  range1=[0.95, 0.5, 0.95]

  dlambda = 1.0d
  pyn = py - ergs.b
  ergq = my_occ(px[bbc], pyn[bbc], ee[bbc], fit0, fit1, range0, range1, dlambda[0], bbc, /double)
  fitqg = triple_occ(px, [ergq.i1, ergq.p1, ergq.w1, ergq.i2, ergq.w2, ergq.i3])
  chisq4 = (1.0d/(n_elements(bbc) - 6.0d)) * total(((pyn[bbc] - fitqg[bbc])/err_ave[bbc])^2)
  fitqg += ergs.b
  qpar1=[ergq.p1, ergq.i1, ergq.w1, ergq.i2, ergq.w2, ergq.i3, chisq4, reform(ergq.sigma)]
  
  ;---------------------------------------
  ; try agin in case of  a failed fit
  ;---------------------------------------
  if (ergq.i1 gt 7d1)or(ergq.w1 gt 8.0/dfac)or(ergq.i3 gt 7.)or(ergq.w1 lt 0.2)or(n_elements(fitqg) lt np)or(chisq4 gt 2.0d1) then begin ; the first check
     py0 = median(py, 5)
     py = gauss_smooth(py0, 1.2, /edge_truncate)
     pyn = py - ergs.b
     
     fit0 = [ergs.i1>2. < o_max, ergs.p1, 5.0/dfac+randomn(seed)*0.3]
     fit1 = [max(fitsg)*0.8d > 1.< 90., 7.d/dfac+randomn(seed)*0.3, max(fitsg)*0.2d > 0.5 < 60.]

     if (c1 gt 15)and(c2 gt 50.) then begin
        o_max = 180.0
        fit0 = [ergs.i1>30. < o_max, ergs.p1, 5.0/dfac+randomn(seed)*0.3]
        fit1 = [max(fitsg)*0.8d > 1.< 180., 7.d/dfac+randomn(seed)*0.3, max(fitsg)*0.2d > 0.5 < 120.]
     endif   
     if (c1 gt 100)and(c2 gt 200.) then begin
        o_max = 540.0
        fit0 = [ergs.i1>100. < o_max, ergs.p1, 5.0/dfac+randomn(seed)*0.3]
        fit1 = [max(fitsg)*0.8d > 1.< 540., 7.d/dfac+randomn(seed)*0.3, max(fitsg)*0.2d > 0.5 < 360.]
     endif   
     
     range0 = [0.2, 0.2, 0.3] 
     range1 = [0.9,  0.3, 0.8]
     
     ergq = my_occ(px[bbc], pyn[bbc], ee[bbc], fit0, fit1, range0, range1, dlambda[0], bbc, /double)
     fitqg = triple_occ(px, [ergq.i1, ergq.p1, ergq.w1, ergq.i2, ergq.w2, ergq.i3])
     chisq4 = (1.0d/(n_elements(bbc) - 6.0d)) * total(((pyn[bbc] - fitqg[bbc])/err_ave[bbc])^2)
     fitqg += ergs.b
     qpar1=[ergq.p1, ergq.i1, ergq.w1, ergq.i2, ergq.w2, ergq.i3, chisq4, reform(ergq.sigma)]
  endif
  ;---------------------------------------
  ; try agin in case of  a failed fit
  ;---------------------------------------
  if (chisq4 gt chisq1*2.) then begin 
     py0 = median(py, 5)
     py = gauss_smooth(py0, 1.2, /edge_truncate)
     pyn = py - ergs.b
     
     fit0 = [ergs.i1>2. < o_max, ergs.p1, 5.0/dfac + randomn(seed)*1.5]
     fit1 = [max(fitsg)*0.1d > 1.< 60., 6.0/dfac + randomn(seed)*0.8, max(fitsg)*0.05d > 0.5 < 40.]
     
     range0 = [0.01, 0.05, 0.05] 
     range1 = [0.3,  0.3, 0.3]
     
     ergq = my_occ(px[bbc], pyn[bbc], ee[bbc], fit0, fit1, range0, range1, dlambda[0], bbc, /double)
     fitqg = triple_occ(px, [ergq.i1, ergq.p1, ergq.w1, ergq.i2, ergq.w2, ergq.i3])
     chisq4 = (1.0d/(n_elements(bbc) - 6.0d)) * total(((pyn[bbc] - fitqg[bbc])/err_ave[bbc])^2)
     fitqg += ergs.b
     qpar1=[ergq.p1, ergq.i1, ergq.w1, ergq.i2, ergq.w2, ergq.i3, chisq4, reform(ergq.sigma)]
     ;oplot, px, fitqg - ergs.b, color=250 & loadct,0,/silent
     ;print, '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  endif
  endif
;print, chisq1, chisq4
;-------------------------------------
;-- plot lines, positions of peaks
;-------------------------------------
if (plt eq 1) then begin
  print, '------------------------------------------------------'
  print, 'single Gaussain: ', spar1[0:4]
  print, 'tripe Gaussain: ', qpar1[0:6]
  print, '------------------------------------------------------'
  plot, px, ppy, /xsty
  loadct, 40, /silent
  oplot, px, fitsg, thick=2, color=50
  if (qpar1[0] ne 0.) then oplot, px, fitqg,color=245, linestyle=2
  loadct, 0, /silent
  xyouts, 0.2, .7, num2string(round(max(fitsg)/rms)), chars=2,/normal
  ans =''  &  read, ans, prompt='Press enter to continue ...'
endif

;----------------------------
;--  Set output parameters
;----------------------------

erg.spar1 = spar1         ; single Gaussian fit parameters
erg.qpar1 = qpar1         ; triple Gaussian fit parameters

erg.sprf  = py            ; line profile
erg.sfit  = fitsg         ; single Gaussian fit to the line profile
erg.qfit  = fitqg         ; triple Gaussian fit to the line profile

return,erg
end
