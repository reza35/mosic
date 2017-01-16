function analyze_c2, iline, ca_ilo, ca_ihi, pos_line, disper, plt, fast

;+
;===============================================================
; function : analyze_c2.pro
;  
; purpose : calculates parameters of a "quad Gaussian fit" to both 
; C II line + a penta Gaussian fit to both C II ansd the Ni II
;
; iline: a line profile
; plt  : plot the profile and fit
; pos_line : position of the second C II line  
; fast : it will skip the last (penta-Gaussian) fit and makes things faster
; list of output parameters: single/double/quad/penta Gaussian fit parameters
;
; May 10, 2016 : two C II and Ni II lines
; Jun 04, 2016 : dispersion of 5 mA/px was added
; Nov 21, 2016 : random initialization of the single Gaussian fits
; Dec 19, 2016 : improved random initialization
;
; R.Rezaei @ IAC                         e-mail:  rrezaei@iac.es      
;===============================================================
;-
;
;on_error, 2
common share_ergm, av_fit

if n_elements(plt) eq 0 then plt = 0 else plt = plt
if n_elements(fast) eq 0 then fast = 0 else fast = fast
if n_elements(av_fit) eq 0 then av_fit = [56., 72, 63., 3.]

prof = iline
pprof = gauss_smooth(prof, 1) ; 19

spar1 = fltarr(9)
dpar1 = fltarr(15)
qpar1 = fltarr(21)
ppar1 = fltarr(27)

;------------------------------------------------
;-- line core position
;------------------------------------------------
np = n_elements(prof)
corepos = pos_line
core_pos = corepos
posmax = corepos > 10.

case disper of
   1.2980: first_line = posmax - 88.
   2.5960: first_line = posmax - 44.
   5.1920: first_line = posmax - 22.
endcase
case disper of
   1.2980: dfac = 1.
   2.5960: dfac = 2.
   5.1920: dfac = 4.
endcase

fitsg = fltarr(np)
fitdg = fltarr(np)
fitqg = fltarr(np)
fitpg = fltarr(np)
erg = {spar1:fltarr(9), dpar1:fltarr(15), qpar1:fltarr(21), ppar1:fltarr(27), sprf:fltarr(np), $
       pfit:fltarr(np), dfit:fltarr(np), sfit:fltarr(np)}

py = prof
px = findgen(np)
bcc = where(py gt (-10.1), bppc)
if (bppc le 10) then return, erg

;----------------------------------------------------------
;-- remove continuum slope, if present
;----------------------------------------------------------
p1 = first_line
tmp = fltarr(np)

if (disper lt 1.3) then begin
  tmp[(p1-40):(p1+40)]=1.
  tmp[(posmax-50):(posmax+40)]=1.
endif
if (disper lt 3.)and(disper gt 2.) then begin
   tmp[(p1-20):(p1+20)]=1.
   tmp[(posmax-25):(posmax+20)]=1.
endif
if (disper lt 6.)and(disper gt 5.) then begin
   tmp[(p1-10):(p1+10)]=1.
   tmp[(posmax-12):(posmax+10)]=1.
endif

qk = where(tmp eq 0) ;py lt abs(median(py)) * 1.5)
kont = good_mean(py(qk)) < 40.0 * dfac  ; kont can be affected by CRs
pym = median(py, 9)
res = poly_fit(px(qk), pym(qk), 2, /double)
base_level = res[0] + res[1]*px + res[2] * px^2 ;+ res[3] * px^3
py = py - base_level + kont

rms = stddev(py(qk))

;---------------------------------
; if one line was hit by CRs
;---------------------------------
dx = 12.0/dfac
max1 = max(py[(first_line-dx):(first_line+dx)]) > .1
max2 = max(py[(posmax-dx):(posmax+dx)]) > .1
if ((max1/max2) gt 2.7) then begin
        x = py[0:(first_line+15/dfac)]
        x(where(x gt (max2*2.))) = kont
        py[0:(first_line+15/dfac)] = x
endif   
if ((max2/max1) gt 3.5) then begin
        x = py[(posmax-15/dfac):*]
        x(where(x gt (max1*2.5))) = kont
        py[(posmax-15/dfac):*] = x
endif   

;----------------------------------------------------------
;-- a single Gassian fit  to the second C II line
;----------------------------------------------------------
  if (min(py) lt 0.) then pye = py - min(py) + 1. else pye = py 
  bcc = where(py gt (-10.1), bppc)
  ;---------------------------------------
  ;-- error in input profiles
  ;--------------------------------------------
  ee = 0.5 * sqrt(abs(py)) > 1.0 
  err_ave = ee < (abs(py)*0.9) 
  err_ave = err_ave > 1.0
  ;---------------------------------------
  tmax = (max(smooth(py[posmax-5:posmax+5], 5)) * 0.8) > 1.
  ;--------------------------------------------
  ;-- initial guess
  ;--------------------------------------------
  fit0 = [kont, tmax,  posmax, 10.0/dfac + randomn(seed)*0.1] * 1.0d 
  range0 = [0.1, 0.6, 0.04, 0.7] 
  ppy = median(py, 3)

  tcc = where((px gt (posmax - 20/dfac))and(px lt (posmax + 60/dfac)<(np-1)))
  tcc = bcc

  ergs = my_sgf(px[tcc], ppy[tcc], ee[tcc], fit0, range0, 1.0d, tcc, /double)
  fitsg = rgauss(px, [ergs.b, ergs.i1, ergs.p1,  ergs.w1])
  chisq1 = (1.0d/(n_elements(tcc) - 4.0d)) * total(((ppy[tcc] - fitsg[tcc])/err_ave[tcc])^2)
  is_bad = evaluate_sgf(fit0, ergs)
  
  if (ergs.w1/disper gt 9.)or(is_bad eq 1.)or(n_elements(fitsg) lt np)or(chisq1 gt 10.0) then begin
     py0 = median(py,5)
     py = gauss_smooth(py0, 1.2, /edge_truncate)
     fit0 = [kont, ergs.i1>4., posmax, av_fit[3]+randomn(seed)*0.3]
     range0 = [0.02, 0.5, 0.02, 0.5]
     ergs = my_sgf(px[tcc], ppy[tcc], ee[tcc], fit0, range0, 1.0d, tcc, /double)
     fitsg = rgauss(px, [ergs.b, ergs.i1, ergs.p1, ergs.w1])
     chisq1 = (1.0d/(n_elements(tcc) - 4.0d)) * total(((ppy[tcc] - fitsg[tcc])/err_ave[tcc])^2)
  endif   


  spar1 = [ergs.b, ergs.p1, ergs.i1,  ergs.w1, chisq1, reform(ergs.sigma)]

  ;------------------------------------------------------------------------------
  ;  if the fit is not satisfactory, then use random initialization 
  ;------------------------------------------------------------------------------
  if (chisq1 gt 1.) then begin
        random_sg_fit, 10., px, py, ppy, ee, err_ave, tcc, 0.6, ergs, spar_new, fitsg_new
        if (spar_new[4] lt chisq1) then begin
           chisq1 = spar_new[4]
           fitsg = fitsg_new
           spar1 = spar_new
        endif
  endif   

  if (abs(spar1[1] - pos_line) gt 10.) then spar1[1] = pos_line
  line_core  = spar1[1] > 5. < (np-6) 
 
;--------------------------------------------------------------------------
;--  a double Guassian = a single Gaussian fit to each C II line 
;--------------------------------------------------------------------------
  if (spar1[3] lt 1.0d)or(spar1[3] gt 14.0/dfac) then spar1[3] = 8.0/dfac
  if (spar1[2] gt tmax*1.1) then spar1[2] = tmax
  if (abs(spar1[1] - posmax) gt 9.) then spar1[1] = posmax

  fit0 = [spar1[2]>4., spar1[1], spar1[3] > .2]
  fit1 = [max(fitsg)*0.82d, spar1[1] - 88./dfac] 
  
  range0 = [0.8, 0.2, 0.2] 
  range1 = [0.8,  0.2]

  dlambda = 1.0d
  pyn = ppy - spar1[0]
  ergd = my_cii1(px[bcc], pyn[bcc], ee[bcc], fit0, fit1, range0, range1, dlambda[0], bcc, /double)
  fitdg = double_cii(px, [ergd.i1, ergd.p1, ergd.w1, ergd.i2, ergd.p2])
  chisq2 = (1.0d/(n_elements(bcc) - 5.0d)) * total(((pyn[bcc] - fitdg[bcc])/err_ave[bcc])^2)
  fitdg += spar1[0]
  dpar1=[ergd.p1, ergd.i1, ergd.w1, ergd.p2, ergd.i2, chisq2, reform(ergd.sigma)]

  if (ergd.i1/ergd.i2 gt 1.5)or(ergd.i1/ergd.i2 lt 0.9)or(n_elements(fitdg) lt np)or(chisq2 gt 5.) then begin
     ; we repeat the fit with amplitude of the second line fixed with  a ratio
  
     fit2 = [spar1[2]>1., spar1[1], spar1[3], fit1[1]]
     range2 = [0.2, 0.2, 0.2, 0.2]
     ergd = my_cii1n(px[bcc], pyn[bcc], ee[bcc], fit2, range2, dlambda[0], bcc, /double)

     fit0 = [ergd.i1, ergd.p1, ergd.w1]
     fit1 = [ergd.i1*0.82d, ergd.p2]
     range0 = [0.1, 0.1, 0.1] 
     range1 = [0.1,  0.1]
     ergd = my_cii1(px[bcc], pyn[bcc], ee[bcc], fit0, fit1, range0, range1, dlambda[0], bcc, /double)
     fitdg = double_cii(px, [ergd.i1, ergd.p1, ergd.w1, ergd.i2, ergd.p2])
     chisq2 = (1.0d/(n_elements(bcc) - 5.0d)) * total(((pyn[bcc] - fitdg[bcc])/err_ave[bcc])^2)
     fitdg += spar1[0]
     dpar1=[ergd.p1, ergd.i1, ergd.w1, ergd.p2, ergd.i1 * 0.82, chisq2, reform(ergd.sigma)]

     ;print, chisq2, '-------------------------------------'
  endif

  ;------------------------------------------------------------------------------
  ;  if the fit is not satisfactory, then use random initialization 
  ;------------------------------------------------------------------------------
  if (chisq2 gt 1.)or(chisq2 lt 0.1)or(ergd.i1/ergd.i2 gt 1.5)or(ergd.i1/ergd.i2 lt 0.9) then begin
        random_dg_cii_fit, 20, px, pyn, ee, err_ave, bcc, 0.2, spar1, ergd, dpar_new, fitdg_new
        if (dpar_new[5] lt chisq2) then begin
           chisq2 = dpar_new[5]
           fitdg = fitdg_new + spar1[0]
           dpar1 = dpar_new
        endif
  endif


  
;  print, chisq1, chisq2
;stop
if (max(fitdg)/rms gt 4.) then begin
;---------------------------------------------------
;--  a quad Guassian fit to both C II lines 
;--  = two Gaussians for each C II line
;---------------------------------------------------
ee =  ir_error(py, /fuv, /dark)
  
  v_peak = max(py[(dpar1[0] - 15/dfac):(dpar1[0] -3/dfac)], v_posmax) & v_posmax = (v_posmax + (dpar1[0] - 15/dfac))> (dpar1[0] - disper*3.) 
  r_peak = max(py[(dpar1[0] +3/dfac):(dpar1[0] + 15/dfac)], r_posmax) & r_posmax = (r_posmax + (dpar1[0] +3/dfac)) > (dpar1[0]) < (np -5.)

  b_peak = max(py[((dpar1[3] - 15/dfac)>0):(dpar1[3] -3/dfac)], b_posmax) & b_posmax = (b_posmax + (dpar1[3] - 15/dfac))> (dpar1[3] - disper*3.) 
  g_peak = max(py[(dpar1[3] +3/dfac):(dpar1[3] + 15/dfac)], g_posmax) & g_posmax = (g_posmax + (dpar1[3] +3/dfac)) > (dpar1[3]) < (np -5.)
 
  fit0 = [dpar1[1]*0.8 > 2., v_posmax, dpar1[2], dpar1[4]*0.8>2., b_posmax]
  fit1 = [dpar1[1]*0.8d> 2., r_posmax, dpar1[2], dpar1[4]*0.8>2., g_posmax] 

  range0 = [0.8, 0.1, 0.7, .9, .1] * 1.0d
  range1 = [0.8, 0.1, 0.7, .9, .1] * 1.0d

  dlambda = 1.0d
  pyn = ppy - spar1[0]
  ergq = my_cii2(px[bcc], pyn[bcc], ee[bcc], fit0, fit1, range0, range1, dlambda[0], bcc, /double)
  fitqg = quad_cii(px, [ergq.i1, ergq.p1, ergq.w1, ergq.i2, ergq.p2, ergq.i3, ergq.p3, ergq.w2, ergq.i4, ergq.p4])
  chisq10 = (1.0d/(n_elements(bcc) - 10.0d)) * total(((pyn[bcc] - fitqg[bcc])/err_ave[bcc])^2)
  fitqg += spar1[0]
  ;print, 'b ',  chisq10
  ;print, [ergq.p1, ergq.i1, ergq.w1, ergq.i2, ergq.p2, ergq.p3, ergq.i3, ergq.w2, ergq.i4, ergq.p4]
  ;print, fit0, fit1
  ;---------------------------------------------------------------------------------------
  if (chisq10 gt 20.0) then begin
  fit0 = [dpar1[1] > 2., dpar1[0], dpar1[2], dpar1[4]>2., dpar1[3]]  ; first Gaussian, both lines
  
  ;if (disper lt 2.) then begin
  fit1 = [dpar1[1]*0.3d, dpar1[0] - 3.0/dfac, dpar1[2] * 1.2, dpar1[4]*0.3d, dpar1[3] - 3.0/dfac] 
  ;endif else begin
  ;   fit1 = [dpar1[1]*0.3d, dpar1[0] - 1.5, dpar1[2] * 1.2, dpar1[4]*0.3d, dpar1[3] - 1.5] 
  ;endelse   

  range0 = [0.2, 0.2, 0.2, 0.2, 0.2] * 1.0d
  range1 = [0.5, 0.1, 0.2, 0.5, 0.1] * 1.0d

  dlambda = 1.0d
  pyn = ppy - spar1[0]
  ergq = my_cii2(px[bcc], pyn[bcc], ee[bcc], fit0, fit1, range0, range1, dlambda[0], bcc, /double)
  fitqg_a = quad_cii(px, [ergq.i1, ergq.p1, ergq.w1, ergq.i2, ergq.p2, ergq.i3, ergq.p3, ergq.w2, ergq.i4, ergq.p4])
  chisq10_a = (1.0d/(n_elements(bcc) - 10.0d)) * total(((pyn[bcc] - fitqg[bcc])/err_ave[bcc])^2)
  fitqg_a += spar1[0]
  ergq_a = ergq
  ;---------------------------------------------------------------------------------------
  if (chisq10_a lt chisq10) then begin
     chisq10 = chisq10_a
     fitqg = fitqg_a
     ergq = ergq_a
  endif
  endif
  
  ;if (ergq.i1/ergq.i2 gt 1.6)or(ergq.i1/ergq.i2 lt 0.8)or(n_elements(fitqg) lt np)or(chisq10 gt 50.) then begin ; something is went wrong
  if (n_elements(fitqg) lt np)or(chisq10 gt 50.) then begin ; something went wrong
                                 ; we repeat the fit with amplitude of the second line fixed with  a ratio
     fit2 = [spar1[2]>1., spar1[1], spar1[3], fit1[1]] * 1.0d
     range2 = [0.2, 0.2, 0.2, 0.2] * 1.0d
     ergq = my_cii1n(px[bcc], pyn[bcc], ee[bcc], fit2, range2, dlambda[0], bcc, /double)
     ;print, ergq.i1, ergq.p1, ergq.w1, ergq.p2
     ; now start with this intitial condition
     fit0 = [ergq.i1, ergq.p1, ergq.w1, ergq.i1*0.82d, ergq.p2] * 1.0d
     fit1 = [ergq.i1*0.3, ergq.p1-3.0, ergq.w1*1.1, ergq.i1*0.3, ergq.p2-3.0] * 1.0d
     range0 = [0.1, 0.1, 0.1, 0.1, 0.1] * 1.0d
     range1 = [0.1, 0.1, 0.1, 0.1, 0.1] * 1.0d
     ergq = my_cii2(px[bcc], pyn[bcc], ee[bcc], fit0, fit1, range0, range1, dlambda[0], bcc, /double)
     fitqg = quad_cii(px, [ergq.i1, ergq.p1, ergq.w1, ergq.i2, ergq.p2, ergq.i3, ergq.p3, ergq.w2, ergq.i4, ergq.p4])
     chisq10 = (1.0d/(n_elements(bcc) - 10.0d)) * total(((pyn[bcc] - fitqg[bcc])/err_ave[bcc])^2)
     fitqg += spar1[0]
     ;print, '-------------------------------------'
      ;print, 'c ',  chisq10
  endif   

  qpar1=[ergq.p1, ergq.i1, ergq.w1, ergq.p2, ergq.i2, $
         ergq.p3, ergq.i3, ergq.w2, ergq.p4, ergq.i4, $
         chisq10, reform(ergq.sigma)]

;----------------------------------------------------------------------
;--  a penta Guassian fit to both C II lines and the Ni II 133.520
;----------------------------------------------------------------------
;###########################################################################################
if (fast eq 0) then begin  
  fit0 = [ergq.i1 > 2., ergq.p1, ergq.w1>.2, ergq.i2>2., ergq.p2]  ; first Gaussian, both lines
  fit0 = [fit0, ergq.i3>2., ergq.p3, ergq.w2>.2, ergq.i4>2., ergq.p4] ; 2nd Gaussian, both line

  ;if (disper lt 2.) then begin
  fit1 = [max(fitsg)*0.1d > 2., spar1[1] - 38./dfac, (dpar1[2] * 0.4)>.2]   ; Ni II parameters
  ;endif else begin
  ;   fit1 = [max(fitsg)*0.1d > 2., spar1[1] - 19., (dpar1[2] * 0.4)>.2]
  ;endelse   
 
  range0 = fltarr(10)+ 0.1 
  range1 = [0.9, 0.05, .3]  ; note the restricted range for the line position
  ee =  ir_error(py, /fuv, /dark)  
  dlambda = 1.0d
  pyn = ppy - spar1[0]
  ;ee[((spar1[1] - 38/dfac)>0):(spar1[1] - 20/dfac)]=100.
  ;print, fit0, fit1
  ergp = my_cii3(px[bcc], pyn[bcc], ee[bcc], fit0, fit1, range0, range1, dlambda[0], bcc, /double)
  fitpg = penta_cii(px, [ergp.i1, ergp.p1, ergp.w1, ergp.i2, ergp.p2, ergp.i3, ergp.p3, ergp.w2, ergp.i4, ergp.p4, ergp.i5, ergp.p5, ergp.w3])
  chisq13 = (1.0d/(n_elements(bcc) - 13.0d)) * total(((pyn[bcc] - fitpg[bcc])/err_ave[bcc])^2)
  fitpg += spar1[0]

  if (ergp.w3 lt 1.)or(ergp.w3 gt (ergp.w1 * .9))or(ergp.i5 lt 0.3)or(chisq13 gt 50.) then begin  ; in case of bad Ni II fit
     fit1[2] = ergp.w1 /2.
     range1 = [0.9, 0.03, .5]   ; note the restricted range for the line position
  
     dlambda = 1.0d
     pyn = ppy - spar1[0]
     ergp = my_cii3(px[bcc], pyn[bcc], ee[bcc], fit0, fit1, range0, range1, dlambda[0], bcc, /double)
     fitpg = penta_cii(px, [ergp.i1, ergp.p1, ergp.w1, ergp.i2, ergp.p2, ergp.i3, ergp.p3, ergp.w2, ergp.i4, ergp.p4, ergp.i5, ergp.p5, ergp.w3])
     chisq13 = (1.0d/(n_elements(bcc) - 13.0d)) * total(((pyn[bcc] - fitpg[bcc])/err_ave[bcc])^2)
     fitpg += spar1[0]
  endif
  ppar1 = [ergp.p1, ergp.i1, ergp.w1, ergp.p2, ergp.i2, $
         ergp.p3, ergp.i3, ergp.w2, ergp.p4, ergp.i4, $
         ergp.p5, ergp.i5, ergp.w3, $  ; Ni II
         chisq13, reform(ergp.sigma)]
  
endif  
;###########################################################################################

endif

;-------------------------------------
;-- plot lines, positions of peaks
;-------------------------------------
if (plt eq 1) then begin
  print, '------------------------------------------------------'
  print, 'single Gaussian:', spar1[0:4]
  print, 'double Gaussain: ', dpar1[0:5]
  print, 'quad Gaussain: ', qpar1[0:10]
  print, 'penta Gaussain: ', ppar1[0:13]
  print, '------------------------------------------------------'
  plot, px, py, psym=-1, /xsty;, yr=[-2,11]
  loadct, 40, /silent
  
  oplot, px, fitsg, thick=2, color=110
  oplot, px, fitdg, thick=2, color=150
  oplot, px, fitqg, thick=2, color=50
  if (fast eq 0) then oplot, px, fitpg,color=245, linestyle=2
  loadct, 0, /silent
  xyouts, 0.15, .7, num2string(round(max(fitdg)/rms)), chars=2,/normal
  ans =''  &  read, ans, prompt='Press enter to continue ...'
endif

;----------------------------
;--  Set output parameters
;----------------------------
erg.spar1 = spar1         ; single Gaussian fit parameters
erg.dpar1 = dpar1         ; double Gaussian fit parameters
erg.qpar1 = qpar1         ; quad Gaussian fit parameters
erg.ppar1 = ppar1         ; penta Gaussian fit parameters

erg.sprf  = py            ; line profile
erg.sfit  = fitsg         ; single Gaussian fit to each C II line profile
erg.dfit  = fitdg         ; double Gaussian fit to each C II line profile
erg.pfit  = fitpg         ; penta Gaussian fit to each C II line profile + Ni II

return,erg
end
