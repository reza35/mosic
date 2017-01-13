function analyze_si, iline, ca_ilo, ca_ihi, pos_line, disper, plt
;+
;===============================================================
; function : analyze_si
;  
; purpose : calculates parameters of the Si IV 1403, S IV 1404.79
;           O IV 1399.97, O IV 1401.16, O IV 1404.82
;
; background : Peter 2001, A&A 360, 761
;
; iline: a line profile
; ca_ilo/ca_ihi : the lower/upper range
; pos_line : line center position
; disper : spectral dispersion in mA/px
; plt=1 plots the line and the fit
;  
;   The workflow:
;   ------------------
;   At first we use a single Gaussian fit. As you can read in the paper
;   mentioned above and see in the data, the Si IV profiles are often
;   more complicated than a single Gaussian, which motivates a
;   multi-Gaussian fitting approach. In complex situations, emission lines
;   show self absorption, so will be like Ca II lines. That again deserves a
;   double or triple Gaussian fit which will be performed right after
;   the single Gaussian fit. Beside the Si IV line, there are three O IV lines
;   in the nearby wavelength band along with some other weaker lines.
;   Therefore, a penta-Gaussian fit is performed after the single and double Gaussian.
;   This fit uses a fixed ratio for the widths of the Si IV and O IV lines,
;   and the wavelength positions are fixed w.r.t. each other.
;   In other words, we only fit the Si IV 1403 with three free parameters
;   and for others only fit their amplitudes. This fit has in total 7 free
;   parameters, equal to number of free parameters in the double Gaussian
;   fit of Si IV line alone (the single Gaussian has four free parameters).
;   Please note that the range of parameters are very restricted so it cannot return nonsense.  
;   The line position and continuum are NOT hardwired in this program.  
;   The continuum is calculated before the fit and is only slghtly
;   modified during the fit. For the penta Gaussian fit, the continuum is
;   subtracted before the fit.
;   After that, we repeat the penta Gaussian fit with two more free parameters
;   for the strong O IV line to have an independent measure of the velocity and line width.
;   Finally, we perform a hexa Gaussian fit including a double Guassian function for the
;   Si IV 1403 line and the rest identical to the second penta Gaussian fit.
;   The last two fits are performed only if the data has a reasonable SNR.
;
;   A similar approach is used for the O I and the two C I lines at 155 nm range.  
;  
;++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
; list of output parameters:
;
; spar[0] :   baseline = continuum -------gaussian  fit
; spar[1] :   line-core position ------- gaussian fit
; spar[2] :   line-core amplitude ------- gaussian fit
; spar[3] :   line-core width ------- gaussian fit  > multiply by
;                                     2. sqrt(2. * alog(2.)) =  FWHM
; spar[4] :   chi^2 -------gaussian  fit
; spar[5:8] : sigma  
;
; dpar[0] :   baseline = continuum -------
; dpar[1] :   line-core position -------1st Gaussian component
; dpar[2] :   line-core amplitude ------- 1st Gaussian component
; dpar[3] :   line-core width ------- 1st Gaussian component
; dpar[4] :   line-core position ------- 2nd Gaussian component
; dpar[5] :   line-core amplitude ------- 2nd Gaussian component
; dpar[6] :   line-core width ------- 2nd Gaussian component
; dpar[7] :   chi^2 ------- double Gaussian  fit
; dpar[8:14] : sigma
  
; qpar[0] :   line-center position -------Si IV 1403
; qpar[1] :   line-core amplitude ------- Si IV 1403
; qpar[2] :   line-core width ------- Si IV 1403
; qpar[3] :   line-core amplitude ------- O IV 1401.16
; qpar[4] :   line-core amplitude ------- O IV 1399.97
; qpar[5] :   line-core amplitude ------- Si IV 1404.79
; qpar[6] :   line-core amplitude ------- O IV 1404.82
; qpar[7] :   chi^2 ------- double Gaussian  fit
; qpar[8:14] : sigma
  
; vpar[0] :   line-core position ------- Voigt fit
; vpar[1] :   Doppler width ------- Voigt fit
; vpar[2] :   damping width ------- Voigt fit
; vpar[3] :   amplitude ------- Voigt fit
; vpar[4] :   area = EQW ------- Voigt fit
; vpar[5] :   baseline = continuum ------- Voigt fit
; vpar[6] :   chi^2 ------- Voigt fit
; vpar[7:10] : sigma
;
; band[0:2]  : similar to H3, H2v, H2r
;
; Dec 08, 2014 : single Gaussian fit
; Dec 09, 2014 : program for double Gaussian fit using MPFIT
;                single Gaussian fit was completely re-written
; Mar 15, 2015 : penta Gaussian fit to Si IV and four other lines
; Dec 25, 2015 : improved documentation
; Apr 25, 2016 : cleaning the routine, main features are the same
; May 23, 2016 : improved hexa Gaussian fit
; Jul 22, 2016 : improving initial guess of the double Gaussian fit
; Nov 21, 2016 : random initialization of the single Gaussian fits
; Dec 19, 2016 : improved random initialization
;
; R.Rezaei @ IAC                         e-mail:  rrezaei@iac.es      
;===============================================================
;-

;on_error, 2
common share_ergm, av_fit
common share_delta, delta

if n_elements(winindex) eq 0 then windex  = 20 else windex = winindex
if n_elements(av_fit) eq 0 then av_fit = [1.1, 9., pos_line, 9.] ; continuum, amplitude, line center, line width

prof = iline
pprof = smooth(prof, 19)
x = abs(deriv(prof))           

band1 = fltarr(3)
spar1 = fltarr(9)
dpar1 = fltarr(15)
qpar1 = fltarr(15)
zpar1 = fltarr(21)
fpar1 = fltarr(27)
case disper of
   1.2980: dfac = 1.
   2.5960: dfac = 2.
   5.1920: dfac = 4.
endcase

;------------------------------------------------
;-- line core position
;------------------------------------------------
np = n_elements(prof)
corepos = pos_line
core_pos = corepos
fitsg = fltarr(np)
fitdg = fltarr(np)
fitqg = fltarr(np)
fitzg = fltarr(np)
fitfg = fltarr(np)

erg = {spar1:dblarr(9), dpar1:dblarr(15), qpar1:dblarr(15), zpar1:dblarr(21), fpar1:dblarr(27), band1:dblarr(3), sprf:dblarr(np), $
       sfit:dblarr(np), dfit:dblarr(np), qfit:dblarr(np), zfit:dblarr(np), ffit:dblarr(np)}

;----------------------------------------------------------
;-- a single Gassian fit to the line
;----------------------------------------------------------
py = reform(prof)
px = findgen(np)
bbc = where(py gt (-5.0), bppc)

if (bppc le 10) then return, erg

ee = 0.5 * sqrt(abs(py)) > 1.0 ;ir_error(py, /fuv, /dark) < py ;0.5 * sqrt(abs(py)) > 0.5
err_ave = ee < (abs(py)*0.9) ;make sure that the error value is larger than 0 and smaller than the data value
err_ave = err_ave > 1.0
posmax = pos_line
tmp = max(py, pos)
if (abs(pos - posmax) gt 5.)and(abs(pos - posmax) lt 20.) then posmax = pos


;-----------------------------
;-- define continuum windows
;-----------------------------
tmp = fltarr(np)
tmp[(posmax-560/dfac)>0:(posmax- 530/dfac)>0] = 1.0
tmp[(posmax-510/dfac)>0:(posmax- 425/dfac)>0] = 1.0
tmp[(posmax-415/dfac)>0:(posmax- 330/dfac)>0] = 1.0
tmp[(posmax-290/dfac)>0:(posmax- 270/dfac)>0] = 1.0
tmp[(posmax-200/dfac)>0:(posmax- 150/dfac)>0] = 1.0
tmp[(posmax-92/dfac)>0:(posmax- 60/dfac)>0] = 1.0
tmp[(posmax+60/dfac)<(np-1):(posmax+140/dfac)<(np-1)] = 1.0
tmp[(posmax+190/dfac)<(np-1):(posmax+210/dfac)<(np-1)] = 1.0
tmp[(posmax+275/dfac)<(np-1):(posmax+295/dfac)<(np-1)] = 1.0

qk = where(tmp eq 1)

kont = good_mean(py(qk)) < 40.0 * dfac            ; kont will be affected by CRs
fitfg = fltarr(np)
;----------------------------------------------------------
;-- remove continuum slope, if present
;----------------------------------------------------------
res = poly_fit(px(qk), py(qk), 2, /double)
base_level = res[0] + res[1]*px + res[2] * px^2 ;+ res[3] * px^3
py = py - base_level + kont

rms = stddev(py(qk))
;---------------------------------------------------------------------
;-- mark position of other lines than Si IV 1403
;-- so the single/double Gaussian fit will not care about them.
;---------------------------------------------------------------------
bad = [findgen(54/dfac)+(posmax - 146.32/dfac), findgen(40/dfac)+(posmax - 255/dfac)]
bad = [bad, findgen(40/dfac)+(posmax + 138.76/dfac), findgen(10)+(posmax + 210/dfac)]
;--------------------------------------------
;-- initial guess
;--------------------------------------------
  fit0 = [kont, max(py) * 0.8,  posmax , 15./dfac] * 1.0d
  range0 = [0.1, 0.1, 0.15, 1.] 

  tcc = where((px gt (posmax - 150/dfac)>0)and(px lt (posmax + 150/dfac)<(np-1)))
  tcc = bbc

  ppy = median(py, 3)
  ee1 = ee
  ee1[bad] = 100.0
  err_ave1 = ee1 > 0.5

  dlambda = 1.0d
  ergs = my_sgf(px[tcc], ppy[tcc], ee1[tcc], fit0, range0, dlambda[0], tcc, /double)
  fitsg = rgauss(px, [ergs.b, ergs.i1, ergs.p1,  ergs.w1])
  chisq1 = (1.0d/(n_elements(py[tcc]) - 4.0d)) * total(((ppy[tcc] - fitsg[tcc])/err_ave1[tcc])^2)
  is_bad = evaluate_sgf(fit0, ergs)

  if (ergs.w1/disper gt 25.)or(is_bad eq 1.)or(n_elements(fitsg) lt np)or(chisq1 gt 50.) then begin 
     py0 = median(py,5)
     py = gauss_smooth(py0, 1.2, /edge_truncate)
     fit0[3] = av_fit[3]+randomn(seed)*0.3
 
     range0 = [0.02,   0.3,  0.1,  0.5] 
     ergs = my_sgf(px[tcc], ppy[tcc], ee1[tcc], fit0, range0, dlambda[0], tcc, /double)
     fitsg = rgauss(px, [ergs.b, ergs.i1, ergs.p1,  ergs.w1])
     chisq1 = (1.0d/(n_elements(py[tcc]) - 4.0d)) * total(((py[tcc] - fitsg[tcc])/err_ave1[tcc])^2)
     is_bad = evaluate_sgf(fit0, ergs)
  endif  
  
  if (ergs.w1/disper gt 25.0)or(is_bad eq 1.)or(n_elements(fitsg) lt np)or(chisq1 gt 20.) then begin 
     py0 = median(py,5)
     py = gauss_smooth(py0, 1.2, /edge_truncate)
     fit0[3] = av_fit[3]
     glx = max(ppy, p_ii)
     if (p_ii gt pos_line * 0.8) then fit0[2] = p_ii
     range0 = [0.02,   0.1,            0.2, 0.5] 
     ergs = my_sgf(px[tcc], ppy[tcc], ee1[tcc], fit0, range0, dlambda[0], tcc, /double)
     fitsg = rgauss(px, [ergs.b, ergs.i1, ergs.p1,  ergs.w1])
     chisq1 = (1.0d/(n_elements(py[tcc]) - 4.0d)) * total(((py[tcc] - fitsg[tcc])/err_ave1[tcc])^2)
  endif  
  
  spar1 = [ergs.b, ergs.p1, ergs.i1,  ergs.w1, chisq1, reform(ergs.sigma)]

  ;------------------------------------------------------------------------------
  ;  if the fit is not satisfactory, then use random initialization 
  ;------------------------------------------------------------------------------
  if (chisq1 gt 1.) then begin
        random_sg_fit, 30, px, py, ppy, ee1, err_ave1, tcc, 0.6, ergs, spar_new, fitsg_new
        if (spar_new[4] lt chisq1) then begin
           chisq1 = spar_new[4]
           fitsg = fitsg_new
           spar1 = spar_new
        endif
  endif   

  line_core  = spar1[1] > (pos_line - 30.) < (pos_line + 30.)
  
;---------------------------------------------------
;--  a double Guassian fit to the line profile
;---------------------------------------------------
  fit0 = [kont, spar1[2], spar1[1], spar1[3]>.2+randomn(seed)*0.3] * 1.0d
  fit1 = [spar1[2] * (0.35 + randomu(seed)*0.1), spar1[1]-4.0/dfac, spar1[3] * 1.8+randomn(seed)*0.3] * 1.0d

  range0 = [0.1, 0.1, 0.1, 0.5] * 1.0d
  range1 = [     0.8, 0.1, 0.4] * 1.0d
  dlambda = 1.0d
  ergd = my_dgf(px[tcc],ppy[tcc], ee1[tcc], fit0, fit1, range0, range1, dlambda[0], tcc, /double)
  fitdg = ergd.b + ergd.i1*exp(-((px - ergd.p1)/ergd.w1)^2) + ergd.i2*exp(-((px - ergd.p2)/ergd.w2)^2)
  chisq2 = (1.0d/(n_elements(tcc) - 7.0d)) * total(((ppy[tcc] - fitdg[tcc])/err_ave1[tcc])^2)
  is_bad = evaluate_dgf(fit0, ergd)
  
  if (ergd.w1 gt 25./dfac)or(is_bad eq 1.)or(n_elements(fitdg) lt np)or(chisq2 gt 50.) then begin 
     py0 = median(py,5)
     py = gauss_smooth(py0, 1.2, /edge_truncate)
     fit0 = [kont,  spar1[2], spar1[1], spar1[3]+randomn(seed)*0.3]  * 1.0d
     fit1 = [spar1[2] * (0.05 + randomu(seed)*0.1),  spar1[1]-2.0/dfac,  spar1[3]* 1.5+randomn(seed)*0.3] * 1.0d  ; 1.2
     range0 = [0.05,   0.3, 0.1, 0.5]  * 1.0d
     range1 = [        0.5, 0.1, 0.5]  * 1.0d ; 0.5; 0.3 sigma
     ergd = my_dgf(px[tcc],ppy[tcc], ee1[tcc], fit0, fit1, range0, range1, dlambda[0], tcc, /double)
     fitdg = ergd.b + ergd.i1*exp(-((px - ergd.p1)/ergd.w1)^2) + ergd.i2*exp(-((px - ergd.p2)/ergd.w2)^2)
     chisq2 = (1.0d/(n_elements(tcc) - 7.0d)) * total(((ppy[tcc] - fitdg[tcc])/err_ave1[tcc])^2)
  endif
  dpar1 = [ergd.b, ergd.p1, ergd.i1, ergd.w1, ergd.p2, ergd.i2, ergd.w2, chisq2, reform(ergd.sigma)]
  
  ;------------------------------------------------------------------------------
  ;  if the fit is not satisfactory, then use random initialization 
  ;------------------------------------------------------------------------------
  if (chisq2 gt 1.) then begin
        random_dg_si_fit, 30, px, py, ppy, ee1, err_ave1, tcc, 0.6, ergd, dpar_new, fitdg_new
        if (dpar_new[7] lt chisq2) then begin
           chisq2 = dpar_new[7]
           fitdg = fitdg_new
           dpar1 = dpar_new
        endif
  endif

    
;-----------------------------------------------------------------------------
;--  a penta Guassian fit to the Si IV line profile + four other lines
;--  using mpfit package with 7 free parameters
;-----------------------------------------------------------------------------
  if (spar1[3] lt 4.0d/dfac)or(spar1[3] gt 52.0d/dfac) then spar1[3] = 10.0d/dfac
  if (abs(spar1[1] - pos_line) gt 20.) then spar1[1] = pos_line

  fit0 = [spar1[2] > 5., spar1[1], spar1[3] > .1]
  fit1 = [max(fitdg)*0.244d, max(fitdg)*0.065d, max(fitdg)*0.065d, max(fitdg)*0.065d] 

  range0=[0.1, 0.1, 0.8] 
  range1=[0.5, 0.5, 0.5, 0.5]

  dlambda = 1.0d
  pyn = ppy - dpar1[0]
  ergq = my_tgf(px[bbc], pyn[bbc], ee[bbc], fit0, fit1, range0, range1, dlambda[0], bbc, /double)
  fitqg = pentagauss(px, [ergq.i1, ergq.p1, ergq.w1, ergq.i2, ergq.i3, ergq.i4, ergq.i5])
  chisq4 = (1.0d/(n_elements(bbc) - 7.0d)) * total(((pyn[bbc] - fitqg[bbc])/err_ave[bbc])^2)
  fitqg += dpar1[0]

  if (ergq.w1 gt 25./dfac)or(ergq.w1/disper lt 0.5)or(n_elements(fitqg) lt np)or(chisq4 gt 50.) then begin 
     py0 = median(py,5)
     py = gauss_smooth(py0, 1.2, /edge_truncate)
     fit0[2] = av_fit[3]+randomn(seed)*0.3
     range0 = [0.05,  0.1, 0.5] * 1.0d
     range1 = [0.3,  0.3, 0.3, 0.3] * 1.0d
     ergq = my_tgf(px[bbc], pyn[bbc], ee[bbc], fit0, fit1, range0, range1, dlambda[0], bbc, /double)
     fitqg = pentagauss(px, [ergq.i1, ergq.p1, ergq.w1, ergq.i2, ergq.i3, ergq.i4, ergq.i5])
     chisq4 = (1.0d/(n_elements(bbc) - 7.0d)) * total(((pyn[bbc] - fitqg[bbc])/err_ave[bbc])^2)
     fitqg += dpar1[0]
  endif

  qpar1=[ergq.p1, ergq.i1, ergq.w1, ergq.i2, ergq.i3, ergq.i4, ergq.i5, chisq4, reform(ergq.sigma)]

;---------------------------------------------------------------------------------
;-- a  checkpoit for LVG (>>>>>>  experimental  <<<<<<)
;---------------------------------------------------------------------------------
  delta = 0.
if 0 then begin
  fit0 = [ergq.i1, ergq.p1, ergq.w1]
  fit1 = [ergq.i2, ergq.i3, ergq.i4, ergq.i5, 12.0 / dfac] ; initial guess

  range = [1.0] 

  dlambda = 0.8d
  pyn = ppy - dpar1[0]
  ergt = lvg(px[bbc], pyn[bbc], ee[bbc], fit0, fit1, range, dlambda[0], bbc, /double)
  fitqg = pentagauss_lvg(px, [ergt.p2])
  chisq4_1 = (1.0d/(n_elements(bbc) - 1.0d)) * total(((pyn[bbc] - fitqg[bbc])/err_ave[bbc])^2)
  fitqg += dpar1[0]
  delta = ergt.p2
endif
  
if (max(fitdg)/rms gt 5.) then begin ; if there is a little bit of signal
;---------------------------------------------------------------------------------
;--  a penta Guassian fit to the Si IV line profile + four other lines
;--  using mpfit package with 10 free parameters (two more for O IV 1401)
;---------------------------------------------------------------------------------
;---------------------------------------------------------------------------------
;-- initial guess for each component
;-- using results of the first penta Gaussian fit
;---------------------------------------------------------------------------------
  o_iv_pos = ergq.p1 - 126.32 / dfac

  o_iv_w = ergq.w1 *1.22 +randomn(seed)*0.3
  fit0 = [ergq.i1, ergq.p1, ergq.w1+randomn(seed)*0.3] * 1.0d
  fit1 = [ergq.i2>2., o_iv_pos, o_iv_w, ergq.i3>1., ergq.i4>0.2, ergq.i5 > 0.2, ergq.i2 * 0.3] * 1.0d ; initial guess

  range0=[0.1, 0.1, 0.1] * 1.0d
  range1=[0.8, 0.1, 0.3, 0.5, 0.5, 0.5, 0.9] * 1.0d
  if (max(pyn) le 10.) then range1 = [0.8, 0.1, 0.1, 0.5, 0.5, 0.5, 0.5] * 1.0d

  dlambda = 1.0d
  pyn = ppy - dpar1[0]
  ergz = my_pgf(px[bbc], pyn[bbc], ee[bbc], fit0, fit1, range0, range1, dlambda[0], bbc, /double)
  fitzg = pentagauss2(px, [ergz.i1, ergz.p1, ergz.w1, ergz.i2, ergz.p2, ergz.w2, ergz.i3, ergz.i4, ergz.i5, ergz.i6])
  chisz5 = (1.0d/(n_elements(bbc) - 10.0d)) * total(((pyn[bbc] - fitzg[bbc])/err_ave[bbc])^2)
  fitzg += dpar1[0]

  if (ergz.w1 gt 25./dfac)or(ergz.w1 lt 0.5/dfac)or(n_elements(fitzg) lt np)or(chisz5 gt 10.) then begin 
     fit0 = [spar1[2], spar1[1],  spar1[3]+randomn(seed)*0.3]
     range1 = [0.5, 0.1, 0.2, 0.3, 0.3, 0.3, 0.6] * 1.0d
     ergz = my_pgf(px[bbc], pyn[bbc], ee[bbc], fit0, fit1, range0, range1, dlambda[0], bbc, /double)
     fitzg = pentagauss2(px, [ergz.i1, ergz.p1, ergz.w1, ergz.i2, ergz.p2, ergz.w2, ergz.i3, ergz.i4, ergz.i5, ergz.i6])
     chisz5 = (1.0d/(n_elements(bbc) - 10.0d)) * total(((pyn[bbc] - fitzg[bbc])/err_ave[bbc])^2)
     fitzg += dpar1[0]
  endif
  
  zpar1=[ergz.p1, ergz.i1, ergz.w1, ergz.p2, ergz.i2, ergz.w2, ergz.i3, ergz.i4, ergz.i5, ergz.i6, chisz5, reform(ergz.sigma)]

;---------------------------------------------------------------------------------
;--  a hexa-Guassian fit to the Si IV line profile + four other lines
;--  using mpfit package with 13 free parameters: 
;--  double Gausian fit for Si IV (6) and single Gaussian fit for
;--  O IV 1401 (3), and three peak intensities for other lines.
;---------------------------------------------------------------------------------
;---------------------------------------------------------------------------------
;-- initial guess for each component
;-- using results of the second penta Gaussian fit
;---------------------------------------------------------------------------------
  fit0 = [dpar1[2], dpar1[1], dpar1[3]>.2+randomn(seed)*0.3, dpar1[5]>2., dpar1[4], dpar1[6]>.2+randomn(seed)*0.3] * 1.0d ; Si IV double Gaussian
  fit1 = [ergz.i2, ergz.p2, ergz.w2>.2+randomn(seed)*0.3, ergz.i3>1., ergz.i4>0.2, ergz.i5>0.2, ergz.i6>0.1]* 1.0d
  range0=[0.1, 0.1, 0.1, 0.1, 0.1, 0.1] * 1.0d
  range1=[0.1, 0.1, 0.5, 0.5, 0.5, 0.5, 0.5] * 1.0d

  dlambda = 1.0d
  pyn = ppy - dpar1[0]
  ergf = my_hgf(px[bbc], pyn[bbc], ee[bbc], fit0, fit1, range0, range1, dlambda[0], bbc, /double)
  fitfg = hexagauss(px, [ergf.i1, ergf.p1, ergf.w1, ergf.i2, ergf.p2, ergf.w2,ergf.i3,ergf.p3,ergf.w3,ergf.i4,ergf.i5,ergf.i6,ergf.i7])
  chisf4 = (1.0d/(n_elements(bbc) - 13.0d)) * total(((pyn[bbc] - fitfg[bbc])/err_ave[bbc])^2)
  fitfg += dpar1[0]

  if (ergf.w1/disper gt 25.)or(ergf.w1/disper lt 0.5)or(n_elements(fitfg) lt np)or(chisf4 gt 10.) then begin ; something is went wrong
     range1=[0.1, 0.1, 0.3, 0.2, 0.2, 0.2, 0.2] * 1.0d
     ergf = my_hgf(px[bbc], pyn[bbc], ee[bbc], fit0, fit1, range0, range1, dlambda[0], bbc, /double)
     fitfg = hexagauss(px, [ergf.i1, ergf.p1, ergf.w1, ergf.i2, ergf.p2, ergf.w2, ergf.i3, ergf.p3, ergf.w3, ergf.i4, ergf.i5, ergf.i6,ergf.i7])
     chisf4 = (1.0d/(n_elements(bbc) - 13.0d)) * total(((pyn[bbc] - fitfg[bbc])/err_ave[bbc])^2)
     fitfg += dpar1[0]
  endif
  
  fpar1=[ergf.p1, ergf.i1, ergf.w1, ergf.p2, ergf.i2, ergf.w2, ergf.p3, ergf.i3, ergf.w3, ergf.i4, ergf.i5, ergf.i6, ergf.i7, chisf4, reform(ergf.sigma)]
endif  

;----------------------------
;--  band parameters
;----------------------------
band1[0] = delta ; total(prof[line_core-4:line_core+4])* disper   ; H3 similar to band intensity in Ca II and Mg II
band1[1] = total(prof[line_core-12:line_core-3])* disper  ; H2v similar to band intensity in Ca II and Mg II
band1[2] = total(prof[line_core+5:((line_core+14)< (np-1))])* disper  ; H2r  similar to band intensity in Ca II and Mg II

;-------------------------------------
;-- plot observed and fit profiles
;-------------------------------------
if (plt eq 1) then begin
  print, '------------------------------------------------------------------------'
  print, 'single Gaussian:', spar1[0:4]
  print, 'double Gaussain: ', dpar1[0:7]
  print, 'penta Gaussain (I): ',qpar1[0:7]
  print, 'penta Gaussain (II):', zpar1[0:9]
  print, 'hexa Gaussain: ',fpar1[0:12]
  print, '------------------------------------------------------------------------'
  plot, px, py, /xsty;, /ylog, yrange = [0.1, max(py)], /ystyle
  loadct, 40, /silent
  oplot, px, fitsg, thick=2, color=70
  oplot, px, fitdg, linestyle=2, thick=2, color=175
  oplot, px, fitqg, color=215
  oplot, px, fitfg, linestyle=2, thick=2, color=254
  loadct, 0, /silent
  xyouts, 0.15, .7, num2string(round(max(fitdg)/rms)), chars=2,/normal
  ans =''  &  read, ans, prompt='Press enter to continue ...'
endif

;----------------------------
;--  Set output parameters
;----------------------------
erg.spar1 = spar1         ; single Gaussian fit with 4 free parameters
erg.dpar1 = dpar1         ; double Gaussian fit with 7 free parameters
erg.qpar1 = qpar1         ; penta Gaussian fit with 7 free parameters
erg.zpar1 = zpar1         ; penta Gaussian fit with 9 free parameters
erg.fpar1 = fpar1         ; penta Gaussian fit with 11 free parameters
erg.band1 = band1         ; 3
erg.sprf  = py            ; observed line profile
erg.sfit  = fitsg         ; single Gaussian fit (4) 
erg.dfit  = fitdg         ; double Gaussian fit (7)
erg.qfit  = fitqg         ; penta Gaussian fit (7)
erg.zfit  = fitzg         ; penta Gaussian fit (10)
erg.ffit  = fitfg         ; hexa Gaussian fit (13)


return,erg
end
