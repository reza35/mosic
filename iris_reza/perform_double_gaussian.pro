function perform_double_gaussian, tag, iline, ca_ilo, ca_ihi, pos_line, disper, plt

;+
;===============================================================
; function : perform_double_gaussian.pro
;
; purpose : performs a single and double Gaussian fit to a single 
;           emisison line 
;
; tag : name of the spectral line, should be either 'C' or 'Si'
; iline: a line profile
; ca_ili/ca_ihi : the lower/upper range
; pos_line : line center position
; disper : spectral dispersion in mA/px
; plt=1 plots the line and the fit
;
;
; - We use Single Gaussian fit. In complex situations, emission lines
;   show self absorbtion, so will be like Ca II lines. That deserves a
;   double or triple Gaussian fit which will be performed right after
;   the singel Gaussian fit.
;   The continuum is calculated befor the fit and is only slghtly
;   modified during teh fit.  
;
; list of output parameters:
;
; spar[0] :   baseline = continuum -------gaussian  fit
; spar[1] :   line-core position ------- gaussian fit
; spar[2] :   line-core amplitude ------- gaussian fit
; spar[3] :   line-core width ------- gaussian fit  > multiply by
;                                     2. sqrt(2. * alog(2.)) =  FWHM
  
; spar[4] :   chi^2 ------- single gaussian  fit
; spar[5] :   error in baseline ------- gaussian fit
; spar[6] :   error in line position -------gaussian  fit
; spar[7] :   error in amplitude ------- gaussian fit
; spar[8] :   error in line width ------- gaussian fit
;
; dpar[0] :   baseline = continuum -------is not a fit parameter
; dpar[1] :   line-core position -------1st Gaussian component
; dpar[2] :   line-core amplitude ------- 1st Gaussian component
; dpar[3] :   line-core width ------- 1st Gaussian component
; dpar[4] :   line-core position ------- 2nd Gaussian component
; dpar[5] :   line-core amplitude ------- 2nd Gaussian component
; dpar[6] :   line-core width ------- 2nd Gaussian component
; dpar[7] :   chi^2 ------- double Gaussian  fit
; dpar[8:14] : sigma of the first six pars.
  
; vpar[0] :   line-core position ------- Voigt fit
; vpar[1] :   Doppler width ------- Voigt fit
; vpar[2] :   damping width ------- Voigt fit
; vpar[3] :   amplitude ------- Voigt fit
; vpar[4] :   area = EQW ------- Voigt fit
; vpar[5] :   baseline = continuum ------- Voigt fit
; vpar[6] :   chi^2 ------- Voigt fit
; vpar[7] :   error in line position ------- Voigt fit
; vpar[8] :   error in Doppler width ------- Voigt fit
; vpar[9] :   error in damping width ------- Voigt fit
; vpar[10] :  error in baseline ------- Voigt fit
;
; band[0:2] : similar to H3, H2v, H2r
;
; Dec 08, 2014 : output for single Gaussian fit
; Dec 09, 2014 : program for double Gaussian fit using MPFIT
;                single Gaussian fit was completely re-written
;
; Apr 22, 2016 : improved continuum estimate
; Apr 25, 2016 : better fits for optically thick lines
; May 01, 2016 : unified scheme, continuum slope  
; Nov 22, 2017 : modified smoothing scheme  
;  
; R.Rezaei @ IAC                         e-mail:  rrezaei@iac.es      
;===============================================================
;-
;on_error, 2

if (n_elements(plt) eq 0) then plt = 0
if (n_elements(tag) eq 0) then begin
   print, 'Line Tag is missing, e.g. Si, or C'
   stop
endif   
if n_elements(winindex) eq 0 then windex  = 20 else windex = winindex

prof = iline
pprof = gauss_smooth(prof, 1.5, /edge_truncate)
x = abs(deriv(prof))           

band1 = fltarr(3)
spar1 = fltarr(9)
dpar1 = fltarr(15)

;------------------------------------------------
;-- line core position
;------------------------------------------------
np = n_elements(prof)
corepos = pos_line
core_pos = corepos

;----------------------------------------------------------
;-- remove continuum slope, if present
;----------------------------------------------------------
py = reform(prof)
px = findgen(np)

q = where(py lt abs(median(py)) * 1.5)
kont = good_mean(py(q)) < 40.0 ; kont will be affected by CRs
res = poly_fit(px(q), py(q), 2, /double)
base_level = res[0] + res[1]*px + res[2] * px^2
py = py - base_level + kont

;----------------------------------------------------------
;-- a single Gassian fit to the line
;----------------------------------------------------------
bpp = where(py gt (-10.), bppc)

ee =  0.5 * sqrt(abs(py)) > 0.5      ; see ir_error.pro for a description
err_ave = ee < (abs(py)*0.9) ;make sure that the error value is larger than 0 and smaller than the data value
err_ave = err_ave > 0.5
posmax = pos_line
 
;--------------------------------------------
;-- initial guess
;--------------------------------------------
if (disper gt 2.)and(tag eq 'C') then fit0 = [kont, max(py) * 0.9,  posmax , 5.] * 1.0d
if (disper gt 2.)and(tag eq 'Si') then fit0 = [kont, max(py) * 0.9,  posmax , 7.] * 1.0d
if (disper lt 1.4)and(tag eq 'C') then fit0 = [kont, max(py) * 0.9,  posmax , 10.] * 1.0d
if (disper lt 1.4)and(tag eq 'Si') then fit0 = [kont, max(py) * 0.9,  posmax , 15.] * 1.0d

fit0[3] += randomn(seed)*0.1
range0 = [0.05, 0.75, 0.2, .99] 
dlambda = 1.0d

ergs = my_sgf(px[bpp], py[bpp], ee[bpp], fit0, range0, dlambda[0], bpp, /double)
fitsg = ergs.b + ergs.i1*exp(-((px - ergs.p1)/ergs.w1)^2)
chisq1 = (1./(n_elements(py[bpp]) - 4.)) * total(((py[bpp] - fitsg[bpp])/err_ave[bpp])^2)
is_bad = evaluate_sgf(fit0, ergs)

if (ergs.w1 gt 25.)or(is_bad eq 1.)or(n_elements(fitsg) lt np)or(chisq1 gt 10.) then begin ; something went wrong
     py0 = median(py, 5)
     py = gauss_smooth(py0, 1.2, /edge_truncate)
     fit0[1] = max(py0) & fit0[3] += randomn(seed)*0.5
     range0=[0.02,   0.3,            0.1, 0.9] 
     ergs = my_sgf(px[bpp], py[bpp], ee[bpp], fit0, range0, dlambda[0], bpp, /double)
endif

fitsg = ergs.b + ergs.i1*exp(-((px - ergs.p1)/ergs.w1)^2)
chisq1 = (1./(n_elements(py[bpp]) - 4.)) * total(((py[bpp] - fitsg[bpp])/err_ave[bpp])^2)

spar1=[ergs.b, ergs.p1, ergs.i1,  ergs.w1, chisq1, reform(ergs.sigma)]
line_core  = ergs.p1 > 12 < (pos_line + 30.)

;---------------------------------------------------
;--  a double Guassian fit to the line profile
;--  using mpfit package
;---------------------------------------------------
;-- initial guess for each component
;--  using results of single Gaussian fit
;--------------------------------------------
test = ergs.p1   &  z = max(smooth(py,7), posmax)  & amp = ergs.i1  & width = ergs.w1
if (abs(test - pos_line) gt 5.) and (abs(posmax - test) gt 5.) and(z gt 15.) then begin
     test = posmax              ; if the single Gaussian fit failed
     amp = max(smooth(py,7))
     if (width gt 12.) then width = 12.
endif

  ; the inital guess for double Guassian fit
if (tag eq 'C') then begin
     fit0 = [1.1, amp * 0.8,  test * 0.98, width * 0.9] * 1.0d
     fit1 = [1.1, amp * 0.1,  test * 1.02, width * 1.3] * 1.0d

     range0=[0.05, 0.5,  0.1, .5]     ; the first component should be close to single Gaussian fit
     range1=[2.0,  0.2, .9]     ; 2.0 for amplitude means it can also be negative (for thick lines)
endif
if (tag eq 'Si') then begin
    fit0 = [1.1, amp * 0.8,  test * 1.01, width * 1.01] * 1.0d
    fit1 = [     amp * 0.2,  test+2,      width * 3.50] * 1.0d ; the second component is broad

    range0=[0.05, 0.75, 0.3, 0.5] 
    range1=[0.8,  0.8, 1.]
endif
  
dlambda = 1.0d
ergd = my_dgf(px[bpp],py[bpp], ee[bpp], fit0, fit1, range0, range1, dlambda[0], bpp, /double)
fitdg = ergd.b + ergd.i1*exp(-((px - ergd.p1)/ergd.w1)^2) + ergd.i2*exp(-((px - ergd.p2)/ergd.w2)^2)
chisq2 = (1./(n_elements(bpp) - 7.)) * total(((py[bpp] - fitdg[bpp])/err_ave[bpp])^2)
dpar1=[ergd.b, ergd.p1, ergd.i1, ergd.w1, ergd.p2, ergd.i2, ergd.w2, chisq2, reform(ergd.sigma)]

;------------------------------------------------------
;--  band parameters
;-- similar to band intensity in Ca II and Mg II
;------------------------------------------------------
if (abs(line_core - pos_line) gt 10./disper) then line_core = pos_line 
band1[0] = total(prof[(line_core-4)>0:line_core+4])* disper   ; H3 
band1[1] = total(prof[(line_core-12)>0:(line_core-3)>0])* disper  ; H2v 
band1[2] = total(prof[((line_core+5)< (np-1)):((line_core+14)< (np-1))])* disper  ; H2r  

;-------------------------------------
;-- plot lines, positions of peaks
;-------------------------------------
if (plt eq 1) then begin
  plot, px, py, /xsty
  loadct, 40, /silent
  oplot, px, fitsg, thick=2, color=70
  oplot, px, fitdg,color=250
  loadct, 0, /silent
  ;read, ans, prompt='Press enter to continue ...'
endif

;----------------------------
;--  Set output parameters
;----------------------------
erg = {spar1:fltarr(9), dpar1:fltarr(15), band1:fltarr(3), sprf:fltarr(np), $
       sfit:fltarr(np), dfit:fltarr(np)}

erg.spar1 = spar1         ; single Gaussian fit parameters
erg.dpar1 = dpar1         ; double Gaussian fit parameters
erg.band1 = band1         ; 3
erg.sprf  = py            ; line profile
erg.sfit  = fitsg         ; single Gaussian fit to the line profile
erg.dfit  = fitdg         ; double Gaussian fit to the line profile

return,erg
end
