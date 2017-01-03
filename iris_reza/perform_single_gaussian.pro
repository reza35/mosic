function perform_single_gaussian, iline,ca_ilo,ca_ihi, pos_line, disper, plt

;+
;===============================================================
; function : perfrm_single_gaussian.pro
;
; purpose : calculates parameters of a single "narrow" emission line
; observed by IRIS, e.g., O I and Cl I line 
;
; iline: a line profile
;
;
; - We use Single Gaussian fit. In complex situations, emission lines
;   show self absorbtion, so will be like Ca II lines. That deserves a
;   double or triple Gaussian fit which exists  in gen/tian folder.
;
; list of output parameters:
;
; hpar[0] :   line-core position ------- gaussian fit
; hpar[1] :   line-core width ------- gaussian fit  > multiply by
;                                     2. sqrt(2. * alog(2.)) =  FWHM
; hpar[2] :   line-core amplitude ------- gaussian fit
; hpar[3] :   line-core area ------- gaussian fit
; hpar[4] :   line-core area (observed) = EWQ
; hpar[5] :   baseline = continuum -------gaussian  fit
; hpar[6] :   1/e width of the Fit profile ------- gaussian fit
; hpar[7] :   chi^2 -------gaussian  fit
; hpar[8] :   error in line position -------gaussian  fit
; hpar[9] :   error in line width ------- gaussian fit
; hpar[10] :  error in amplitude ------- gaussian fit
; hpar[11] :  error in baseline ------- gaussian fit
;
; band[0:2] : H3, H2v, H2r, each 9 pixel wide
;
; Aug 12, 2014 : output for single Gaussian fit
; Jun 18, 2015 : updated and cleaned extra parameters
; May 01, 2016 : unified scheme, continuum slope  
;
; R.Rezaei @ IAC                         e-mail:  rrezaei@iac.es      
;===============================================================
;-


on_error, 2
common share_ergm, av_fit

if (n_elements(plt) eq 0) then plt = 0
if n_elements(winindex) eq 0 then windex  = 20 else windex = winindex

prof = iline
pprof = smooth(prof, 19)
x = abs(deriv(prof))           

band1 = fltarr(3)
spar1 = fltarr(9)
dpar1 = fltarr(15)

;------------------------------------------------
;-- line core position
;------------------------------------------------
np = n_elements(prof)
corepos = pos_line
posmax = pos_line
;----------------------------------------------------------
;-- a single Gassian fit to the line
;----------------------------------------------------------
py = reform(prof)
np = n_elements(py)
px = findgen(np)

;----------------------------------------------------------
;-- remove continuum slope, if present
;----------------------------------------------------------
tmp = fltarr(np)
if (disper lt 2.) then begin
  tmp[(posmax-20):(posmax+20)]=1.
endif else begin
   tmp[(posmax-10):(posmax+10)]=1.
endelse

q = where(tmp eq 0) 
kont = good_mean(py(q))< 40.  ; kont will be affected by CRs
res = poly_fit(px(q), py(q), 2, /double)
base_level = res[0] + res[1]*px + res[2] * px^2
py = py - base_level + kont
  
bpp = where(py gt (-10.0), bppc)

ee =  0.5 * sqrt(abs(py))  > 0.5    ; see ir_error.pro for a description
err_ave = ee < (abs(py)*0.9) ;make sure that the error value is larger than 0 and smaller than the data value
err_ave = err_ave > 0.5


;--------------------------------------------
;-- initial guess
;--------------------------------------------
if (disper gt 2.) then fit0 = [kont, max(py) * 0.9,  posmax , 2. + randomn(seed)*.2] * 1.0d
if (disper lt 1.4) then fit0 = [kont, max(py) * 0.9,  posmax , 3.1 + randomn(seed)*.2] * 1.0d
  
range0 = [0.05, 0.2, 0.2, 0.8] 
dlambda = 1.0d
  
ergs = my_sgf(px[bpp], py[bpp], ee[bpp], fit0, range0, dlambda[0], bpp, /double)
fitsg = ergs.b + ergs.i1*exp(-((px - ergs.p1)/ergs.w1)^2)
chisq1 = (1.0d/(n_elements(py[bpp]) - 4.0d)) * total(((py[bpp] - fitsg[bpp])/err_ave[bpp])^2)
is_bad = evaluate_sgf(fit0, ergs)
if (abs(ergs.p1 - posmax) gt 10.) then is_bad = 1.  ; we know that Cl I cannot have large Doppler shift

if (ergs.w1 gt 7.)or(is_bad eq 1.)or(n_elements(fitsg) lt np)or(chisq1 gt 10.) then begin ; something went wrong
     py0 = median(py, 5)
     py = gauss_smooth(py0, 1.2, /edge_truncate)
     fit0 = [av_fit[0], av_fit[1]>2., av_fit[2],  av_fit[3]+randomn(seed)*0.1]
     range0=[0.02,   0.3,            0.1, 0.5] 
     ergs = my_sgf(px[bpp], py[bpp], ee[bpp], fit0, range0, dlambda[0], bpp, /double)
endif   

fitsg = ergs.b + ergs.i1*exp(-((px - ergs.p1)/ergs.w1)^2)
chisq1 = (1.0d/(n_elements(py[bpp]) - 4.0d)) * total(((py[bpp] - fitsg[bpp])/err_ave[bpp])^2)

spar1=[ergs.b, ergs.p1, ergs.i1,  ergs.w1, chisq1, reform(ergs.sigma)]
;print, spar1[0:5]
;-------------------------------------
;-- plot lines, positions of peaks
;-------------------------------------
if (plt eq 1) then begin
  print, '------------------------------------------------------'
  print, 'single Gaussain: ', spar1[0:4]
  print, '------------------------------------------------------'
  plot, px, py, /xsty
  loadct, 40, /silent
  oplot, px, fitsg, thick=2, color=250
  loadct,0,/silent
endif

;----------------------------
;--  Set output parameters
;----------------------------
erg = {spar:fltarr(9), sprf:fltarr(np), sfit:fltarr(np)}

erg.spar = spar1        ; single Gaussian fit parameters
erg.sprf = py           ; line profile
erg.sfit = fitsg        ; single Gaussian fit to the line profile

return,erg
end
