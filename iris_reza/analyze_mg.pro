pro cmin_cmax, zmins, zmaxs, cmin, cmax, xmin, xmax

  tmp_xmin = reform(zmins[*,0])
  tmp_ymin = reform(zmins[*,1])
  tmp_ymin = tmp_ymin(sort(tmp_xmin))
  tmp_xmin = tmp_xmin(sort(tmp_xmin))
  i=0
  while (tmp_xmin[i] eq 0)and(i lt (n_elements(tmp_xmin)-1)) do i+=1
  tmp_xmin = tmp_xmin[i:*]
  u = fltarr(n_elements(tmp_xmin), 2)
  u[*,0] = tmp_xmin
  u[*,1] = tmp_ymin[i:*]
  zmins = u

  tmp_xmax = reform(zmaxs[*,0])
  tmp_ymax = reform(zmaxs[*,1])
  tmp_ymax = tmp_ymax(sort(tmp_xmax))
  tmp_xmax = tmp_xmax(sort(tmp_xmax))
  i=0
  while (tmp_xmax[i] eq 0)and(i lt (n_elements(tmp_xmax)-1)) do i+=1
  tmp_xmax = tmp_xmax[i:*]
  u = fltarr(n_elements(tmp_xmax), 2)
  u[*,0] = tmp_xmax
  u[*,1] = tmp_ymax[i:*]
  zmaxs = u

  tmp_xmin = reform(zmins[*,0])
  aux = where(tmp_xmin gt 0., count)
  if (count ge 1) then begin
     u = tmp_xmin(aux)
     xmin = u(sort(u))
     cmin = count
  endif else begin
     cmin = 0
     xmin = 0
  endelse
  
  tmp_xmax = reform(zmaxs[*,0])
  aux = where(tmp_xmax gt 0., count)
  if (count ge 1) then begin
     u = tmp_xmax(aux)
     xmax = u(sort(u))
     cmax = count
  endif else begin
     cmax = 0
     xmax = 0
  endelse
end
;===============================================================




function analyze_mg, line_profile, konti, ca_ilo, ca_ihi, master, $
                     disper, zeta, plt=plt, do_H_line=do_H_line, do_gauss=do_gauss
;+
;===============================================================
; function : analyze_mg.pro
;  
; purpose : calculates  Mg II h/k line parameters observed by IRIS
;
; history : It was inspired by my Ca II H analysis tools written for
;           POLIS data. For a definition of terms like H2v, etc, see e.g.
;           Cram & Dame 1983, ApJ 272, 355
;           Lites, Rutten, & Kalkofen 1993, ApJ 414, 345
;           Rezaei, Schlichenmaier, Beck, et al. 2007, A&A 466, 1131
;
;
;
; line_profile: observed line profile
; konti: a wing intensity
;
; ca_ilo / ca_ihi : two parameters array of line positions
;                   it contains locations of both h and k lines
; The program by default analyzes the k line. to analyze the h line,
; use the /do_H_line keyword.
;  
; master : array with ranges of k1v and k1r or the same for the h line
; plt = 1, plot profile and fits  
; plt = 2, also print out processing steps for debugging
; disper : spectral dispersion in mA
; gauss = 1 performs a Gaussian fit to the line core region
;  
; list of output parameters:
;
; hpar[0] :   H1 emission width (pm)
; hpar[1] :   H2 emission width (pm)
; hpar[2] :   W-B emission width (pm)
; hpar[3] :   W-B mean intensity in unit of wing intensity
; hpar[4] :   H2v/H2r    Violet/Red ratio
; hpar[5] :   H2v/H3     H2 emission strength
;----------------------------------------
;- Gaussian parameters
; spar : single Gaussian
; dpar : double Gaussian
; tpar : triple Guassian 
;----------------------------------------
; spar[0] :   baseline = continuum  ------- gaussian fit
; spar[1] :   line-core position ------- gaussian fit  
; spar[2] :   line-core amplitude ------- gaussian fit
; spar[3] :   line-core width ------- gaussian fit > multiply by 2. sqrt(2. * alog(2.)) =  FWHM
; spar[4] :   chi^2 ------- single Gaussian  fit
; spar[5:8] : sigma  
;
; dpar[0] :   line-core position -------1st Gaussian component
; dpar[1] :   line-core amplitude ------- 1st Gaussian component
; dpar[2] :   line-core width ------- 1st Gaussian component
; dpar[3] :   line-core position ------- 2nd Gaussian component
; dpar[4] :   line-core amplitude ------- 2nd Gaussian component
; dpar[5] :   line-core width ------- 2nd Gaussian component
; dpar[6] :   chi^2 ------- double Gaussian  fit
; dpar[8:14] : sigma
  
; tpar[0] :   line-center position ------- H3
; tpar[1] :   line-core amplitude ------- H3
; tpar[2] :   line-core width ------- H3
; tpar[3] :   emisison peak amplitude ------- H2v
; tpar[4] :   emission peak position ------- H2v
; tpar[5] :   emission peak width ------- H2v
; tpar[6] :   emisison peak amplitude ------- H2r
; tpar[7] :   emission peak position ------- H2r
; tpar[8] :   chi^2 ------- double Gaussian  fit
; tpar[9:18] : sigma
;
; velpos :   position of two photospheric lines
; fe_int :   core intensity of the photospheric lines
;
; band[0:5]   : different wing intensities
; band[6:8]   : H3, H2v, H2r
; band[9:14]  : H-index(1.0 & 0.5 A), other intensity parameters
; band[15]    : center-of-gravity in 1.0 A
; band[16]    : intensity of center-of-gravity in 1.0 A
; band[17]    : center-of-gravity in 0.5 A
; band[18]    : intensity of center-of-gravity in 0.5 A
; and the same order for h line from 19 to 29, if it was set.
;  
; emb        : blends around Mg II k line which 'can' go into emission   
; emb[0]     : continuum for Mg II 2798
; emb[1]     : core intensity for Mg II 2798, either emission or absorption
; emb[2]     : EWQ in units of I_c * pm
; emb[3]     : continuum for  2793.8
; emb[4]     : core  intensity
; emb[5]     : EWQ
; emb[6]     : continuum for  2797.0
; emb[7]     : core  intensity
; emb[8]     : EWQ
; emb[9]     : continuum for  2799.2
; emb[10]    : core  intensity
; emb[11]    : EWQ
; embv       : line-core position of the four lines  
;--------------------------------------------------------------
; Nov 03, 2014 : bug fix for very dark umbra
; Dec 09, 2014 : single Gaussian fit for all profiles
;                fixed wavelength range for the h line in band parameters
; Dec 14, 2014 : improved analysis of Mg lines for different spectral sampling
; Feb 09, 2015 : bug fix for redundant maxima
; Dec 25, 2015 : the dispersion also affects selection of the spectral range.
; Dec 27, 2015 : analysis of weak lines which go into emission like Mg II 2798.0
; Jan 14, 2016 : improved for case of dispersion = 5 pm
; May 19, 2016 : bug fix in Gaussian fit
; Jun 17, 2016 : remove the wing slope before the Gaussian fitting
; Jul 01, 2016 : the Gaussian fit was improved. The distribution of
;                the reduced chi_square forms the expected curve around one.
; Nov 12, 2016 : random initialization of Gaussian fits was improved.
; Sep 21, 2017 : improved k3 position in noisy data.
; Nov 20, 2017 : zeta controls the random processes.
;  
; R.Rezaei @ IAC                         e-mail:  rrezaei@iac.es      
;===============================================================
;-
;on_error, 2

if n_elements(plt) eq 0 then plt=0

iline = reform(line_profile)
iline2 = reform(line_profile)

q=where(~finite(iline), count) &  if (count ge 1) then stop

if n_elements(do_H_line) eq 0 then do_H_line = 0 else do_H_line = 1
if n_elements(zeta) eq 0 then zeta = 1
if (do_gauss ne 0) then do_gauss = 1 
do_improve = 1

iline = iline / konti
prof = iline 
band = fltarr(30)
hpar = fltarr(6)

minh1 = fltarr(2,2)
emb = fltarr(12)
embv = fltarr(4)
spar1 = fltarr(9)
dpar1 = fltarr(15)
tpar1 = fltarr(21)

case disper of
   2.544: dfac = 1.
   5.088: dfac = 2.
   10.176: dfac = 4.
endcase

threshold = 0.7
;---------------------------------------------------------
;--  Find reference Fe I lines in the observed spectrum 
;---------------------------------------------------------
velpos = fltarr(n_elements(ca_ilo))
fe_int = fltarr(n_elements(ca_ilo))
for k = 0, n_elements(ca_ilo)-1 do begin
  tmp1 = prof[ca_ilo[k]:ca_ihi[k]]
  lpff, tmp1, pos
  if (~finite(pos)) then pos = (ca_ihi[k] - ca_ilo[k]) * 0.5  ; in case of dark umbra with no continuum
  pos = pos + ca_ilo[k]
  tmp = round(pos) + findgen(5) - 2
  linecut = prof(tmp)
  res = poly_fit(tmp, linecut, 2, /double, yfit=yfit)
  fe_int[k] = res[0] + res[1]*pos +  res[2]*(pos^2)
  ;---------------------------------
  velpos[k] = pos
endfor
mh1v1 = master[0]    &    mh1v2 = master[1]
mh1r1 = master[2]    &    mh1r2 = master[3]
nk1v = mh1v2 - mh1v1 + 1.
nk1r = mh1r2 - mh1r1 + 1.

x = smooth(prof, 3)
lpff, x[mh1r1:mh1r2],pos
if (~finite(pos)) then pos = (mh1r2 - mh1r1) * 0.5
minh1[1,0] = pos + float(mh1r1)
vk = interpol(reform(x[mh1r1:mh1r2]), findgen(nk1r), findgen(nk1r*20.)/20.)
minh1[1,1] =  vk[round(pos * 20.)] * konti
;-------------------------------------------------
p=mh1v2                    
while (x[p-1] lt x[p])and(x[p-1] gt  minh1[1,0]) do p -= 1 
lpff, x[(p-3):mh1v2],ppos  &  ppos=ppos + p-3
if (~finite(ppos)) then ppos = p-3. 
if (abs(ppos - p) gt 3.) then minh1[0,0] = p else minh1[0,0] = ppos
vk = interpol(reform(x[(mh1v1-20):mh1v2]), findgen(nk1v+20), findgen((nk1v+20.)*20.)/20.)
minh1[0,1] =  vk[(round((minh1[0,0] - mh1v1 + 20)  * 20.))] * konti
;-------------------------------------------------
cagravity, minh1[0,0], minh1[1,0], x, gravity1, amp1
if (abs(gravity1 - 106.) gt 10.) then gravity1 = 110 

if (do_H_line eq 1) then begin  ; for the Mg II h line
   cagravity, minh1[0,0], minh1[1,0], x, gravity1, amp1

    ; flatten the wing to allow better Gaussian fitting, h line
   u = [findgen(45/dfac)+gravity1 + 20/dfac, findgen(38/dfac)+gravity1 - 58/dfac]
endif  else begin
    ; flatten the wing to allow better Gaussian fitting, k line
   u = [findgen(25/dfac)+gravity1 - 60/dfac, findgen(25/dfac)+gravity1 + 25/dfac, findgen(5/dfac)+gravity1 - 25/dfac]
   u = [u, findgen(30/dfac)+gravity1 - 95/dfac, findgen(8/dfac)+gravity1 + 55/dfac]
endelse

if (konti gt 10.) then begin  ; skip this correction off disk
   iline3 = median(iline2, 3)
   res = poly_fit(u, iline3[u], 3, /double)
   xx = findgen(n_elements(line_profile))
   yfit = res[0] + xx * res[1] + xx^2 * res[2] + xx^3 * res[3]
   erg = where( finite(yfit,/NAN) , count)
   if (count gt 0)or(max(abs(yfit)) gt 1.0d4) then begin
      res = poly_fit(u, iline3[u], 2, /double)
      xx = findgen(n_elements(line_profile))
      yfit = res[0] + xx * res[1] + xx^2 * res[2]   
   endif   
   iline2 = (iline2 - yfit) > 0.
endif

subtracted_wing = iline * konti -  iline2
   ; if one adds this curve to the 'py' (below), it produces the 'prof'.
   ; py is the wing-subtracted profile for the Gaussian fitting while prof
   ; is the original profiel for the line analysis. The emission peaks
   ; measured by the two methods have a little offset due to this
   ; curve which can be compensated. 


;------------------------------------------------
;-- Mg II core parameters: profile analysis
;------------------------------------------------
corepos = gravity1 ;110     ; rough core position
core_pos = corepos
core0 = corepos
if (disper lt 3.)and(disper gt 2.) then int_area = prof(corepos-29:corepos+29)
if (disper gt 5.)and(disper lt 6.) then int_area = prof(corepos-14:corepos+14) ; Mg II core region
;int_area = smooth(int_area,3,/edge_truncate)
np = n_elements(int_area)


;------------------------------------------------
;-- cut-out the line rofile
;------------------------------------------------
if (disper lt 3.)and(disper gt 2.) then begin
   int_area2 = iline2[corepos-29:corepos+29]
   subtracted_wing = subtracted_wing[corepos-29:corepos+29]
   py = reform(int_area2[6:np-6])
   subtracted_wing = subtracted_wing[6:np-6]
endif

if (disper gt 5.)and(disper lt 6.) then begin
   int_area2 = iline2[corepos-14:corepos+14]
   subtracted_wing = subtracted_wing[corepos-14:corepos+14]
   py = reform(int_area2[3:np-3]) 
   subtracted_wing = subtracted_wing[3:np-3]
endif

;int_area2 = smooth(int_area2, 3, /edge_truncate)
;if (np gt  40) then  py = reform(int_area2[6:np-6]) else py = reform(int_area2[3:np-3]) ; the Gaussian fit range
;------------------------------------------------
;-- Mg II core parameters: Gaussian fit
;------------------------------------------------
nnp = n_elements(py)
fitsg = fltarr(nnp)
fitdg = fltarr(nnp)
fittg = fltarr(nnp)
spar1 = fltarr(9)
dpar1 = fltarr(15)
tpar1 = fltarr(21)

erg = {type:intarr(1), hpar:fltarr(6), spar:fltarr(9), dpar:fltarr(15),  tpar:fltarr(21), $
       fe_int:fltarr(n_elements(ca_ilo)), band:fltarr(30), eml:fltarr(12), ems:fltarr(4), $
       xmins:fltarr(4,2), xmaxs:fltarr(3,2), xmh1:fltarr(2,2), velpos:fltarr(n_elements(ca_ilo)), $
       sub_wing:fltarr(nnp), sfit:fltarr(nnp), sprf:fltarr(nnp), dfit:fltarr(nnp), tfit:fltarr(nnp)}

bbc = where(int_area gt (-1.0), bppc)

if (bppc le 0) then return, erg
;------------------------------------------------
;------------------------------------------------
;-- blend analysis around Mg II k line
;------------------------------------------------
;------------------------------------------------
if (disper lt 3.) then begin 
   u = iline[185:225]
   emb[0] = (mean(u[0:3]) > mean(u[37:*])) > 0.
   lpff, u[15:25], posm  & posm+=15.   &   if (~finite(posm)) then posm=20

   if (u[posm] lt emb[0]) then begin ; absorption line
      minimum, posm, u, tmp
      emb[1]=  tmp[4] > 0.
      embv[0] = tmp[3] + 185.
   endif else begin
      maximum, posm, u, tmp
      emb[1]=  tmp[4] > 0.
      embv[0] = tmp[3] + 185.
   endelse
   emb[2] = total(iline[185:225] - emb[0])* disper /konti ; observed equivalent width
   ;---------------------------
   u = iline[40:50]
   emb[3] = (mean(u[0:1]) > mean(u[8:9]))>0.
   lpff, u, posm    &   if (~finite(posm)) then posm=5
   emb[4] = iline[posm + 40.] >0.
   emb[5] = total(u - emb[3])* disper/konti ; observed equivalent width
   embv[1] = posm + 40.
   ;---------------------------
   u = iline[165:175]
   emb[6] = (mean(u[0:1]) > mean(u[8:9]))>0.
   lpff, u, posm    &   if (~finite(posm)) then posm=5
   emb[7] = iline[posm + 165.]>0.
   emb[8] = total(u - emb[6])* disper/konti ; observed equivalent width
   embv[2] = posm + 165.
   ;---------------------------
   if (np gt 263) then begin
      u = iline[252:262]
      emb[9] = (mean(u[0:1]) > mean(u[8:9]))>0.
      lpff, u, posm    &   if (~finite(posm)) then posm=5
      emb[10] = iline[posm + 252.] > 0.
      emb[11] = total(u - emb[9])* disper /konti ; observed equivalent width
      embv[3] = posm + 252.
   endif
   ;---------------------------
endif
if (disper gt 3.)and(disper lt 6.) then begin 
   u = iline[145:172]
   emb[0] = (mean(u[0:3]) > mean(u[27:*]))>0.
   lpff, u[12:21], posm  & posm+=12.   &   if (~finite(posm)) then posm=15

   if (u[posm] lt emb[0]) then begin ; absorption line
      minimum, 4., u[12:21], tmp
      emb[1]=  tmp[4] > 0.
      embv[0] = tmp[3] + 154.
   endif else begin
      maximum, 4., u[12:21], tmp
      emb[1]=  tmp[4] > 0.
      embv[0] = tmp[3] + 154.
   endelse
   emb[2] = total(iline[145:172] - emb[0])* disper /konti ; observed equivalent width
endif

;---------------------------

if (plt eq 1) and (do_gauss eq 0) then begin
   plot, prof, xstyle=1
   if (disper lt 3.) then begin  ;------------------- 2.544 pm/px
      oplot, [corepos-24,corepos-24],[0,8]
      oplot, [corepos+25,corepos+25],[0,8]
      vline, emb[1]
      vline, emb[3]
      vline, emb[5]
      vline, emb[7]
   endif
   if (disper gt 4.) then begin  ;------------------- 5.088 pm/px
      oplot, [corepos-14,corepos-14],[0,8]
      oplot, [corepos+14,corepos+14],[0,8]
   endif
   
endif


lobes = fltarr(2)           
i_peaks = fltarr(2)
i_core = 0.
type = 0       ; profile type

xmax = fltarr(8)
xmin = fltarr(8)
cmax = 0
cmin = 0
if (disper lt 3.) then ux = corepos-29.  else ux = corepos-14.


if (disper gt 4.) then begin
   i = 3
   max_loc = 6
   check_limit = 25
endif else begin
   i = 5
   max_loc = 12
   check_limit = 50
endelse
;###################################################################################
while (i le check_limit)and(cmin lt 8)and(cmax lt 8) do begin
  y1 = int_area[i]
  if  (cmax ge 2) then m1 = max(int_area[i-3:i+3]) else m1 = max(int_area[i-1:i+2])
  if  (cmax ge 1)and(cmin eq 1) then m2 = min(int_area[i-2:((i+1)< check_limit)]) else m2 = min(int_area[i-2:((i+2)< check_limit)])
  if (i eq 5) then m2 = int_area[i]

  if ((y1) ge m1)AND(i lt (check_limit - 6))and(i gt 9)and(y1 ge threshold) then begin
	xmax[cmax] = i
	cmax  = cmax + 1
        ;i = i + 1
        if (cmax eq 1) then i += 1
        if (cmax eq 0)and(cmin eq 1) then i += 3
        if (cmax eq 3)or(cmin eq 3) then i += 3
        if (cmax eq 1)and(xmax[0] lt max_loc) then i += 1
  endif 
  if (y1 eq m2) then  begin 
	xmin[cmin] = i
	cmin  = cmin + 1
        i = i + 1
  endif 
  i = i + 1
endwhile
;###################################################################################

if ((xmin[1]-xmin[0]) le 4.) then begin
  cmin -= 1
  xmin = xmin[1:*]
endif
if (cmax eq 3) then begin
  if (int_area[xmax[0]] lt int_area[xmax[1]])and(int_area[xmax[0]] lt int_area[xmax[2]]) then begin
    cmax = 2
    xmax = xmax[1:2]
 endif
endif
if (cmax eq 3) then begin
  if (int_area[xmax[2]] lt int_area[xmax[1]])and(int_area[xmax[2]] lt int_area[xmax[0]]) then begin
    cmax = 2
    xmax = xmax[0:1]
  endif 
endif 

if (cmax gt 0) then xmax = xmax(0:cmax-1) else xmax = fltarr(1)
if (cmin gt 0) then xmin = xmin(0:cmin-1) else xmin = fltarr(1)

if (cmin eq 5) then begin
  if (abs(xmin[4] - xmin[3]) lt 9)and(xmin[3] ge 39.) then begin
    cmin -=1
    xmin = xmin[0:3]
  endif
endif 
if (cmin ge 5) then begin
  if (abs(xmin[1] - xmin[0]) lt 5)and(xmin[1] lt 10) then begin
    cmin -=1
    xmin = xmin[1:*]
  endif
endif 
if (cmin eq 5) then begin
  if (abs(xmin[3] - xmin[2]) lt 5)and(xmin[2] ge 40.) then begin
    cmin -=2
    xmin = xmin[0:3]
  endif
endif 
if (cmin eq 4) then begin
  if ((abs(xmin[3] - xmin[2]) lt 9.)and(xmin[2] gt 30.)) then begin
    cmin -=1
    xmin = xmin[0:2]
  endif 
endif 

if (cmin eq 4) then begin
  if ((abs(xmin[3] - xmin[2]) lt 9.)and(xmin[2] gt 40.)) then begin
    cmin -=1
    xmin = xmin[0:2]
  endif 
endif 
if (cmin eq 4) then begin
  if ((abs(xmin[3] - xmin[2]) lt 19.)and(xmin[2] gt 35.))and(max(int_area[xmin[2]:xmin[3]]) lt 0.7) then begin
    cmin -=1
    xmin = xmin[0:2]
  endif 
endif 
if (cmin eq 4)and(cmax eq 1) then begin
  if ((abs(xmin[2] - xmin[1]) lt 9.)and(xmin[1] ge 39.)) then begin
    cmin -=1
    xmin = [xmin[0], xmin[1], xmin[3]]
  endif 
endif 
if (cmin eq 4)and(cmax eq 2) then begin
  if ((xmin[2]-xmin[1]) eq 1.) then begin
    cmin -=1
    xmin = [xmin[0], xmin[1], xmin[3]]
  endif 
endif 
if (cmin eq 4)and(cmax eq 1) then begin
  if ((xmin[1]-xmin[0]) le 11.)and(xmin[1] lt 12) then begin
    cmin -=1
    xmin = xmin[1:*]
  endif 
endif 
if (cmin eq 4)and(cmax le 2) then begin
  if ((abs(xmin[3] - xmin[2]) le 10.)and(xmin[2] ge 39.)) then begin
    cmin -=1
    xmin = [xmin[0], xmin[1], xmin[2]]
  endif 
endif 


if (cmin le 3)and(cmax gt 3) then begin  ; in case of redundant maxima
  if ((xmax[2]-xmax[1]) le 4.)and((xmax[3]-xmax[2]) le 4.) then begin
    cmax = 2
    xmax = [xmax[0], xmax[1]]

    if (cmin eq 2) then begin 
       if (xmin[0] lt xmax[0])and(xmin[1] gt xmax[1]) then begin
         cmin = 3
         newmin = min(int_area[xmax[0]:xmax[1]], nnewpos)
         xmin =[xmin[0], nnewpos+xmax[0], xmin[1]]
       endif
    endif

    if (cmin eq 1)and(xmin[0] lt xmax[0]) then begin ; only first minimum
       cmin = 3
       newmin = min(int_area[xmax[0]:xmax[1]], nnewpos)
       xmin =[xmin[0], nnewpos+xmax[0]]
       newmin = min(int_area[xmax[1]:xmax[1]+20], nnewpos)
       xmin =[xmin[0], xmin[1], nnewpos+xmax[1]]
    endif
    if (cmin eq 1)and(xmin[0] gt xmax[1]) then begin; only third minimum
       cmin = 3
       newmin = min(int_area[(xmax[0]-9)>0:xmax[0]], nnewpos)
       xmin =[xmin[0], nnewpos+xmax[0]-9.]
       newmin = min(int_area[xmax[0]:xmax[1]], nnewpos)
       xmin =[xmin[0], xmin[1], nnewpos+xmax[0]]
    endif

  endif 
endif 

if (n_elements(xmin) gt 0)and(n_elements(xmin) ne cmin) then cmin = n_elements(xmin)
if (n_elements(xmax) gt 0)and(n_elements(xmax) ne cmax) then cmax = n_elements(xmax)


if (disper lt 3.) then begin ;###############################################
if (cmin eq 3)and(cmax eq 1) then begin
  if ((abs(xmin[2] - xmin[1]) lt 9.)and(xmin[1] gt 40.)) then begin
    cmin -=1
    xmin = xmin[0:1]
  endif 
endif 

if (cmin eq 2)and(cmax eq 3) then begin
  if (xmax[2] gt 40.)and(int_area[xmax[2]] lt 1.) then begin
    cmax -=1
    xmax = xmax[0:1]
    cmin += 1
    xn = round(total(float(xmax))*0.5)
    xmin = [xmin[0], xn, xmin[1]]
  endif 
endif 

if (cmin ge 2) then begin
if (xmin[0] le 10.) and (xmin[1] le 12.) then begin
  cmin -= 1
  xmin = xmin[1:*]
endif
endif

if (cmin eq 1)and(cmax ge 1) then begin
  if ((xmin[0] lt 12.)) then begin
    cmin +=1
    xn = min(int_area[38:*],pos)
    xmin = [xmin[0], pos+38.]
  endif 
endif 

if (cmin eq 0)and(cmax ge 1) then begin
    cmin = 2
    xn = min(int_area[4:12],pos)
    xmin = [pos+4.]
    xn = min(int_area[38:*],ppos)
    xmin = [xmin, ppos+38.]
endif 

if (cmin eq 2)and(cmax eq 2) then begin
  if ((int_area[xmax[0]] gt 0.7)and(int_area[xmax[1]] gt 0.7)) then begin
    cmin +=1
    xn = min(int_area[xmax[0]:xmax[1]], pos)
    xmin = [xmin[0], pos+xmax[0], xmin[1]]
  endif 
endif 

endif else begin  ;##########################################
if (cmin eq 3)and(cmax eq 1) then begin
  if ((abs(xmin[2] - xmin[1]) lt 5.)and(xmin[1] gt 20.)) then begin
    cmin -=1
    xmin = xmin[0:1]
  endif 
endif 


if (cmin eq 2)and(cmax eq 3) then begin
  if (xmax[2] gt 20.)and(int_area[xmax[2]] lt 1.) then begin
    cmax -=1
    xmax = xmax[0:1]
    cmin += 1
    xn = round(total(float(xmax))*0.5)
    xmin = [xmin[0], xn, xmin[1]]
  endif 
endif 

if (cmin ge 2) then begin
if (xmin[0] le 5.) and (xmin[1] le 6.) then begin
  cmin -= 1
  xmin = xmin[1:*]
endif
endif

if (cmin eq 1)and(cmax ge 1) then begin
  if ((xmin[0] lt 6.)) then begin
    cmin +=1
    xn = min(int_area[19:*],pos)
    xmin = [xmin[0], pos+19.]
  endif 
endif 

if (cmin eq 0)and(cmax ge 1) then begin
    cmin = 2
    xn = min(int_area[2:6],pos)
    xmin = [pos+2.]
    xn = min(int_area[19:*],ppos)
    xmin = [xmin, ppos+19.]
endif 

if (cmin eq 2)and(cmax eq 2) then begin ;-- either umbral profile or plage
  if ((int_area[xmax[0]] gt 0.7)and(int_area[xmax[1]] gt 0.7)) then begin
    cmin +=1
    xn = min(int_area[xmax[0]:xmax[1]], pos)
    xmin = [xmin[0], pos+xmax[0], xmin[1]]
  endif 
endif

if (cmax le 1) then begin
    cmax =2
    xn = max(int_area[7:14], pos1)
    xn = max(int_area[17:20], pos2)
    xmax = [pos1+7, pos2+17]
    if ((xmax[1]-xmax[0]) lt 2.) then begin
       xn = max(int_area, pos)
       if (int_area[pos-1] lt int_area[pos])and(int_area[pos+1] lt int_area[pos]) then begin
         xmax = pos
         cmax = 1
      endif else begin
         cmax = 0
         xmax = 0
      endelse   
    endif
endif   
endelse  ;################################################################   

if (n_elements(xmin) gt 0)and(n_elements(xmin) ne cmin) then cmin = n_elements(xmin)
if (n_elements(xmax) gt 0)and(n_elements(xmax) ne cmax) then cmax = n_elements(xmax)

if (cmin eq 2)  then begin
 if (disper lt 3.) then begin  ;----------------- dispersion = 2.5 pm
  if (xmin[1] le 30.) then begin
    if (xmin[1] le 11.) then begin
      xmin = xmin[1]
      xn = min(int_area[30:50], pos)
      xmin = [xmin, pos+30.]
   endif else begin
      xmin = xmin[0]
      xn = min(int_area[30:50], pos)
      xmin = [xmin, pos+30.]
  endelse
  endif else begin
    if (xmin[0] gt 11.) then begin
      xmin = xmin[1]
      xn = min(int_area[3:12], pos)
      xmin = [xmin, pos+3.]
   endif 
  endelse
 endif else begin;----------------- dispersion = 5 pm
  if (xmin[1] le 19.) then begin
    if (xmin[1] le 7.) then begin
      xmin = xmin[1]
      xn = min(int_area[15:25], pos)
      xmin = [xmin, pos+15.]
   endif else begin;--- missing K1r minima
      xmin = xmin[0]
      xn = min(int_area[17:*], pos)
      xmin = [xmin, pos+17.]
      cmin += 1
  endelse
  endif else begin
    if (xmin[0] gt 7.) then begin;--- missing K1v minima
      ;xmin = xmin[1]
      xn = min(int_area[1:6], pos)
      xmin = [pos+1., xmin]
      cmin += 1
    endif
    if (xmin[0] le 7.)and(xmin[1] ge 19.)and(cmax eq 2) then begin;--- missing central minima
       xn = min(int_area[xmax[0]:xmax[1]], pos) &  nm = pos
       if ((xmax[1]-nm) lt 2.) then begin
          xmax[1] += 1
          nm -= 1
       endif
       if ((nm - xmax[0]) lt 2.) then begin
          xmax[1] -= 1
          nm += 1
       endif
       xmin = [pos, xmin]
       cmin += 1
    endif
    
  endelse   
 endelse  
endif

if (n_elements(xmin) gt 0)and(n_elements(xmin) ne cmin) then cmin = n_elements(xmin)
if (n_elements(xmax) gt 0)and(n_elements(xmax) ne cmax) then cmax = n_elements(xmax)
;---------------------------
;-- omit extra minimum
;---------------------------
if (cmin ge 4)and(cmax ge 2) then begin
  if (xmin[1] lt xmax[0])then begin 
    xmin = xmin[1:3]
    cmin -= 1
  endif
  if (xmin[2] gt (xmax[1]+2)) then begin
    xmin = xmin[0:2]
    cmin = 3
  endif
endif
if (cmin ge 4)and(cmax eq 1) then begin
  if (xmin[1] lt xmax[0])then begin 
    xmin = xmin[1:3]
    cmin -= 1
  endif
  if (xmin[1] lt (xmax[0])) then begin
    xmin = xmin[1:*]
    cmin -=1
  endif
endif

if (plt gt 1) then print, cmin, cmax

;---------------------------
;-- omit extra maximum
;---------------------------
if (cmin eq 3)and(cmax eq 4) then begin
  if (xmin[2] lt xmax[3])then begin 
    xmax = xmax[0:2]
    cmax -= 1
  endif
  if (xmin[0] gt xmax[0]) then begin
    xmin = xmin[1:*]
    cmax -= 1
  endif
  if (cmax eq 4) then begin
    u = xmax
    v = int_area[xmax]
    vv = v(sort(v))
    uu = u(sort(v))
    xmax = uu[2:3]
    cmax = 2
  endif

endif



if (cmin eq 3)and(cmax eq 2) then begin
   if (xmax[0] gt xmin[1]) then begin
     xn = min(int_area[0:xmax[0]], pos1)
     xn = min(int_area[xmax[1]:*], pos2)  &  pos2 += xmax[1]
     xmin = [pos1, xmin[1], pos2]
  endif 
endif
if (plt gt 1) then print, 'n_min, n_max: ', cmin, cmax ;, xmax
;=========================================================

if (cmin gt 1.) then xmin = xmin(sort(xmin))
if (xmin[0] eq 5)and(int_area[4] lt int_area[5]) then xmin[0]=4
if (xmin[0] eq 4)and(int_area[3] lt int_area[4]) then xmin[0]=3

if (cmin eq 3) then begin
  if (xmin[1] le 5.) then begin
     xmin = [xmin[0], xmin[2]]
     cmin = 2
  endif
endif

if (cmin eq 2) then begin
  if (n_elements(int_area) gt 50)and(xmin[1] lt 40)and(xmin[1] gt 20) then begin
     xn = min(int_area[xmin[1]:*], pos)
     xmin = [xmin, pos + xmin[1]]
     cmin +=1
  endif
endif

;--------------------------------------------------------------------------------------------
;    >>>>>>>>>> at this level, we add the fix offset to the line positions <<<<<<<<<<<<<
;++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
line_offset = 0
if (disper lt 3.)and(disper gt 2.) then begin
 line_offset = corepos - 29
endif   
if (disper gt 5.) then begin
 line_offset = corepos - 14
endif   
xmax = xmax  + line_offset
xmin = xmin  + line_offset


if (plt gt 1) then print, '-----------------------------------'
if (plt gt 1) then print, 'ff ', cmin, cmax

;-----------------------------------
;-- if it failed to find the core 
;-----------------------------------
t_check = 18.
if (disper gt 4.) then t_check = 9.
if (cmin eq 0) then begin
  if (cmax eq 2) then begin                ; either there is a minimum
    if ((xmax[1]-xmax[0]) gt 5) then begin;-----------------------------
      s0 = min(prof[xmax[0]:xmax[1]], s)                          ;- 
      s  = s + xmax[0]                                                ;-
      lpff, prof(s-2:s+2),pos_core                                ;-

       if (pos_core le xmax[0])OR(pos_core ge xmax[1]) then begin     ;-
         if ((xmax[1]-xmax[0])ge t_check) then begin                      ;-
         pos_core = s                                                 ;-
         i_core = prof(s)                                         ;- 
         core_pos =  s                                                ;- 
         cmin = 1                                                     ;-
        endif                                                         ;-
      endif else begin                                               ;-
        pos_core = (xmax[0]+xmax[1])*.5                              ;-
        i_core = prof(pos_core + s - 2)                          ;-
        core_pos =  pos_core + s - 2                                 ;-
        cmin = 1                                                     ;-
      endelse                                                        ;-
    endif else begin                                                 ;-
      s = round((xmax[0] + xmax[1]) * 0.5)                           ;-
      i_core = prof(s)                                           ;-
      core_pos = s                                                   ;-
      cmin = 0    &  cmax = 1                                        ;-
    endelse   ;--------------------------------------------------------
  endif    

  ; or we have just one maximum like umbral profiles
  if (cmax eq 1)AND((core0 - xmax[0]) le 6.) then begin
    i_core = i_peaks[0]
    core_pos = lobes[0]
    cmin = 0    &  cmax = 1
  endif
endif
if (n_elements(xmin) gt 0)and(n_elements(xmin) ne cmin) then cmin = n_elements(xmin)
if (n_elements(xmax) gt 0)and(n_elements(xmax) ne cmax) then cmax = n_elements(xmax)
;------------------------------------ 
;-find max/min position/amplitude
;-
;-and store them in xmaxs and xmins 
;------------------------------------
if (plt gt 1.)then print, 'aa ', cmin, cmax, xmax
xmin = xmin(sort(xmin))
xmax = xmax(sort(xmax))
q=where(xmin lt ux, count) & if (count gt 0) then xmin(q) += ux
q=where(xmax lt ux, count) & if (count gt 0) then xmax(q) += ux


if (cmin gt 1) then xmins = fltarr(cmin,2) else xmins = fltarr(1,2)
if (cmax gt 1) then xmaxs = fltarr(cmax,2) else xmaxs = fltarr(1,2)

;------------------------------------
;--- for minimums
;------------------------------------
if (cmin ge 1) then begin
  for i=0, cmin-1 do begin
    tmp =  xmin[i]
    lpff, prof(tmp-2:tmp+2),pos_core
    if (~Finite(pos_core)) then pos_core = 0.
  pos_core += tmp - 2.
  tmp = round(pos_core(0)) + findgen(5) - 2
  linecut = prof(tmp)
  res = poly_fit(tmp, linecut, 2, /double)
  xmins[i,0] = pos_core(0)
  xmins[i,1] = res[0] + res[1]*pos_core[0] +  res[2]*(pos_core[0]^2)

  if (abs(xmins[i,0] - xmin[i]) ge 1.) then begin 
    tmp =  xmin[i]
    lpff, prof(tmp-3:tmp+3),pos_core
    if (~Finite(pos_core)) then pos_core = 0.
    pos_core += tmp - 3.
    tmp = round(pos_core(0)) + findgen(7) - 3
    linecut = prof(tmp)
    res = poly_fit(tmp, linecut, 2, /double)
    xmins[i,0] = pos_core(0)
    xmins[i,1] = res[0] + res[1]*pos_core[0] +  res[2]*(pos_core[0]^2)
  endif
  endfor            
endif
if (n_elements(xmin) gt 0)and(n_elements(xmin) ne cmin) then cmin = n_elements(xmin)
if (n_elements(xmax) gt 0)and(n_elements(xmax) ne cmax) then cmax = n_elements(xmax)
;------------------------------------
if (plt gt 1.)then print, 'ee ', cmin, cmax, xmaxs
xmin = xmin(sort(xmin))
xmax = xmax(sort(xmax))
;------------------------------------
;---- and maximums
;------------------------------------
if (cmax ge 1) then begin  ; emission peak(s)
  for i=0, cmax-1 do begin
    tmp =  xmax[i]
    lpff, prof(tmp-2:tmp+2),pos_core
    if (~Finite(pos_core)) then pos_core = 0.0
  pos_core += tmp - 2.
  tmp = round(pos_core(0)) + findgen(5) - 2
  linecut = prof(tmp)
  res = poly_fit(tmp, linecut, 2, /double)
  xmaxs[i,0] = pos_core(0)
  xmaxs[i,1] = res[0] + res[1]*pos_core[0] +  res[2]*(pos_core[0]^2)

  if (abs(xmaxs[i,0] - xmax[i]) ge 1.) then begin 
    tmp =  xmax[i]
    lpff, prof(tmp-4:tmp+4),pos_core
    if (~finite(pos_core)) then pos_core = core_pos 
    pos_core += tmp - 4.
    tmp = round(pos_core(0)) + findgen(9) - 4
    linecut = prof(tmp)
    res = poly_fit(tmp, linecut, 2, /double)
    xmaxs[i,0] = pos_core(0)
    xmaxs[i,1] = res[0] + res[1]*pos_core[0] +  res[2]*(pos_core[0]^2)
  endif
  endfor            
endif
if (plt gt 1.)then print, 'cc', cmin, cmax, xmax

cmin_cmax, xmins, xmaxs, cmin, cmax, xmin, xmax

;####################################################################################
;-------------------------------------
;-- to check for very weak max/mins
;-------------------------------------
crit = 0.0013 ; 0.0013
if (cmin eq 3 and cmax gt 2) then begin
  if ((xmaxs[0,1]-xmins[0,1]) lt crit) then begin 
   if (plt gt 1.)then print, 'aaa', xmaxs[0,1]-xmins[0,1]
     cmax = 1
     xmaxs = xmaxs[1,*]
     xmax = xmax[1]
 endif else begin
  if ((xmaxs[1,1]-xmins[2,1]) lt crit) then begin 
    if (plt gt 1.)then print, 'bbb', xmaxs[1,1]-xmins[2,1]
     cmax = 1
     xmaxs = xmaxs[0,*]
     xmax = xmax[0]
 endif
 endelse
endif
if (plt gt 1.)then print, 'ddd ', cmin, cmax, xmax
xmin = xmin(sort(xmin))
xmax = xmax(sort(xmax))
;------------------------------
if (cmin eq 3 and cmax gt 2) then begin ; in case of extra maximum
  ; print, xmaxs
  if ((xmaxs[0,1]-xmins[1,1]) lt (-0.5)) then begin 
     cmax -= 1
     xmaxs = xmaxs[1:*,*]
     xmax = xmax[1:*]
  endif 
endif
;------------------------------
if (cmin eq 3 and cmax gt 2) then begin ; in case of extra maximum
  if ((xmaxs[2,1]-xmins[1,1]) lt (-0.5)) then begin 
     cmax -= 1
     xmaxs = xmaxs[0:1,*]
     xmax = xmax[0:1]
  endif 
endif
;------------------------------
if (cmin eq 3)and(cmax gt 2)and(xmax[0] lt xmin[0]) then begin ; in case of wrong maximum
     cmax -= 1
     xmaxs = xmaxs[1,*]
     xmax = xmax[1]
endif
;------------------------------
cmin_cmax, xmins, xmaxs, cmin, cmax, xmin, xmax
;####################################################################################

if (cmin eq 3) and (cmax eq 1) then begin ; in case of plage profiles
   if (abs(xmins[1,0]-xmaxs[0,0]) lt abs(xmins[0,0]-xmaxs[0,0]))and(abs(xmins[1,0]-xmaxs[0,0]) lt abs(xmins[2,0]-xmaxs[0,0])) then begin
     cmax += 1
     if (xmins[1,0] gt xmaxs[0,0]) then begin ; only H2v
        xmins[1,0] -= 1  &  xmins[1,1] = prof(xmins[1,0])
        xn =xmins[1,0] + abs(xmins[1,0]-xmaxs[0,0])
        xmax = [xmax, xn]
        xmaxs = fltarr(2,2)  & xmaxs[*,0] = xmax &  xmaxs[*,1] = prof(xmax)
     endif else begin ; only H2r
        xmins[1,0] += 1  &  xmins[1,1] = prof(xmins[1,0])
        xn =xmins[1,0] - abs(xmins[1,0]-xmaxs[0,0])
        if (xn gt (xmaxs[0,0]+8.)) then begin
           xmax = [xn, xmax]
           xmaxs = fltarr(2,2)  & xmaxs[*,0] = xmax &  xmaxs[*,1] = prof(xmax)
        endif   else begin
           cmax -=1
        endelse
        saba = max(prof[xmin[0]:xmin[1]], posmax)
        saba = reform(prof[xmin[0]:xmin[1]])
        lpff, saba, pos1
        pos1 += xmin[0]
        if (finite(pos1))and(pos1 lt (n_elements(prof)-10)) then begin
           jk = min(saba, xk)
           maximum, xk, saba, ergebnis
           xmax = [ergebnis[3] + xmin[0], xmax[0]]   &   xmax = xmax(sort(xmax))
           xmaxs = fltarr(2,2)
           xmaxs[*,0] = xmax
           xmaxs[*,1] = [prof[xmax]]
           cmax = 2
        endif
     endelse
  endif else begin  ;--- if not, then remove extra minimum

     cmin -=1
     xmin = [xmin[0], xmin[2]]
     u = [xmins[0,0], xmins[2,0]]
     v = [xmins[0,1], xmins[2,1]]
     xmins = fltarr(2,2)  & xmins[*,0] = u &  xmins[*,1] = v
  endelse   
endif
if (plt gt 1.)then print, 'km ', cmin, cmax

if (cmax eq 2) then begin
  if(xmaxs[0,0] ge (n_elements(prof)-10))or(xmaxs[0,0] le 50) then begin
    cmax -=1
    xmaxs = xmaxs[1,*]
    xmax = xmax[1]
  endif
  if(xmaxs[1,0] ge (n_elements(prof)-10))or(xmaxs[1,0] le 50) then begin
    cmax -=1
    xmaxs = xmaxs[0,*]
    xmax = xmax[0]
  endif
endif

cmin_cmax, xmins, xmaxs, cmin, cmax, xmin, xmax
cmax = n_elements(xmaxs)/2
cmin = n_elements(xmins)/2
if (plt gt 1.)then print, 'kn ', cmin, cmax
;####################################################################################
;-------------------------------------
;------    some final checks
;-------------------------------------

if (abs(xmax[0] - xmaxs[0,0]) gt 1.) then begin 
   xmaxs[0,0] = xmax[0] 
   xmaxs[0,1] = prof(xmax[0])
endif
if (xmax[0] lt xmin[0])and(cmax eq 1)and(cmin ge 2) then begin 
        saba = max(prof[xmin[0]:xmin[1]], posmax)
        saba = reform(prof[xmin[0]:xmin[1]])
        lpff, saba, pos1
        pos1 += xmin[0]
        if (finite(pos1))and(pos1 lt (n_elements(prof)-10)) then begin
           jk = min(saba, xk)
           maximum, xk, saba, ergebnis
           xmax = ergebnis[3] + xmin[0]
           xmaxs[0,0] = xmax
           xmaxs[0,1] = prof[xmax]
           cmax = 1
        endif   
endif

if (cmin gt 0) then begin ; fake maximum ?
   for f=0, cmin-1 do if (xmins[f, 1] lt 0.) then xmins[f, *] = 0
endif    
if (cmax gt 0) then begin ; fake maximum ?
   for f=0, cmax-1 do if (xmaxs[f, 1] lt 0.) then xmaxs[f, *] = 0
endif    

if (cmax gt 1) then begin 
  if (abs(xmax[1] - xmaxs[1,0]) gt 1.) then begin 
    xmaxs[1,0] = xmax[1] 
    xmaxs[1,1] = prof(xmax[1])
  endif
endif

if (cmin gt 1) then begin 
  if (abs(xmin[1] - xmins[1,0]) gt 1.) then begin 
    xmins[1,0] = xmin[1] 
    xmins[1,1] = prof(xmin[1])
  endif
endif
cmin_cmax, xmins, xmaxs, cmin, cmax, xmin, xmax
if (plt gt 1.)then print, 'kb ', cmin, cmax
if (cmin gt 4) then begin
   xmins = xmins[0:3,*]
   cmin_cmax, xmins, xmaxs, cmin, cmax, xmin, xmax
endif   
if (cmax gt 3) then begin
   xmaxs = xmaxs[0:2,*]
   cmin_cmax, xmins, xmaxs, cmin, cmax, xmin, xmax
endif   
;####################################################################################
;----------------------
;-- find profile type
;----------------------
type = -1
if (cmax eq 0) then type = 0 ; absorption profile
if (cmax eq 2)AND(cmin eq 3) then type = 1 ; normal profile
if (cmax eq 3)AND(cmin ge 3) then type = 2 ; 3-lobe profile
if (cmax eq 1)AND(cmin eq 3) then type = 3 ; plage profile
if (cmax eq 1)AND(cmin eq 1) then type = 4 ; very asymm. profile
if (cmax eq 1)AND(cmin eq 2) then type = 5 ; umbral profile
;-------------------------------------
;---- improve the accuracy of H1/2/3
;-------------------------------------
if (type eq 1)or(type eq 2) then begin
  
  if (prof[xmax[0]-1] gt prof[xmax[0]]) then xmax[0]-=1
  mh=5
  mksp = prof[xmax[0]-mh:xmax[0]+2] 
  maximum, mh, mksp, ergebnis
  xmax[0]=ergebnis[3] + xmax[0]- float(mh)

  if (prof[xmax[1]-1] gt prof[xmax[1]]) then xmax[1]-=1  &  mh = 2
  mksp = prof[xmax[1]-mh:xmax[1]+5]
  maximum, mh, mksp, ergebnis
  xmax[1]=ergebnis[3] + xmax[1]- float(mh)

  mksp = prof[xmin[1]-3:xmin[1]+3]  &  minimum, 3, mksp, ergebnis
  if (abs(ergebnis[3]-3.) lt 2.) then xmin[1]=ergebnis[3] + xmin[1]-3.   ; if the minimum is reasonalbe

endif
if (type eq 3) then begin
  if (prof[xmax-1] gt prof[xmax]) then xmax-=1
  mksp = prof[xmax-5:xmax+4]  &  maximum, 5, mksp, ergebnis
  xmax=ergebnis[3] + xmax-5.

  mksp = prof[xmin[1]-4:xmin[1]+4]  &  minimum, 4, mksp, ergebnis
  if (abs(ergebnis[3]-4.) lt 2.) then xmin[1]=ergebnis[3] + xmin[1]-4.   ; if the minimum is reasonalbe

endif
if (type eq 5)or(type eq 4) then begin
   if (xmax gt (n_elements(prof) - 3)) then begin
      temp = max(prof, xmax)
  endif   
  if (prof[xmax-1] gt prof[xmax]) then xmax-=1
  mksp = prof[xmax-5:xmax+4]  &  maximum, 5, mksp, ergebnis
  xmax=ergebnis[3] + xmax-5.
endif
;####################################################################################
;-------------------------------------
;-- plot lines, positions of peaks
;-------------------------------------
if (plt ge 1) and (do_gauss eq 0) then begin
  plot,int_area, psym=-1, /xstyle;,xrange=[60,140]
  print, cmin, cmax
  print, xmin, xmax
  print, '----------------------------------------------------'
  print,  abs(xmins[*,0] - xmin)
  print,  abs(xmaxs[*,0] - xmax)
  tm = max(int_area)*1.2
  loadct,40,/silent
  if (disper lt 3.) then begin
     for i=0,cmin-1 do oplot,[xmin[i]-corepos+29,xmin[i]-corepos+29],[0,tm],color=70
     for i=0,cmin-1 do oplot,[xmins[i,0]-corepos+29,xmins[i,0]-corepos+29],[0,tm], linestyle=1,color=70
     for i=0,cmax-1 do oplot,[xmax[i]-corepos+29,xmax[i]-corepos+29],[0,tm],linestyle=2, color=245
     for i=0,cmax-1 do oplot,[xmaxs[i,0]-corepos+29,xmaxs[i,0]-corepos+29],[0,tm],linestyle=1, color=245
  endif
  if (disper gt 4.) then begin
     for i=0,cmin-1 do oplot,[xmin[i]-corepos+14,xmin[i]-corepos+14],[0,tm],color=70
     for i=0,cmin-1 do oplot,[xmins[i,0]-corepos+14,xmins[i,0]-corepos+14],[0,tm], linestyle=1,color=70
     for i=0,cmax-1 do oplot,[xmax[i]-corepos+14,xmax[i]-corepos+14],[0,tm],linestyle=2, color=245
     for i=0,cmax-1 do oplot,[xmaxs[i,0]-corepos+14,xmaxs[i,0]-corepos+14],[0,tm],linestyle=1, color=245
  endif
  loadct,0, /silent
  if (cmax eq 8) then stop
  print, 'Press enter to continue ...'
  ans = get_kbrd()
endif

;####################################################################################
;--------------------------------------------
;-- only for umbral profiles
;--------------------------------------------
if type eq 5 then begin

  ew1 = 0.
  sh3 =  xmins[0,1]>xmins[1,1]
  hpar[0] = (xmins[1,0] - xmins[0,0]) < (254./disper) ; H1 emission width
  lobe_width_ca, 0, prof, xmaxs[0,0], xmaxs[0,1], sh3, fwhm1, ew1,/mg
  hpar[1] = fwhm1 < (38./disper)                   ; H2 emission width

  xwb1 = (xmins[0,0]+xmaxs[0,0])/2.
  tmp = ceil(xwb1)
  ywb1 = (prof[tmp]-prof[tmp-1])*(xwb1 - tmp*1.+1.)+ prof[tmp-1]
 
  xwb2 = (xmins[1,0]+xmaxs[0,0])/2.
  tmp = ceil(xwb2)
  ywb2 = (prof[tmp]-prof[tmp-1])*(xwb2 - tmp*1.+1.)+ prof[tmp-1]

  hpar[2] = xwb2 - xwb1       ; WB emission width
  hpar[3] = (ywb2 + ywb1)*0.5 ; mean WB intensity
endif
;####################################################################################
;--------------------------------------------
;-- only for normal profiles
;--------------------------------------------
if type eq 1 then begin
  hpar[0] = (xmins[2,0] - xmins[0,0]) < (254./disper) ; H1 emission width
  hpar[1] = (xmaxs[1,0] - xmaxs[0,0]) < (254./disper) ; H2 emission width

  if (xmaxs[0,0] gt n_elements(prof))or(abs(xmaxs[0,0] - xmax[0]) gt 10.) then begin
     xmaxs[0,0] = corepos  < ((xmins[0,0] + xmins[1,0])*0.5)
     xmaxs[0,1] = prof(xmaxs[0,0])
  endif
  if (xmaxs[1,0] lt xmins[1,0])or(xmaxs[1,0] gt n_elements(prof)) then begin
     xmaxs[1,0] = (xmins[2,0] + xmins[1,0])*0.5
     xmaxs[1,1] = prof(xmaxs[1,0])
  endif

  if (prof(round(xmins[1,0]-1)) lt prof(round(xmins[1,0]))) then begin
     xmins[1,0] -=1
     xmins[1,1] = prof(xmins[1,0])
  endif

  if (xmaxs[0,0] lt 20.) then stop

  xwb1 = (xmins[0,0]+xmaxs[0,0])/2.
  tmp = ceil(xwb1)
  ywb1 = (prof[tmp]-prof[tmp-1])*(xwb1 - tmp*1.+1.)+ prof[tmp-1]
 
  xwb2 = (xmins[2,0]+xmaxs[1,0])/2.
  tmp = ceil(xwb2)
  ywb2 = (prof[tmp]-prof[tmp-1])*(xwb2 - tmp*1.+1.)+ prof[tmp-1]

  hpar[2] = xwb2 - xwb1       ; WB emission width
  hpar[3] = (ywb2 + ywb1)*0.5 ; mean WB intensity

  hpar[4] = xmaxs[0,1]/xmaxs[1,1]   ; H2v/H2r
  hpar[5] = (xmaxs[0,1]/xmins[1,1]) < 10.   ; H2v/H3
endif
;####################################################################################



;########################################################################
;--  band intensities, similar to Ca II line
;########################################################################
if (disper ge 5.)and(disper lt 6.) then begin  ;-------------------------- 5.088 pm/px

  band[0] = total(prof[120:125])    ; ~ K1r
  band[1] = total(prof[85:94])      ; K-line blue wing
  band[2] = total(prof[65:78])      ; K-line blue wing
  band[3] = total(prof[145:154])    ; K-line red wing
  band[4] = total(prof[156:160])    ; Mg II 2798  --- K-line red wing
  band[5] = total(prof[140:149])    ; W4  --- K-line red wing

  band[6] = total(prof[110:112])    ; K3  3 pixel
  band[7] = total(prof[104:109])    ; K2v 6 pixel
  band[8] = total(prof[113:118])    ; K2r 6 pixel

  band[9] = total(prof[102:122])    ; K-index 1.0 A filter
  band[10] = total(prof[107:117])   ; K-index 0.506 A filter
  band[11] = total(prof[109:114])   ; K-index 0.26 A filter
  band[12] = total(prof,/double)    ; total intensity

  band[13] = total(prof[97:102])              ; ~ K1v
  band[14] = total(prof[127:136])   ; K-line red wing
  ;----------------------------------------
  ;-- center-of-gravity for the Ca profile 
  ;-- in 1.0 & 0.5 nm bands 
  ;----------------------------------------
  cagravity, 101, 121, prof, gravity1, k_amp1 ; 1 A   K-line
  cagravity, 106, 116, prof, gravity2, k_amp5 ; 0.5 A K-line
  band[15] = gravity1
  band[16] = k_amp1
  band[17] = gravity2
  band[18] = k_amp5
endif
  
                 
if (disper lt 3.)and(disper gt 1.5) then begin  ;-------------------------- 2.544 pm/px

  band[0] = total(prof[130:139])   ; ~ K1
  band[1] = total(prof[10:19])     ; W1  --- K-line blue wing
  band[2] = total(prof[30:39])     ; W2  --- K-line blue wing
  band[3] = total(prof[60:69])     ; W3  --- K-line blue wing
  band[4] = total(prof[200:209])   ; Mg II 2798  --- K-line red wing
  band[5] = total(prof[170:189])   ; W4  --- K-line red wing

  band[6] = total(prof[109:113])   ; K3  5 pixel 0.1 A 
  band[7] = total(prof[101:108])   ; K2v 8 pixel 0.2 A 
  band[8] = total(prof[114:121])   ; K2r 8 pixel 0.2 A 

  band[9] = total(prof[72:148])    ; K-index 1.0 A filter
  band[10] = total(prof[90:129])   ; K-index 0.506 A filter
  band[11] = total(prof[100:120])  ; K-index 0.26 A filter
  band[12] = total(prof,/double)   ; total intensity

  band[13] = total(prof[81:90])    ; ~ K1v
  band[14] = total(prof[134:152])  ; K-line red wing
  ;----------------------------------------
  ;-- center-of-gravity for the Ca profile 
  ;-- in 1.0 & 0.5 nm bands 
  ;----------------------------------------
  cagravity, 90, 130, prof, gravity1, k_amp1 ; 1.0 A K-line
  cagravity, 100, 120, prof, gravity2, k_amp5  ; 0.5 A K-line

  band[15] = gravity1
  band[16] = k_amp1
  band[17] = gravity2
  band[18] = k_amp5
endif


if (do_H_line eq 1)and(disper lt 3.) then begin  ;-------------- 2.544 pm/px
  band[19] = total(prof[390:394])    ; H3  0.1 A = 5 pixel
  band[20] = total(prof[382:389])    ; H2v 0.2 A = 8 pixel
  band[21] = total(prof[395:403])    ; H2r 0.2 A = 8 pixel

  band[22] = total(prof[360:420])    ; H-index 1.0 A filter
  band[23] = total(prof[380:400])    ; H-index 0.506 A filter
  band[24] = total(prof[385:395])    ; H-index 0.26 A filter
  band[25] = total(prof[360:370])    ; ~ H1

  cagravity, 372, 412, prof, gravity3, h_amp1  ; 1.0 A K-line
  cagravity, 382, 402, prof, gravity4, h_amp5  ; 0.5 A K-line
  band[26] = gravity3
  band[27] = h_amp1
  band[28] = gravity4
  band[29] = h_amp5
endif


if (do_H_line eq 1)and(disper gt 5.) then begin
  band[19] = total(prof[249:251])    ; H3  0.1 A = 3 pixel
  band[20] = total(prof[245:249])    ; H2v 0.2 A = 4 pixel
  band[21] = total(prof[252:256])    ; H2r 0.2 A = 4 pixel

  band[22] = total(prof[235:265])    ; H-index 1.0 A filter
  band[23] = total(prof[245:255])    ; H-index 0.506 A filter
  band[24] = total(prof[248:253])    ; H-index 0.26 A filter
  band[25] = total(prof[235:240])    ; ~ H1

  ;km = round((352. -110.) * 2.544 / disper + 110.) & kn = round((428. -110.) * 2.544 / disper + 110.)
  cagravity, 240, 260, prof, gravity3, h_amp1  ; 1.0 A K-line
  cagravity, 245, 255, prof, gravity4, h_amp5  ; 0.5 A K-line

  band[26] = gravity3
  band[27] = h_amp1
  band[28] = gravity4
  band[29] = h_amp5
endif
;###############################################################################################



if (do_gauss eq 1) then begin
;###############################################################################################
; perform  a single/double/triple Gaussian fit to the line core
; so far works only for disp = 2.5 and 5 mA
;###############################################################################################
  if (np gt  40) then  new_offset = line_offset + 6. else new_offset = line_offset + 3.
  kont = konti < abs(min(py))
  if (abs(abs(min(py)) - konti)/konti gt 0.5) then kont = konti
  kont = kont < 200. > 0.0
  np = n_elements(py)
  px = findgen(np)
  bbc = where(py gt (-10.), bppc)

  ee = ir_error(py, /nuv, /dark) ;0.29 * sqrt(abs(py))     ; see ir_error.pro for a description
  err_ave = ee ;< (abs(py)*0.9) ;make sure that the error value is larger than 0 and smaller than the data value

  cagravity, 0, np-1, py, gravity, h
  if (gravity lt 9.) then begin ; the case where Mg window is so small that having 110 pixels before the Mg II K core results
                                ; in mixing 280 nm data with Mg II data
     s = reform(py[0:9])
     sn = good_mean(py[10:*])
     chk_prf = where(s gt sn, count_s)
     if (count_s ge 1) then begin
        s(chk_prf) = sn
        py[0:9] = s
        cagravity, 0, np-1, py, gravity, h
     endif else begin
        gravity = np/2
     endelse   
  endif
  if (gravity gt np*0.75)or(~finite(gravity)) then gravity = np/2.
  posmax = gravity
  ;--------------------------------------------
  ;-- initial guess
  ;--------------------------------------------
  fit0 = [kont, max(py) * 0.8,  posmax , 24./dfac + randomn(seed)*0.1] * 1.0d

  if (do_H_line eq 0) then fit0[2] = gravity1 - new_offset else fit0[2] = gravity3 - new_offset
  fit0[2] = fit0[2] > (np/3.) < (np/1.3)
  kont = py[0] < 20. > 0.0
  fit0[0] = kont   ; the base level 
  
  ppy = py 
  range0=[0.5, 0.5, 0.2, 0.5] 
  dlambda = 1.0d

  ergs = my_sgf(px[bbc], ppy[bbc], ee[bbc], fit0, range0, dlambda[0], bbc, /double)
  fitsg = rgauss(px, [ergs.b, ergs.i1, ergs.p1,  ergs.w1])
  chisq1 = (1.0d/(n_elements(bbc) - 4.0d)) * total(((ppy[bbc] - fitsg[bbc])/err_ave[bbc])^2)
  is_bad = evaluate_sgf(fit0, ergs)
  if (abs(ergs.p1 - posmax) gt 15.) then is_bad = 1.  ; we know that Mg II cannot have large Doppler shift
  if (ergs.b lt 1.) then is_bad = 1.
  
  ;--------------------------------------------------------
  ;------------------- the second loop  ------------------- 
  ;--------------------------------------------------------
  if (ergs.w1/disper gt 8.)or(is_bad eq 1.)or(n_elements(fitsg) lt np)or(chisq1 gt 90.0) then begin
     py0 = median(py, 3)
     py = gauss_smooth(py0, 1.1, /edge_truncate)
     fit0[3] = (fit0[3] + randomn(seed)) > (fit0[3] - 2.0)
     fit0[1] = max(py) > 2.0 < 1.0d4
     range0=[0.2, 0.6, 0.6, 0.6] 
     ergs_new = my_sgf(px[bbc], ppy[bbc], ee[bbc], fit0, range0, dlambda[0], bbc, /double)
     fitsg_new = rgauss(px, [ergs_new.b, ergs_new.i1, ergs_new.p1,  ergs_new.w1])
     chisq1_new = (1.0d/(n_elements(bbc) - 4.0d)) * total(((ppy[bbc] - fitsg[bbc])/err_ave[bbc])^2)
     if (chisq1_new lt chisq1) then begin
        chisq1 = chisq1_new
        fitsg = fitsg_new
        ergs = ergs_new
     endif   
  endif
  
  spar1=[ergs.b, ergs.p1, ergs.i1, ergs.w1, chisq1, reform(ergs.sigma)]
  line_core  = ergs.p1

  ;------------------------------------------------------------------------------
  ;  if the fit is not satisfactory, then use random initialization 
  ;------------------------------------------------------------------------------
  if (chisq1 gt 1.)and(do_improve eq 1) then begin
     if (chisq1 le 800.) then begin
        random_sg_fit, (50*zeta), px, py, ppy, ee, err_ave, bbc, 0.6, ergs, spar_new, fitsg_new
        if (spar_new[4] lt chisq1) then begin
           chisq1 = spar_new[4]
           fitsg = fitsg_new
           spar1 = spar_new
        endif
     endif
     if (chisq1 gt 800.) then begin
        random_sg_fit, (50*zeta), px, py, ppy, ee, err_ave, bbc, 0.9, ergs, spar_new, fitsg_new
        if (spar_new[4] lt chisq1) then begin
              chisq1 = spar_new[4]
              fitsg = fitsg_new
              spar1 = spar_new
        endif
     endif
  endif   

  
;---------------------------------------------------
;--  a double Guassian fit to the line profile
;---------------------------------------------------
  ergs.w1 = ergs.w1 > 1. < 12.      ; in any case
  v_peak = max(py[0:(posmax-1)], v_posmax)  &  v_posmax = v_posmax > (ergs.p1 - disper*3.) 
  r_peak = max(py[(posmax+1):*], r_posmax)  &  r_posmax = (r_posmax + (posmax+1)) > (ergs.p1) < (np -5.) 
  width = (ergs.w1 * 0.25) > 3.0
  range0 = [0.8, 0.4, 0.8] 
  range1 = [0.8, 0.4, 0.8] 
  fit0 = [v_peak, v_posmax, width] * 1.0d
  fit1 = [r_peak, r_posmax, width] * 1.0d
   
  dlambda = 1.50d
  ergd = my_dgf_mg(px[bbc],ppy[bbc], ee[bbc], fit0, fit1, range0, range1, dlambda[0], bbc, /double)
  fitdg = ergd.i1*exp(-((px - ergd.p1)/ergd.w1)^2) + ergd.i2*exp(-((px - ergd.p2)/ergd.w2)^2)
  chisq2 = (1.0d/(n_elements(bbc) - 6.0d)) * total(((ppy[bbc] - fitdg[bbc])/err_ave[bbc])^2)
  dpar1=[ergd.p1, ergd.i1, ergd.w1, ergd.p2, ergd.i2, ergd.w2, chisq2, reform(ergd.sigma)]

  ;--------------------------------------------------------
  ;------------------- the second loop  ------------------- 
  ;--------------------------------------------------------
  is_bad = evaluate_dgf_mg(fit0, ergd)
  if ((ergd.w1 < ergd.w2) lt 3.)and(disper lt 3.) then is_bad = 1
  if (ergd.w1/disper gt 20.)or(is_bad eq 1.)or(n_elements(fitdg) lt np)or(chisq2 gt 30.) then begin
     ;print, '----------'
     py0 = median(py,3)
     py = gauss_smooth(py0, 1.1, /edge_truncate)

     if (disper lt 3.) then width = width > 4.
     fit0 = [v_peak, v_posmax, width+0.1] * 1.0d
     fit1 = [r_peak, r_posmax, width] * 1.0d

     range0 = [0.6, 0.3, 0.6]  * 1.0d
     range1 = [0.6, 0.2, 0.6]  * 1.0d
     ergd_new = my_dgf_mg(px[bbc],py[bbc], ee[bbc], fit0, fit1, range0, range1, dlambda[0], bbc, /double)
     fitdg_new = ergd_new.i1*exp(-((px - ergd_new.p1)/ergd_new.w1)^2) + ergd_new.i2*exp(-((px - ergd_new.p2)/ergd_new.w2)^2)
     chisq2_new = (1.0d/(n_elements(bbc) - 6.0d)) * total(((py[bbc] - fitdg_new[bbc])/err_ave[bbc])^2)
     if (chisq2_new lt chisq2) then begin
        chisq2 = chisq2_new
        fitdg = fitdg_new
        ergd = ergd_new
     endif   
  dpar1=[ergd.p1, ergd.i1, ergd.w1, ergd.p2, ergd.i2, ergd.w2, chisq2, reform(ergd.sigma)]
  endif

  ;--------------------------------------------------------
  ;------------------- the third loop  -------------------- 
  ;--------------------------------------------------------
  if (chisq1 lt chisq2)or(chisq2 gt 20.) then begin  ;----- if single Gaussian was a better fit
     ;print, '+++++++++++++++++++++++++++++++++'
     fit0 = [ergs.i1, ergs.p1, ergs.w1]
     fit1 = [ergd.i1 * 0.5, ergs.p1+randomn(seed), ergs.w1+randomn(seed)*0.1]
     range0 = [0.3, 0.05, 0.5]
     range1 = [0.9, 0.2, 0.5]
     ergd_new = my_dgf_mg(px[bbc],ppy[bbc], ee[bbc], fit0, fit1, range0, range1, dlambda[0], bbc, /double)
     fitdg_new = ergd_new.i1*exp(-((px - ergd_new.p1)/ergd_new.w1)^2) + ergd_new.i2*exp(-((px - ergd_new.p2)/ergd_new.w2)^2)
     chisq2_new = (1.0d/(n_elements(bbc) - 6.0d)) * total(((ppy[bbc] - fitdg_new[bbc])/err_ave[bbc])^2)
     if (chisq2_new lt chisq2) then begin
        chisq2 = chisq2_new
        fitdg = fitdg_new
        ergd = ergd_new
     endif   
  dpar1=[ergd.p1, ergd.i1, ergd.w1, ergd.p2, ergd.i2, ergd.w2, chisq2, reform(ergd.sigma)]
  endif

  ;--------------------------------------------------------
  ;------------------- the fourth loop  ------------------- 
  ;--------------------------------------------------------
  if (chisq2 gt 100.) then begin ; --- in case of wide emission lobes
     fit0 = [v_peak, v_posmax, width*1.1] * 1.0d
     fit1 = [r_peak, r_posmax, width*1.1] * 1.0d

     range0 = [0.99, 0.9, 0.99]  * 1.0d
     range1 = [0.99, 0.7, 0.99]  * 1.0d

     ergd_new = my_dgf_mg(px[bbc],ppy[bbc], ee[bbc], fit0, fit1, range0, range1, dlambda[0], bbc, /double)
     fitdg_new = ergd_new.i1*exp(-((px - ergd_new.p1)/ergd_new.w1)^2) + ergd_new.i2*exp(-((px - ergd_new.p2)/ergd_new.w2)^2)
     chisq2_new = (1.0d/(n_elements(bbc) - 6.0d)) * total(((ppy[bbc] - fitdg_new[bbc])/err_ave[bbc])^2)
     if (chisq2_new lt chisq2) then begin
        chisq2 = chisq2_new
        fitdg = fitdg_new
        ergd = ergd_new
     endif   
  dpar1=[ergd.p1, ergd.i1, ergd.w1, ergd.p2, ergd.i2, ergd.w2, chisq2, reform(ergd.sigma)]
  endif
  
  ;------------------------------------------------------------------------------
  ;  if the fit is not satisfactory, then use random initialization 
  ;------------------------------------------------------------------------------
  if (chisq2 gt 1.)and(do_improve eq 1) then begin
     if (chisq2 le 10.) then begin
        new_mg2_fit, (100*zeta), px, py, ppy, ee, err_ave, bbc, 0.3, ergd, dpar_new, fitdg_new
        if (dpar_new[6] lt chisq2) then begin
           chisq2 = dpar_new[6]
           fitdg = fitdg_new
           dpar1 = dpar_new
        endif
     endif
     if (chisq2 gt 10.)and (chisq2 le 40.) then begin
        new_mg2_fit, (100*zeta), px, py, ppy, ee, err_ave, bbc, 0.4, ergd, dpar_new, fitdg_new
        if (dpar_new[6] lt chisq2) then begin
           chisq2 = dpar_new[6]
           fitdg = fitdg_new
           dpar1 = dpar_new
        endif
     endif
     if (chisq2 gt 40.) then begin
        new_mg2_fit, (150*zeta), px, py, ppy, ee, err_ave, bbc, 0.9, ergd, dpar_new, fitdg_new
        if (dpar_new[6] lt chisq2) then begin
               chisq2 = dpar_new[6]
              fitdg = fitdg_new
              dpar1 = dpar_new
        endif
     endif
  endif


;---------------------------------------------------
;--  a triple Guassian fit to Mg II h/k
;---------------------------------------------------
; broad component core position, amplitude, and width
  fit0 = [max(py) * 0.5,  posmax, 7.8 *disper/dfac] 
  ; initial guess for amplitudes, line centers, and widths of V and R peaks
  fit1 = [max(py) * 0.9, gravity - 6./dfac, 3.0 * disper/dfac, max(py) * 0.8, gravity + 6./dfac, 3.0 * disper/dfac]
 
  range0 = [0.8, 0.3, 0.8] ; 0.5, 0.2, 0.8
  range1 = [0.8, 0.3, 0.4, 0.8, 0.2, 0.4]
  fit0_0 = fit0
  fit1_0 = fit1
  
  if (chisq2 lt 50.) then begin  ;----- if double Gaussian was a good fit
    fit0[1] = (ergd.p1 + ergd.p2) * 0.5
    fit1 = [ergd.i1, ergd.p1, ergd.w1 < fit0[2]*0.7, ergd.i2, ergd.p2, ergd.w1 < fit0[2]*0.7]
  endif   
 
  dlambda = 1.0d
  pyn = ppy 
  ergt = my_tgf_mg(px[bbc], pyn[bbc], ee[bbc], fit0, fit1, range0, range1, dlambda[0], bbc, /double)
  fittg = tgf_mg(px, [ergt.i1, ergt.p1, ergt.w1, ergt.i2, ergt.p2, ergt.w2, ergt.i3, ergt.p3, ergt.w3])
  chisq4 = (1.0d/(n_elements(bbc) - 9.0d)) * total(((pyn[bbc] - fittg[bbc])/err_ave[bbc])^2)
  tpar1=[ergt.p1, ergt.i1, ergt.w1, ergt.p2, ergt.i2, ergt.w2, ergt.p3, ergt.i3, ergt.w3, chisq4, reform(ergt.sigma)]
  w = ergt.w1 < ergt.w2




  ;--------------------------------------------------------
  ;------------------- the second loop  ------------------- 
  ;--------------------------------------------------------
  bad = 0.
  if ((ergt.i2/ergt.i3) gt 10.)or((ergt.i3/ergt.i2) gt 10.)or(w le 2./dfac) then bad = 1.
  if (ergt.w1 lt 22./disper)or(n_elements(fittg) lt np) then bad = 1.
  if (chisq4 gt chisq2)or(chisq4 gt 4.) then bad = 1.
  if (bad ne 0.) then begin
     ;---------------------------------------------
     syn_k3_aux = reform(fittg)
     dg_x1 = ergd.p1 > 2.
     dg_x2 = ergd.p2 < (np-2)
     if (dg_x1 gt dg_x2) then begin  ; sort the two Gaussians if required
        dg_x1 = ergd.p2 > 2. < (np/2)
        dg_x2 = ergd.p1 < (np-2) 
     endif
     if (dg_x1 ge dg_x2) then stop
     syn_k3_aux = reform(syn_k3_aux[dg_x1:dg_x2])
     syn_aux = min(syn_k3_aux, pos_aux)
     if (pos_aux eq 0) then pos_aux = n_elements(syn_k3_aux)/2.
     k3w_init = pos_aux
     k3_aux = syn_k3_aux[pos_aux] * 0.7
     ;---------------------------------------------
     py0 = median(py,3)
     py = gauss_smooth(py0, 1.1, /edge_truncate)
     fit0[0] = k3w_init ;(ergd.p1 + ergd.p2) * 0.5d
     fit0[1] = k3_aux
     fit1[0] = max(py) * 0.75
     fit1[3] = max(py) * 0.85
     range0 = [0.6, 0.3, 0.75]
     range1 = [0.95, 0.4, 0.75, 0.95, .3, 0.75]
     
     if (chisq2 lt chisq4) then begin  ;----- if double Gaussian was a better fit
       fit1 = [ergd.i1, ergd.p1, ergd.w1*0.9, ergd.i2, ergd.p2, ergd.w1*0.9]
       range1 = fltarr(6) + 0.6
       range0[1] = 0.6
     endif  

     
     if ((chisq1 lt chisq4)and(chisq1 lt chisq2))or(chisq1 lt 10.) then begin  ;----- if single Gaussian was the best fit
       fit0 = [ergs.i1, ergs.p1, ergs.w1]
       fit1 = [ergd.i1*0.2, ergd.p1+randomn(seed), ergd.w1*0.5, ergs.i1*0.2, ergs.p1+randomn(seed), ergd.w1*0.5]
       range0 = [0.3, 0.2, 0.25]
       range1 = [0.5, 0.1, 0.5, 0.5, 0.1, 0.5]
    endif
     
     dlambda = 1.0d
     ergt_new = my_tgf_mg(px[bbc], pyn[bbc], ee[bbc], fit0, fit1, range0, range1, dlambda[0], bbc, /double)
     fittg_new = tgf_mg(px, [ergt_new.i1, ergt_new.p1, ergt_new.w1, ergt_new.i2, ergt_new.p2, ergt_new.w2, ergt_new.i3, ergt_new.p3, ergt_new.w3])
     chisq4_new = (1.0d/(n_elements(bbc) - 9.0d)) * total(((pyn[bbc] - fittg_new[bbc])/err_ave[bbc])^2)
     if (chisq4_new lt chisq4) then begin
        chisq4 = chisq4_new
        fittg = fittg_new
        ergt = ergt_new
     endif   
     tpar1=[ergt.p1, ergt.i1, ergt.w1, ergt.p2, ergt.i2, ergt.w2, ergt.p3, ergt.i3, ergt.w3, chisq4, reform(ergt.sigma)]
  endif
  
  if 0 then begin

  ;--------------------------------------------------------
  ;------------------- the third loop  -------------------- 
  ;--------------------------------------------------------
  if (chisq4 gt 10.) then begin  ;---- in case of pathological profiles 
     fit0 = fit0_0  &  fit0[2] *= 2.0
     fit1 = fit1_0
     range0=[1.2, 0.8, 0.6] 
     range1=[0.9, 0.8, 0.6, 0.9, 0.8, 0.6]
     ergt_new = my_tgf_mg(px[bbc], pyn[bbc], ee[bbc], fit0, fit1, range0, range1, dlambda[0], bbc, /double)
     fittg_new = tgf_mg(px, [ergt_new.i1, ergt_new.p1, ergt_new.w1, ergt_new.i2, ergt_new.p2, ergt_new.w2, ergt_new.i3, ergt_new.p3, ergt_new.w3])
     chisq4_new = (1.0d/(n_elements(bbc) - 9.0d)) * total(((pyn[bbc] - fittg_new[bbc])/err_ave[bbc])^2)

     if (chisq4_new lt chisq4) then begin
        chisq4 = chisq4_new
        fittg = fittg_new
        ergt = ergt_new
     endif   
     tpar1=[ergt.p1, ergt.i1, ergt.w1, ergt.p2, ergt.i2, ergt.w2, ergt.p3, ergt.i3, ergt.w3, chisq4, reform(ergt.sigma)]  
  endif
  ;--------------------------------------------------------
  ;------------------- the fourth loop  -------------------- 
  ;--------------------------------------------------------
  if (chisq4 gt 10.) then begin  ;---- in case of pathological profiles 
    fit0 = [max(py) * 0.5,  (ergd.p1 + ergd.p2)*0.5, 7.8 *disper/dfac] 
    fit1 = [max(py) * 0.9, gravity - 6./dfac, 3.0 * disper/dfac, max(py) * 0.8, gravity + 6./dfac, 3.0 * disper/dfac]
    range0 = [0.8, 0.3, 0.8] ; 0.5, 0.2, 0.8
    range1 = [0.8, 0.3, 0.4, 0.8, 0.2, 0.4]

     fit0 = fit0_0  &  fit0[2] = (10.+randomn(seed))/dfac
     fit1 = fit1_0
     fit1[0] *= 0.75 & fit1[3] *= 0.75 &  fit0[0] *= 0.75
     fit0[1] = posmax
     range0 = [1.2, 0.2, 0.6] 
     range1 = [0.9, 0.3, 0.6, 0.9, 0.2, 0.6]
     ergt_new = my_tgf_mg(px[bbc], pyn[bbc], ee[bbc], fit0, fit1, range0, range1, dlambda[0], bbc, /double)
     fittg_new = tgf_mg(px, [ergt_new.i1, ergt_new.p1, ergt_new.w1, ergt_new.i2, ergt_new.p2, ergt_new.w2, ergt_new.i3, ergt_new.p3, ergt_new.w3])
     chisq4_new = (1.0d/(n_elements(bbc) - 9.0d)) * total(((pyn[bbc] - fittg_new[bbc])/err_ave[bbc])^2)
     if (chisq4_new lt chisq4) then begin
        chisq4 = chisq4_new
        fittg = fittg_new
        ergt = ergt_new
     endif   
     tpar1=[ergt.p1, ergt.i1, ergt.w1, ergt.p2, ergt.i2, ergt.w2, ergt.p3, ergt.i3, ergt.w3, chisq4, reform(ergt.sigma)]
  endif
  ;--------------------------------------------------------
  ;------------------- the fifth loop  -------------------- 
  ;--------------------------------------------------------
  if (chisq4 gt 10.)and(max(py) gt 600.) then begin  ;---- in case of network/active profiles 
    fit0 = [max(py) * 045,  posmax, 7.8 *disper/dfac] ; core component
    fit1 = [max(py) * 0.9, gravity - 7./dfac, 3.1 * disper/dfac, max(py) * 0.8, gravity + 7./dfac, 3.1 * disper/dfac] ; H2v H2v
    range0 = [0.6, 0.3, 0.8] 
    range1 = [0.6, 0.3, 0.5, 0.6, 0.2, 0.5]

    fit0 = fit0_0  &  fit0[2] = (10.+randomn(seed))/dfac
    fit1 = fit1_0
    fit1[0] *= 0.75 & fit1[3] *= 0.75 &  fit0[0] *= 0.75
    fit0[1] = posmax
    range0 = [1.2, 0.2, 0.6] 
    range1 = [0.9, 0.3, 0.6, 0.9, 0.2, 0.6]
    ergt_new = my_tgf_mg(px[bbc], pyn[bbc], ee[bbc], fit0, fit1, range0, range1, dlambda[0], bbc, /double)
    fittg_new = tgf_mg(px, [ergt_new.i1, ergt_new.p1, ergt_new.w1, ergt_new.i2, ergt_new.p2, ergt_new.w2, ergt_new.i3, ergt_new.p3, ergt_new.w3])
    chisq4_new = (1.0d/(n_elements(bbc) - 9.0d)) * total(((pyn[bbc] - fittg_new[bbc])/err_ave[bbc])^2)
    if (chisq4_new lt chisq4) then begin
        chisq4 = chisq4_new
        fittg = fittg_new
        ergt = ergt_new
    endif   
    tpar1=[ergt.p1, ergt.i1, ergt.w1, ergt.p2, ergt.i2, ergt.w2, ergt.p3, ergt.i3, ergt.w3, chisq4, reform(ergt.sigma)]  
  endif
  
  ;--------------------------------------------------------
  ;------------------- the sixth loop  ------------------- 
  ;--------------------------------------------------------
  if (ergt.w1 lt ergt.w2)and(ergt.w1 lt 6./dfac)and(ergt.w2 gt 7./dfac) then begin ;---- in case of too narrow k3 peak
    fit0 = [ergt.i1, ergt.p1, ergt.w2] 
    fit1 = [ergt.i2, ergt.p2, ergt.w1, ergt.i3, ergt.p3, ergt.w3]
    range0 = [0.2, 0.2, 0.3] 
    range1 = [0.2, 0.2, 0.3, 0.2, .2, 0.3]
    ergt_new = my_tgf_mg(px[bbc], pyn[bbc], ee[bbc], fit0, fit1, range0, range1, dlambda[0], bbc, /double)
    fittg_new = tgf_mg(px, [ergt_new.i1, ergt_new.p1, ergt_new.w1, ergt_new.i2, ergt_new.p2, ergt_new.w2, ergt_new.i3, ergt_new.p3, ergt_new.w3])
    chisq4_new = (1.0d/(n_elements(bbc) - 9.0d)) * total(((pyn[bbc] - fittg_new[bbc])/err_ave[bbc])^2)

    if (chisq4_new lt chisq4) then begin
        chisq4 = chisq4_new
        fittg = fittg_new
        ergt = ergt_new
    endif   
    tpar1=[ergt.p1, ergt.i1, ergt.w1, ergt.p2, ergt.i2, ergt.w2, ergt.p3, ergt.i3, ergt.w3, chisq4, reform(ergt.sigma)]  
  endif
endif
  
  
  ;------------------------------------------------------------------------------
  ;  if the fit is not satisfactory, then use random initialization 
  ;------------------------------------------------------------------------------
  if (chisq4 gt 1.0)and(do_improve eq 1) then begin
     if (chisq4 le 5.) then begin
        new_mg3_fit, (100*zeta), px, py, ppy, ee, err_ave, bbc, 0.3, ergt, tpar_new, fittg_new
        if (tpar_new[9] lt chisq4) then begin
           chisq4 = tpar_new[9]
           fittg = fittg_new
           tpar1 = tpar_new
        endif
     endif
     if (chisq4 gt 5.) then begin
        new_mg3_fit, (100*zeta), px, py, ppy, ee, err_ave, bbc, 0.4, ergt, tpar_new, fittg_new
        if (tpar_new[9] lt chisq4) then begin
           chisq4 = tpar_new[9]
           fittg = fittg_new
           tpar1 = tpar_new
        endif
     endif
     if (chisq4 gt 10.) then begin
        new_mg3_fit, (250*zeta), px, py, ppy, ee, err_ave, bbc, 0.9, ergt, tpar_new, fittg_new
        if (tpar_new[9] lt chisq4) then begin
              chisq4 = tpar_new[9]
              fittg = fittg_new
              tpar1 = tpar_new
        endif
     endif
  endif
  
endif

;################################### end of Gaussian fit #######################################
;###############################################################################################
if (plt ge 1)and(do_gauss eq 1) then begin
  print, '------------------------------------------------------'
  print, 'single Gaussain: ', spar1[0:4]
  print, 'double Gaussain: ', dpar1[0:6]
  print, 'triple Gaussain: ', tpar1[0:9]
  print, '------------------------------------------------------'
  plot,  px, py, /xst, psym=-1
  loadct,40,/silent
  oplot, px, fitsg, color=75, thick=2
  oplot, px, fitdg, color=175, thick=2
  oplot, px, fittg, color=245, thick=2
  loadct,0,/silent
  ans =''  &  read, ans, prompt='Press enter to continue ...'
endif
;############################################################################################
;############################################################################################
;----------------------------
;--  Set output parameters
;----------------------------

erg.type = type           ; 1 parameter
erg.hpar = hpar           ; 6
erg.eml = emb             ; 12
erg.ems = embv            ; 4
erg.spar = spar1          ; 9
erg.dpar = dpar1          ; 15
erg.tpar = tpar1          ; 21

erg.sprf = py
erg.sfit = fitsg
erg.dfit = fitdg
erg.tfit = fittg

erg.sub_wing = subtracted_wing


erg.velpos = velpos       ; 6
erg.fe_int = fe_int       ; 5
erg.band   = band         ; 11
erg.xmins  = xmins        ; 4X2
erg.xmaxs  = xmaxs        ; 3X2
erg.xmh1  = minh1         ; 2X2


return,erg
end
