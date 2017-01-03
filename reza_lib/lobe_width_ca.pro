pro lobe_width_ca, pr, prof, x_lobe, y_lobe, ref, fwhm, xl, xr, mg=mg

;+
;===============================================================
; procedure: lobe_width_ca.pro
;
; purpose: find FWHM of an emission peak in Ca H profile
;
; prof: profile
; x_lobe/y_lobe: position/amplitude of a single emission peak
; ref: H3 minimum intensity, or generally the continuum level
; fwhm: FWHM in pixel unit
;
; Jan 14, 2007 : created  
; May 08, 2014 : a new switch to be used in case of Mg II line profiles
; May 13, 2014: bug fix for the first last point of a lobe
; May 14, 2014: added 1/e width
;
; R.Rezaei @ KIS                         e-mail:  rrezaei@iac.es      
;===============================================================
;-       


if (n_elements(mg) le 0) then mg = 0 else mg = 1

; 1/e = 0.367879
efac =  0.367879
prof = smooth(prof, 3)

if (x_lobe lt (n_elements(prof)/4.))and(mg eq 0) then begin
   r=max(prof, pos)
   x_lobe = pos
   y_lobe = prof[pos]
endif else begin
   r=max(prof[(x_lobe - 5):(x_lobe + 5)], pos)
   x_lobe = round(pos + x_lobe - 5)
   y_lobe = prof[pos]
   if (prof[x_lobe-1] gt prof[x_lobe]) then begin
     x_lobe -= 1
     y_lobe = prof[x_lobe] 
   endif
endelse


  peak = round(x_lobe)
  nt = n_elements(prof)
;pr = 0
  ;w = max(prof, ss)
  left = peak  &  right = peak    
  ;print, left, right, prof[peak]

  xrf= (y_lobe-ref)*0.5 + ref
  i = peak
  while (prof[i] gt xrf)and(i gt 2)and(prof[i-1] lt prof[i]) do  begin 
     ;print, i, prof[i]
     i -= 1   
  endwhile
  left = i
  ;print, left
  if ((peak-left) le 1) then left -= 1
  if ((peak-left) gt 15.) then left = peak - 2.
  ;----------------------------------------------
  ;print, 'aaa', left, right
  i = peak
  while (prof[i] gt xrf)and(i lt nt-2)and(prof[i+1] lt prof[i]) do begin 
     ;print, i, prof[i]
     i += 1
  endwhile
  right = i-1
  ;print, right
  if ((right-peak) eq 1) then right += 1
  if ((right-peak) gt 15.) then right = peak + 2.
  ;--------------------------
  x1 = left & x2 = x1 + 1
  xl = (xrf - prof[x1])/(prof[x2]-prof[x1]) + float(x1) 
  if (abs(xl - left) gt 1.) then xl = left 
  ;------------------------------------------------
  x2 = right & x1 = x2 -1.
  xr = (xrf - prof[x1])/(prof[x2]-prof[x1]) + float(x1) 
  if (abs(xr - right) gt 1.) then xr = right 
  fwhm = xr - xl
  ;print, left, right, xl, xr, fwhm
  ;print, '-----------------------------------------'
  ;----------------------------------------------



if (fwhm ge 50.) then fwhm = 0.01
if (fwhm lt 0.) then fwhm = 100.

if pr then print, xl, xr, left, right

if pr and (mg eq 1) then begin ; Mg
!p.multi=0
pst = 80
  plot, findgen(n_elements(prof[pst:*])), prof[pst:*],/yno
  oplot, findgen(n_elements(prof[pst:*])), prof[pst:*],psym=1,symsize=.7
  hline, xrf
  vline, xr-pst
  vline, xl-pst
  ;plots, right-pst, prof[right], psym=3
  ;plots, left-pst, prof[left], psym=3
  oplot, [right,right] - pst, [0,10], linestyle=1
  oplot, [left,left] - pst, [0,10], linestyle=1

wait, .5
;paused
endif
if pr and (mg eq 0) then begin  ; Ca
!p.multi=0
pst = 210
  plot, findgen(n_elements(prof[pst:310])), prof[pst:310],/yno
  oplot, findgen(n_elements(prof[pst:310])), prof[pst:310],psym=1,symsize=.7
  hline, xrf
  vline, xr-pst
  vline, xl-pst
  plots, right-pst, prof[right], psym=3
  plots, left-pst, prof[left], psym=3
  oplot, [right,right] - pst, [0,10], linestyle=1
  oplot, [left,left] - pst, [0,10], linestyle=1

wait, .5
;paused
endif

if 0 then begin
  xrf= (y_lobe-ref)*efac + ref
  i = peak
  while (prof[i] gt xrf)and(i gt 2) do    i -= 1   
  left = i
  if ((peak-left) le 1) then left -= 1
  if ((peak-left) gt 20.) then left = peak - 2.
  ;----------------------------------------------
  i = peak
  while (prof[i] gt xrf)and(i lt nt-2) do i += 1
  right = i
  if ((right-peak) eq 1) then right += 1
  if ((right-peak) gt 20.) then right = peak + 2.
  ;------------------------
  x1 = left & x2 = x1 + 1
  xl = (xrf - prof[x1])/(prof[x2]-prof[x1]) + float(x1) 
  ;----
  x2 = right & x1 = x2 -1.
  xr = (xrf - prof[x1])/(prof[x2]-prof[x1]) + float(x1) 
  ew = xr - xl
  ;----------------------------------------------
endif


end
