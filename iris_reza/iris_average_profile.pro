pro iris_average_profile, inpath, nsl=nsl

;+
;===============================================================
; function : iris_average_profile.pro
;  
; purpose : to create average in each wavelength band
;
; Note: the printed exposure time in case of flattend timeseries is NOT correct.
;  
; Sep 25, 2014 : created
; May 08, 2016 : updated for multi file structure  
; Jun 28, 2016 : read the data cube slice-by-slice, so doesn't
;                 matter how big it is
;
; R.Rezaei @ IAC                         e-mail:  rrezaei@iac.es      
;===============================================================
;-
on_error, 2
  
if (n_elements(nsl) eq 0)  then nsl=0  else nsl=1

check = read_dir(inpath , filter='aver_prof_*.sav')
num = n_elements(check.files)
if ( num ge 1)and(strlen(check.files[0]) gt 5) then return
  
files_d = read_dir(inpath , filter='iris*_im.fits')
num = n_elements(files_d.files)

for kapa=0, num-1 do begin
  data = inpath + files_d.files[kapa]

  ch = strsplit(files_d.files[kapa], '., _', /extract)
  term = ch[2]+'_'+ch[3]+'_'+ch[6]


  a=readfits(data, nslice=0) * 1.0d
  q1 = where(~finite(a), count) &  if (count gt 0) then a(q1) = -1.0  ; for individual 2D rasters
  hed = headfits(data)
  n_slit = fxpar(hed, 'naxis3')
  print, n_slit
  for i=1,n_slit-1 do begin
     b=readfits(data, nslice=i, /silent)
     q1 = where(~finite(b), count) &  if (count gt 0) then b(q1) = -1.0  ; for individual 2D rasters
     a += b
     print, i
  endfor
  a /= float(n_slit)
  if (kapa eq 0) then begin
    exp = fxpar(hed, 'CDELT4')/float(fxpar(hed, 'NAXIS3 '))
    print, 'exposure   x0    y0   '
    x0 = fxpar(hed, 'XCEN')
    y0 = fxpar(hed, 'YCEN')
    print, exp, x0, y0,format='(F6.3,F8.1,F9.1,/)'
  endif

  tay = size(a)
 
  a = total(a, 2)/float(tay[2]) ; average spectrum along the slit

  avprof = a
  plot, a, /xstyle, psym=-1
  ;stop
  save, filename=inpath+'aver_prof_'+term+'.sav', avprof, hed
endfor


end
