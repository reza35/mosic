pro iris_timeseri, filepath

;+
;===============================================================
; procedure : iris_timeseri.pro
;  
; purpose : to flatten timeseries to have lambda as first axis
;
; Sep 25, 2014 : created
;
; R.Rezaei @ IAC                         e-mail:  rrezaei@iac.es      
;===============================================================
;-
  
inpath = filepath
files_d = read_dir(inpath, filter='iris*_im.fits')

num_files = n_elements(files_d.files)
data = inpath + files_d.files[0]

a = readfits(data, hed)

vtf = size(a)
print, vtf

if (vtf[0] eq 4) then begin
  cx = vtf[1]*vtf[4]
  b = fltarr(cx, vtf[2], vtf[3])
  for i=0,vtf[4]-1 do b[(vtf[1]*i):(vtf[1]*(i+1)-1),*,*] = a[*,*,*,i]
  a = b
  b = 0.
endif else begin
  cx = vtf[2]
  b = fltarr(cx, vtf[1], vtf[3])
  for i=0,vtf[3]-1 do b[*,*, i] = transpose(a[*,*,i])
  a = b
  b = 0.
endelse


writefits, data, a, hed

end
  
