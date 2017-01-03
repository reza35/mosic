pro reorder, inpath

;+
;===============================================================
; procedure :  reorder.pro
;
; purpose : to transpose IRIS Level 3 data to have lambda as first axis, slit as
; second, and the scan as the third. If the data has 4 dimensions, it should be splitted in advance.
;
; reorder, '/file/path/'
;
; Apr 15, 2016 : created  
; Jun 28, 2016 : read the data cube slice-by-slice, so doesn't
;                 matter how big it is
; Sep 10, 2016 : bug fix
; Oct 05, 2016 : improved documentation
;  
; R.Rezaei @ IAC                         e-mail:  rrezaei@iac.es      
;===============================================================
;-

files_d = read_dir(inpath, filter='iris_l3*_im.fits')

num_files = n_elements(files_d.files)
if (num_files lt 1) then begin
   print, 'The calibrated file does not exist.'
   stop
endif   
if (num_files gt 1) then begin
   print, 'There are multiple calibrated files in the directory.'
   stop
endif   

file = inpath+files_d.files[0]

print, file

hed = headfits(file)
n3 = fxpar(hed, 'naxis3')
n1 = fxpar(hed, 'naxis1')
n2 = fxpar(hed, 'naxis2')
print, n1, n2, n3
a = fltarr(n3, n2, n1)

for i=0,n3-1 do begin
     b = readfits(file, nslice=i, /silent)
     sb = where(~finite(b), count)  &  if (count gt 0) then b(sb) = -1.  ; for individual 2D rasters
     a[i,*,*] = transpose(b)
     print, n3-i-1, string(13b), format='(I,a1,$)'
endfor

writefits, file, a, hed

end

