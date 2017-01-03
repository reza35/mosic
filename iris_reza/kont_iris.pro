pro kont_iris, polis_path, only_mg = only_mg, normal=normal

;
;+
;===============================================================  
; procedure: kont_iris.pro
;  
; purpose : to create a simple qusi-continuum map from IRIS data,
;  so we get an idea what is in the data !
;  
; 2015/12/26 : created
; 2016/01/10 : works with splitted L3 data, still incomplete !
; 2016/05/25 : works only with reordered data structure [lambda,slit, scan] 
; 2016/06/01 : takes C II range to have some idea about off-disk area
; 2016/06/15 : read the data slice-by-slice
; 2016/11/14 : big fix for Jpeg continuum map
;
; R.Rezaei @ IAC                         e-mail:  rrezaei@iac.es      
;===============================================================
;-
  
common share_timeseri, ts


h = check_if_file_exists(polis_path+'kont.sav')
if (string2num(h[1]) eq 1) then return


if (n_elements(normal) ne 0) then normal=1 else normal=0
if (n_elements(only_mg) ne 0) then only_mg=1 else only_mg=0

if (only_mg eq 0)or(normal eq 1) then files_d = read_dir(polis_path, filter='iris_l3*_all_im.fits')
if (only_mg eq 1) then files_d = read_dir(polis_path, filter='iris*_MgIIk2796_im.fits')

infile= polis_path + files_d.files[0]

print, only_mg
print, infile
print, '-------------------------'

hed = headfits(infile)
n1 = fxpar(hed, 'NAXIS1')
n2 = fxpar(hed, 'NAXIS2')
n3 = fxpar(hed, 'NAXIS3')
kont_map = fltarr(n3, n2)
a = readfits(infile, nslice=0, /silent)

for i=1, n3-1 do begin
     b = readfits(infile, nslice=i, /silent)
     q1 = where(~finite(b), count) &  if (count gt 0) then b(q1) =  0.
     a += b
     kont_map[i, *] = total(median(b[round(n1*0.5):*, *], 3), 1) / float(n1/2.)
     print, i
endfor
a /= float(n3)

;stop
print, infile
q1 = where(~finite(a), count) &  if (count gt 0) then a(q1) =  0.
tay = size(a)
v = a ;total(a[w1:w2, *, *], 1) > 0.



stdev_despike, kont_map, 10, 4
s = percentiles(kont_map, value=[0.05, 0.5, 0.999])
v = kont_map > s[0] < s[2]
;v = transpose(v)
s = strsplit(infile, '.',/extract)
write_jpeg, s[0]+'.jpg', quality=100, bytscl(v)
print, s[0]+'.jpg'


kont = transpose(v)
u = total(kont, 2)/10000.
print, 'select the proper y-range for data analysis'
window, xs=800, ys=600
read_plot_click, 2, u, u, pos
n = n_elements(u)
if (pos[0] lt 0) then pos[0] = 0
if (pos[1] gt (n-1)) then pos[1] = n-1

if (pos[0] eq 0)and(pos[1] eq (n-1)) then begin
   print
   print, 'Full slit height was selected !'
   print
   wait, 1.
endif

cfile = polis_path + 'kont.sav'  ; the first rough overall map
save, filename=cfile, kont, pos

end
