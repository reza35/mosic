pro iris_average_area, inpath

;+
;===============================================================
; function : iris_average_area.pro
;
; purpose : to create average profile of a selected area
; the user marks area of three regions (e.g., quiet sun, umbra, ..)
; and then the prgram calculated two separate average profile for each region.
;  
; 16 June 2016 : created
;
; R.Rezaei @ IAC                         e-mail:  rrezaei@iac.es      
;===============================================================
;-

on_error, 2

restore, inpath + 'iris_cont_analyzed_scan_0.sav'

tay = size(cont)
window, xs=tay[1], ys=tay[2]
tvscl, cont
print
print, 'mark an area for the first region'
q1 = defroi(tay[1], tay[2])
print
print, 'mark an area for the second region'
q2 = defroi(tay[1], tay[2])
print
print, 'mark an area for the third region'
q3 = defroi(tay[1], tay[2])
print

mask = cont - cont
mask(q1) = 1
mask(q2) = 2
mask(q3) = 3


files_d = read_dir(inpath , filter='iris*_im.fits')
data = inpath + files_d.files[0]
a=readfits(data, hed, nslice=0)
av = size(a)
prof1 = dblarr(av[1])
prof2 = dblarr(av[1])
prof3 = dblarr(av[1])

for kapa=0, tay[1]-1 do begin
  
  a=readfits(data, hed, nslice=kapa, /silent)
  a(where(~finite(a)))=-1.
  for j = 0, tay[2]-1 do begin
     case mask[kapa,j] of
        1: prof1 += reform(a[*,j])
        2: prof2 += reform(a[*,j])
        3: prof3 += reform(a[*,j])
        0: x=0
     endcase   
  endfor
  print, tay[1]-1-kapa
endfor

c1 = float(n_elements(q1))
c2 = float(n_elements(q2))
c3 = float(n_elements(q3))

prof1 /= c1
prof2 /= c2
prof3 /= c3
print
print, 'statistics of profiles:'
print, '#1  : ', round(c1)
print, '#2  : ', round(c2)
print, '#3  : ', round(c3)

save, filename=inpath+'area_av_prof.sav', prof1, prof2, prof3, c1, c2, c3, a

end
