pro iris_qs_display, inpath

;+
;===============================================================  
; procedure : iris_qs_display.pro
;  
; purpose: to create a Jpeg file showing the average rofile.
;          It should be run after iris_dispersion.pro
;  
; 2016/10/07 : created
;  
; R.Rezaei @ IAC                         e-mail:  rrezaei@iac.es      
;===============================================================
;-

print, inpath
files_d = read_dir(inpath , filter='aver_prof*.sav')
num = n_elements(files_d.files)
if (num gt 1) then begin
   files_d = read_dir(inpath , filter='aver_prof*CII*.sav')
endif   

restore, inpath + files_d.files[0]

f_d = read_dir(inpath , filter='iris_*im.fits')
new_file = strsplit(f_d.files[0], '.,_',/extract)
infile= inpath + new_file[0]+'_'+new_file[1]+'_'+new_file[2]+'_'+new_file[3]+'_'+'avprof.sav'
outfile= inpath + new_file[0]+'_'+new_file[1]+'_'+new_file[2]+'_'+new_file[3]+'_'+'avprof.jpg'

h = check_if_file_exists(infile)
if (string2num(h[1]) ne 1.) then begin 
   iris_dispersion, inpath
endif

print
print, outfile
print, '------------------'
restore, infile

plot, avprof, /xstyle, charsize=1.6, xthick=1.5, ythick=1.5, thick=1.5
xyouts, 0.2, 0.8, files_d.files[0], /normal, charsize=2, charthick=1.2
write_jpeg, outfile, tvrd(), quality=100


end  
