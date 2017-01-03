pro iris_benchmark, inpath, do_si=do_si, do_cii=do_cii, do_cl=do_cl, do_o=do_o, do_mg=do_mg, do_gauss=do_gauss

;+
;===============================================================
; procedure :  iris_benchmark.pro
;  
; purpose : to run the full analysis program only for the average profile.
;
; It is arbitrary to run it prior to the main analysis. No other
; program  depends on this program. It is just for visualization.  
;  
; May 09, 2016 : created. currently only works for one single average profile.
; May 19, 2016 : improved documentation, bug fix
;  
; R.Rezaei @ IAC                         e-mail:  rrezaei@iac.es      
;===============================================================
;-

common share_disp, dispersion_fuv, dispersion_nuv
common share_line, the_line                       

if n_params() lt 1 then begin
	print
	print, "usage:  iris_benchmark, /path/to/data/, keyword"
	print
	print, "	Analyze and plot IRIS line parameters of the average profile"
	print, "	different keys show different spectral lines"
	print
	print, "	/do_o, /do_cl, /do_cii, /do_si"
	print
	print, "	/do_mg with or without /do_gauss"
	print
        print, " e.g., iris_benchmark, '/data/obs/iris/spots/spot_20141024c/',/do_mg, /do_gauss "
	return
endif

if (n_elements(do_mg) eq 0) then do_mg=0 else do_mg=1
if (n_elements(do_gauss) eq 0) then do_gauss=0 else do_gauss=1
if (n_elements(do_cii) eq 0)  then do_cii=0  else do_cii=1
if (n_elements(do_o) eq 0)  then do_o=0  else do_o=1
if (n_elements(do_si) eq 0) then do_si=0 else do_si=1
if (n_elements(do_cl) eq 0) then do_cl=0 else do_cl=1

speed_of_light = 299792.458     ; km/s in vacuum

files_d = read_dir(inpath, filter='iris*_im.fits')
num_files = n_elements(files_d.files)
  
check = read_dir(inpath , filter='aver_prof_*.sav')
num = n_elements(check.files)
if ( num ge 1)and(strlen(check.files[0]) gt 5) then begin
   print, 'First calculate the average profile....'
   iris_average_profile, inpath
endif
check = read_dir(inpath , filter='aver_prof_*.sav')

restore, inpath + check.files[0];,/v

check = read_dir(inpath , filter='*avprof.sav')
print, check.files[0]
restore, inpath + check.files[0];,/v

dispersion_fuv = disp_fuv
dispersion_nuv = disp_nuv

;-------------------------------------------------------------------------------------
if (do_o eq 1) then begin
the_line = 'O'
if (num_files eq 1) then vv = avprof[OI_range[0]:OI_range[1]] else vv = avprof
if (min(vv) lt 0.) then vv -= (min(vv)) -1.
s = max( gauss_smooth(vv, 3, /edge_truncate), posmax)
master_o1 = posmax

del_o1 = ((139400.0d * 15.0 / speed_of_light)/ dispersion_fuv) 
np = n_elements(vv)

window, 10, xs=700, ys=550, title='O I profile'
erg = analyze_occ(vv, 0, np-1, master_o1, dispersion_fuv, 1)
endif
;-------------------------------------------------------------------------------------


;-------------------------------------------------------------------------------------
if (do_cl eq 1) then begin
the_line = 'Cl'
if (num_files eq 1) then vv = avprof[ClI_range[0]:ClI_range[1]] else vv = avprof
x = max(vv, posmax)
master_si = posmax
lpff, vv[master_si-9: master_si+9], ppos
master_cl = ppos + master_si-9

pos1 = (master_cl - 150.) > 0.
pos2 = (master_cl + 150.) < (ClI_range[1] - ClI_range[0])

erg = perform_single_gaussian(vv, pos1, pos2, master_cl, dispersion_fuv, 1)
endif
;-------------------------------------------------------------------------------------



if (do_si eq 1) then begin
;-------------------------------------------------------------------------------------
the_line = 'Si'   
if (num_files eq 1) then vv = avprof[SiIV_range1[0]:SiIV_range1[1]] else vv = avprof
guess_si_pos, inpath, avprof, dispersion_fuv, SiIV_range1, SiIV_range2
restore, inpath+'si_lines_position.sav'

ergm = perform_double_gaussian('Si', vv, 0, n1-1, master_si1, dispersion_fuv, 1)
;-------------------------------------------------------------------------------------
ans=''   &   read, ans, prompt='Press enter to continue ...'
;-------------------------------------------------------------------------------------
if (num_files eq 1) then vv = avprof[SiIV_range2[0]:SiIV_range2[1]] else vv = avprof
guess_si_pos, inpath, avprof, dispersion_fuv, SiIV_range1, SiIV_range2
restore, inpath+'si_lines_position.sav'

ergm = analyze_si(vv, 0, n2-1, master_si2, dispersion_fuv, 1)
;-------------------------------------------------------------------------------------

endif


;-------------------------------------------------------------------------------------
if (do_cii eq 1) then begin
the_line = 'C'
if (num_files eq 1) then vv = avprof[CII_range[0]:CII_range[1]] else vv = avprof
guess_cii_pos, inpath, vv, dispersion_fuv
restore, inpath+'cii_lines_position.sav'

ergm = analyze_c2(vv, ca_ilo[1], ca_ihi[1], master_c2, dispersion_fuv, 1)
endif
;-------------------------------------------------------------------------------------

if (do_mg eq 1) then begin
;-------------------------------------------------------------------------------------
the_line = 'Mg'
if (num_files eq 1) then vv = avprof[Mg_range[0]:Mg_range[1]] else vv = avprof
guess_mg_pos, inpath, vv, dispersion_nuv, /do_h_line
restore, inpath+'mg_lines_position.sav'
master = [master_k1v1, master_k1v2, master_k1r1, master_k1r2]

if (num_files eq 1) then vv = avprof[Mg_range[0]:Mg_range[1]] else vv = avprof
wing = mean(vv[60:64])
erg =  analyze_mg(vv, wing, ca_ilo, ca_ihi, master, dispersion_nuv, plt=9, do_gauss) ;

stop

if (num_files eq 1) then vv = avprof[Mg_range[0]:Mg_range[1]] else vv = avprof
wing = mean(vv[60:64])
master = [master_h1v1, master_h1v2, master_h1r1, master_h1r2]
erg =  analyze_mg(vv, wing, ca_ilo, ca_ihi, master, dispersion_nuv, plt=1, /do_H_line, do_gauss)
;-------------------------------------------------------------------------------------
endif


end  
