pro analyze_iris, filepath, do_mg=do_mg, do_gauss=do_gauss, do_cii=do_cii, do_si=do_si,do_o1=do_o1,do_cl=do_cl,do_cont=do_cont, do_h=do_h, zfac=zfac, fast=fast, do_1394=do_1394, do_plot=do_plot, do_quiet=do_quiet

;+
;===============================================================
; procedure : analyze_iris
; purpose : to analyze IRIS data (all channels)
;
; May 29, 2014 : Mg II profile analysis
;
; Nov 03, 2014 : bug fix for very dark umbra;  
;
; Dec 09, 2014 : add Gaussian fit to Mg II analysis
;                some bug fix in Mg II analysis
;                improved single & double Gaussian fit to all FUV lines
;
; Dec 14, 2014 : improved analysis of Mg lines for different spectral sampling
;                we keep the normalization factor for each profile in mg_wing array
;                bug fix for rest wavelength
;
;
; June 19, 2015 : improved documentation about Mg II lines
;                 analyze_occ.pro: triple Gaussian fit for O I 1356 and the two Cl I lines
;
; Dec 25, 2015 : improved modification for Si IV multi line fit in presence of 
;                different dispersions. Many hardwired numbers in my_tgf.pro, triplegaussian.pro, 
;                and analyze_si,pro should be changed.
;
; Dec 27, 2015 : new routines to calculate spatial and spectral range
;                in each data and also to find the correct dispersion in each channel  
;                new routine to evaluate weak emission lines in the wing of Mg II k 
;
; Dec 28, 2015 : new procedure to calculate a safe mean
;
; Dec 30, 2015 : new routine to find line positions of Mg II range
;
; Apr 20, 2016 : new routines to find line positions of Si IV and C II,
;                bug fix for continuum map of the O I line
;
; May 31, 2016 : program works on off-limb data
;  
; Jun 01, 2016 : dummy O I velocity gradient, as sometimes this channel is missing in the data
;
; Jan 27, 2017 : keep 280 nm continuum intensity in DN units
;
; Aug 21, 2017 : bug fix for the dummy orbital velocity.  
;  
; Nov 20, 2017 : more control of the random steps.
;  
; Nov 29, 2017 : limb flag  
;  
; Dec 04, 2017 : the output filename keeps track of the input filename (date and time)
;
; Dec 11, 2017 : creates new overview Jpeg output from the Mg II data
;  
; R.Rezaei @ IAC                         e-mail:  rrezaei@iac.es      
;===============================================================
;-
common share_disp, dispersion_fuv, dispersion_nuv
common share_timeseri, do_ts
common share_ergm, best_fit
common share_line, the_line                       

if n_params() lt 1 then begin
	print
	print, "usage:  analyze_iris, /path/to/data/, /keyword"
	print
	print, "	This is the primary IRIS data analysis program."
	print, "	The users mainly needs to interact with this program and it communicates "
	print, "	with all subroutines internally."
	print
	print, "        The fits file should have the following structure: "
	print, "        wavelength as first axis, slit as second axis, and scan as the third axis."
        print, "        If you have IRIS Level3 data, you can re-arrange it to this shape using the following command:"
        print, "        .r reorder"
        print, "        reorder, '/full/path/to/your/file/.._all_im.fits'"
        print
        print, "        it moves the original file to a new file by adding _orig and then creates the desired data cube."
        print, "        The main reason to do so is that readfits can only cut one slice in the last axis."
	print
	print, "        Spectral line keywords:"
	print, "	/do_cont   two photospheric lines in 280 nm range"
        print, "	/do_mg     the Mg II k line"
        print, "	/do_h      the Mg II h line"
        print, "	/do_gauss   can be added to /do_mg to perform the Gaussian fitting for the lines."
	print
	print
        print, "        The following are all multi-Gaussian fitting procedures to any of the " 
        print, "        IRIS FUV lines."
        print
	print, "        /do_o1     the  O I and the two C I lines at 156 nm range"
	print, "        /do_cl     the  Cl I line at 156 nm range"
	print, "        /do_1394   the Si IV 139.4 nm line" 
        print, "        /do_si     the Si IV 140.3 nm line, three O IV lines, and one S IV line"
        print
	print, "        /do_cii    the C II 133 nm line pair" 
        print, "	/fast      can be added to /do_cii to skip fitting a penta-Gaussian fit."
	print
        print, "	The analysis of a data set should be started by running the /do_cont."
        print, "        The program checks several parameters and interactively ask the user some questions." 
	print, "        The range of spectral lines are determined by the user in this step."
	print, "        It is advisable o keep as much continuum as you have since the UV continuum is very weak."
	print
        print, "        The O I analysis should be performed prior to any other FUV data as we use the position of this line "
	print, "        to get the orbital velocity of the spacecraft."
	print
	print, "        Other keywords:"
	print, "	/zfac         the inverse of the zoom factor (the larger the zfac, the smaller the displayed image)"
        print, "	/do_plot      to plot the results for each individual profile"
	print
	print
        print, " e.g., , analyze_iris, '/data/obs/iris/spots/spot_20141024c/',/do_cont, zfac=1 "
        print
        print, "         At this point, you have the average profile of the data, spectral ranges, etc."
        print, "         So it is a good idea to first run the code for your average profile."
        print, "         check documentations of iris_benchmark.pro"
        print, "         If the fits are alright, then you can proceed to run for the whole data."
        print
        print, "         analyze_iris, '/data/obs/iris/spots/spot_20141024c/',/do_mg, /do_gauss, zfac=2, /do_h "
        print
        print, "         analyze_iris, '/data/obs/iris/spots/spot_20141024c/',/do_o1 "
        print
        print, "         analyze_iris, '/data/obs/iris/spots/spot_20141024c/',/do_cl "
        print
        print, "         analyze_iris, '/data/obs/iris/spots/spot_20141024c/',/do_cii, /fast, zfac=2"
        print
        print, "         analyze_iris, '/data/obs/iris/spots/spot_20141024c/',/do_si, /do_1394, zfac=2"
        print
        print, "         More than one line keyword can be set but make sure RAM of your system allows to do so."
        print, "         Gaussian fitting of big maps create output files >~ 1GB."
        return
endif

if (n_elements(do_mg) eq 0)  then do_mg=0  else do_mg=1
if (n_elements(do_gauss) eq 0)  then do_gauss=0  else do_gauss=1  ; perform a Gaussian fit for Mg II h/k
if (n_elements(do_h) eq 0)  then do_h=0  else if (do_mg eq 1) then do_h=1

if (n_elements(do_cii) eq 0) then do_cii=0 else do_cii=1
if (n_elements(do_o1) eq 0)  then do_o1=0  else do_o1=1
if (n_elements(do_si) eq 0)  then do_si=0  else do_si=1
if (n_elements(do_1394) eq 0)  then do_1394=0  else do_1394=1
if (n_elements(do_cl) eq 0)  then do_cl=0  else do_cl=1
if (n_elements(do_cont) eq 0) then do_cont=0 else do_cont=1
if (n_elements(fast) eq 0)  then fast=0  else fast=1
if (n_elements(do_plot) eq 0)  then do_plot=0

if (n_elements(do_quiet) eq 0)  then do_quiet=0  else do_quiet=1
if (do_quiet eq 1) then do_plot=0

if (n_elements(zfac) eq 0)  then zfac=1  else zfac=2
; the larger the zoom factor, the smaller the diplayed image <<

speed_of_light = 299792.458 ; km/s in vacuum
;---------------------------------------------------
; create an average profile from the data
;---------------------------------------------------
iris_average_profile, filepath

inpath = filepath
files_d = read_dir(inpath, filter='iris_l3*_im.fits')

num_files = n_elements(files_d.files)
;------------------------------------------------------
; set name of input files, if there are more than one
;------------------------------------------------------
if (num_files gt 1) then begin  ; open the right file, if different bands are stored in different files
   if (do_cont eq 1) then files_d = read_dir(inpath, filter='iris*_2832_im.fits')
   if (do_mg eq 1) then files_d = read_dir(inpath, filter='iris*_Mg*_im.fits')
   if (do_cii eq 1) then files_d = read_dir(inpath, filter='iris*_CII*_im.fits')
   if (do_si eq 1) then files_d = read_dir(inpath, filter='iris*_SiIV*_im.fits')
   if (do_o1 eq 1) then files_d = read_dir(inpath, filter='iris*_OI*_im.fits') 
endif  

data = inpath + files_d.files[0]

;---------------------------------------------------
; create a quasi-continuum map
;---------------------------------------------------
avfile =  inpath + 'kont.sav'
h = check_if_file_exists(avfile)

if (string2num(h[1]) ne 1.) then begin
   if (num_files ne 1) then kont_iris, inpath,/only_mg
   if (num_files eq 1) then kont_iris, inpath,/normal
endif 
restore, inpath+'kont.sav'
;------------------------------
; range aling the slit
;------------------------------
z0 = pos[0]
z1 = pos[1]

hed = headfits(data)
data_size = [fxpar(hed, 'NAXIS3'), fxpar(hed, 'NAXIS2'), fxpar(hed, 'NAXIS1'), fxpar(hed, 'NAXIS4')]
print, 'observed location in heliographic coordinate:   x0    y0   '
x0 = fxpar(hed, 'XCEN')
y0 = fxpar(hed, 'YCEN')
print, x0, y0,format='(F8.1,F9.1,/)'
if (abs(x0) gt 850.)or(abs(y0) gt 850.) then at_limb=1 else at_limb=0

vtf = size(a)

for scannum = 0,0 do begin  ;####################   loop over all repeated scans   ##################
cycle=scannum
cc = strtrim(num2string(cycle),1)
;------------------------------------------------------
; set name of output files
;------------------------------------------------------
s = strsplit(files_d.files[0], '.,_',/extract)
name_term = inpath + s[0] +'_'+ s[1] +'_'+ s[2] +'_'+ s[3]+'_'

outfile_cont = name_term + 'cont_analyzed.sav'
outfile_cii  = name_term + 'cii_analyzed.sav'
outfile_si   = name_term + 'si_analyzed.sav'
outfile_oi   = name_term + 'oi_analyzed.sav'
outfile_cl   = name_term + 'cl_analyzed.sav'
outfile_mg   = name_term + 'mg_analyzed.sav'
outfile_mg_em = name_term + 'mg_em_analyzed.sav'
outfile_nuv_temp = name_term + 'NUV_temporal.sav'
outfile_fuv_temp = name_term + 'FUV_temporal.sav'
avfile = name_term +'avprof.sav'

;------------------------------
; the spectral dispersion
;------------------------------
h = check_if_file_exists(avfile)
if (string2num(h[1]) ne 1.) then begin 
   save, filename=avfile, avprof, hed
   iris_dispersion, inpath
endif

restore, avfile
 
dispersion_nuv = disp_nuv       ; spectral sampling in pm/px
dispersion_fuv = disp_fuv       ; spectral sampling in pm/px
print,'+++++++++++++++++++++++++++++++++++++++'
print, 'spectral dispersion NUV/FUV = ', dispersion_nuv, '     / ', dispersion_fuv, ' pm/px'
print, 'my estimated disp   NUV/FUV = ', disp4, '     / ', disp3
print,'+++++++++++++++++++++++++++++++++++++++'
print,'spectral ranges were selected !        '
print,'+++++++++++++++++++++++++++++++++++++++'
wait, 1.


nx = data_size[0]  ; &   ny = tay[2]
ny = pos[1]-pos[0]+1  ; analyze only the selected range along the slit

;------------------------------------------------
;--  create arrays for output variables
;------------------------------------------------
if (do_cont eq 1) then begin
  cont = fltarr(nx, ny)
  kont = fltarr(nx, ny)
  fe = fltarr(nx, ny, 2)
endif
if (do_mg eq 1) then begin
   kpar = fltarr(nx, ny, 6)    ; K line  parameters
   k_spar = fltarr(nx, ny, 9)   ; K line single Gauusian fit parameters
   k_dpar = fltarr(nx, ny, 15)  ; K line double Gauusian fit parameters
   k_tpar = fltarr(nx, ny, 21)  ; K line triple Gauusian fit parameters
   fe_par = fltarr(nx, ny, 2)   ; Fe I wavelength positions 
   fe_int = fltarr(nx, ny, 2)   ; Fe I core intensities
   type_k = fltarr(nx, ny)      ; profile type = [0,1,2,3,4,5,6]
   mg_bnd = fltarr(nx, ny, 30)  ; band intensity parameters
   mg_wing = fltarr(nx, ny)     ; normalization factor for each pixel ~ continuum
   xmins_k = fltarr(nx, ny, 8)  ; position/amplitude of mins
   xmaxs_k = fltarr(nx, ny, 6)  ; position/amplitude of maxs
   min1v_k = fltarr(nx, ny, 2)  ; position/amplitude of h1v minima
   min1r_k = fltarr(nx, ny, 2)  ; position/amplitude of h1r minima

   em_lines = fltarr(nx, ny, 12) ; cont /core / EQW of four weak lines that go into emission in flare
   em_vel = fltarr(nx, ny, 4)    ; line-core position of the lines
;-------------------------------------------------------------------
   if (do_h eq 1) then begin
      hpar = fltarr(nx, ny, 6) ; H line  parameters
      h_spar = fltarr(nx, ny, 9) ; H line single Gauusian fit parameters
      h_dpar = fltarr(nx, ny, 15) ; H line double Gauusian fit parameters
      h_tpar = fltarr(nx, ny, 21) ; KHline triple Gauusian fit parameters
      type_h = fltarr(nx, ny)     ; profile type = [0,1,2,3,4,5,6]
      xmins_h = fltarr(nx, ny, 8) ; position/amplitude of mins
      xmaxs_h = fltarr(nx, ny, 6) ; position/amplitude of maxs
      min1v_h = fltarr(nx, ny, 2) ; position/amplitude of h1v minima
      min1r_h = fltarr(nx, ny, 2) ; position/amplitude of h1r minima
;-------------------------------------------------------------------
   endif
endif

if (do_cii eq 1) then begin
;------------------------------------------------------------------------
   c_ii_kont =  fltarr(nx, ny)         ; C II continuum from profile
   cii_cont =  fltarr(nx, ny)          ; C II continuum from Gaussian fit
   cii_s = fltarr(nx, ny, 9)           ; C II 1336 A line single Gaussian fitting
   cii_d = fltarr(nx, ny, 15)          ; C II 1335+1335 A line single Gaussian fitting
   cii_q = fltarr(nx, ny, 21)          ; C II 1335+1336 A line double Gaussian fitting
   cii_p = fltarr(nx, ny, 27)          ; C Ii lines  +  Ni II
   cii_fit_gauss = fltarr(nx, ny, cii_range[1]-cii_range[0]+1, 5) ; fits
;------------------------------------------------------------------------
endif

if (do_1394 eq 1) then begin
;------------------------------------------------------------------------
  si_iv1 = fltarr(nx, ny, 9)   ; Si IV 1394 A line single Gaussian fitting
  si_iv1_v = fltarr(nx, ny, 15) ; Si IV 1394  A line double Gaussian fitting
  si_iv1_bnd = fltarr(nx, ny, 3) ; band intensity parameters
  si_iv_tot = fltarr(nx, ny)
  si1_fit_gauss = fltarr(nx, ny, SiIV_range1[1]-SiIV_range1[0]+1, 3) ; to keep results of single & double Gaussian fits
;------------------------------------------------------------------------
endif


if (do_si eq 1) then begin
;------------------------------------------------------------------------
  ; the Si 1403 A fit also includes three O IV and a S IV lines
  si2_fit_gauss = fltarr(nx, ny, SiIV_range2[1]-SiIV_range2[0]+1, 6) ; to keep results of all fits
  si_iv_tot = fltarr(nx, ny)
  
  si_iv2 = fltarr(nx, ny, 9)     ; Si IV 1403 A line single-Gaussian fit with 4 free parameters
  si_iv2_v = fltarr(nx, ny, 15)  ; Si IV 1403 A line double-Gaussian fit with 7 free parameters
  si_iv2_q = fltarr(nx, ny, 15)  ; Si IV 1403 A line penta-Gaussian fit with 7 free parameters
  si_iv2_z = fltarr(nx, ny, 21)  ; Si IV 1403 A line penta-Gaussain fit with 9 free parameters
  si_iv2_f = fltarr(nx, ny, 27)  ; Si IV 1403 A line hexa-Gaussain fit with 12 free parameters
  si_iv2_bnd = fltarr(nx, ny, 3) ; band intensity parameters
;------------------------------------------------------------------------
endif

if (do_o1 eq 1) then begin
;------------------------------------------------------------------------
   occ = fltarr(nx, ny, 9)       ; O I + C I + C I lines at 1565 A : single Gaussian fit
   occ_q = fltarr(nx, ny, 14)    ; O I + C I + C I lines at 1565 A : triple Gaussian fit

   occ_fit = fltarr(nx, ny, OI_range[1]-OI_range[0]+1, 3) ; to keep results of single Gaussian fits  ; 106
   o_i_kont = fltarr(nx, ny)
   
   o_i1 = fltarr(nx, ny, 9)         ; O I 1356 A line  single Gaussian fitting
   oi_fit_gauss = fltarr(nx, ny, 41, 2) ; to keep results of single Gaussian fits
;------------------------------------------------------------------------
endif

if (do_cl eq 1) then begin
;------------------------------------------------------------------------
   cl = fltarr(nx, ny, 9)        ; Cl I 1352 A  single gaussian fit parameters
   cl_fit_gauss = fltarr(nx, ny, ClI_range[1]-ClI_range[0]+1, 2) ; to keep results of single Gaussian fits
;------------------------------------------------------------------------
endif






;------------------------------------------------------------------
;++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;------------------------ photosphere -----------------------------
;++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;------------------------------------------------------------------
if (do_cont eq 1) then begin

ca_ilo = [30, 78]  ; dummy variable
ca_ihi = [50, 88]  ; dummy variable
the_line = 'Fe'

ctd = float(ct_range[1] - ct_range[0]) + 1.  ; width of the continuum window

window, 10, xs=700, ys=550
if (num_files eq 1) then vv = avprof[ct_range[0]:ct_range[1]] else vv = avprof
print
print,'+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
print, 'click Carefully on the left/right borders of 2 photospheric blend.'
print, ' the first line should be the deeper one.'
print,'+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
print
read_plot_click, 4, vv, vv, pos_si
ca_ilo[0] = pos_si[0]   &  ca_ihi[0] = pos_si[1]
ca_ilo[1] = pos_si[2]   &  ca_ihi[1] = pos_si[3]

del_mg = ((298000.0d * 20.0 / speed_of_light)/ dispersion_nuv) ; velocities within 20 km/s
;--------------------------------------------------------------------
;- we use lpff.pro to find the line minimum.
;- lpff is not as accurate as a Gaussian fit but is faster. 
;--
;- in case of off-limb data, most of the continuum map returns empty.
;--------------------------------------------------------------------
for i=0, nx-1 do begin
   v = readfits(data, hed, nslice=i, /silent) 
   q = where(~finite(v), count)
   if (count gt 0) then v(q) = -10.

   for j=0, ny-1 do begin
      tmp = reform(v[ct_range[0]:ct_range[1], j + z0])

      kont[i,j] = max(median(tmp,5))
      cont[i,j] = good_mean(tmp[0:round(ctd/4.)])
      chk = where(tmp le (-10.0), count)
      ;-----------------------------
      if (median(tmp) gt 1.)and(count lt 70.) then begin   ;---- if it is a valid profile
         for k = 0, 1 do begin
           py = reform(tmp[ca_ilo[k]:ca_ihi[k]])  
           px = findgen(n_elements(py))
           lpff, py, pos_line
           if (~finite(pos_line)) then pos_line = (ca_ihi[k] - ca_ilo[k]) * 0.5
           fe[i,j,k] = pos_line + ca_ilo[k]
         endfor
      endif 
   endfor
   print, nx-1-i
endfor

aux = percentiles(kont/max(kont), value=[0.99, 1.])
kont = kont < (max(kont) * aux[0])
aux = percentiles(cont/max(cont), value=[0.995, 1.])
cont = cont < (max(cont) * aux[0])

;stdev_despike, cont, 10, 4
u = reform(fe[*,*,1])  &   stdev_despike, u, 10, 3   &  fe[*,*,1] = u
u = reform(fe[*,*,0])  &   stdev_despike, u, 10, 3   &  fe[*,*,0] = u

;--------------------------------------------
;-- maps of photospheric velocity/continuum
;--------------------------------------------
window, 0, xs=nx*3/zfac, ys=ny/zfac
s=median(fe[*,*,0])
tvsclm, fe[*,*,0] > (s-1.) < (s+1.), zm=zfac, cb=99
s=median(fe[*,*,1])
tvsclm, fe[*,*,1] > (s-1.) < (s+1.), xp=1*nx/zfac, yp=0, zm=zfac, cb=99
tvsclm, cont,      xp=2*nx/zfac, yp=0, zm=zfac
wait, 2.
;stop
s = percentiles(kont, value=[0.6, .92, .99])
q = where((kont gt s[1])and(kont lt s[2]))
i_cont = good_mean(cont(q))
cont /= i_cont
if (max(cont) gt 10.)or(max(cont lt 0.9)) then begin
   print, 'consider renormalization of the continuum map'
   ;stop
endif   
cont = cont < 10.0 > 0.
;-------------------------------------
;--  280 nm continuum properties
;-------------------------------------
print, 'Mean continuum counts [ADU] = ', i_cont
print, 'RMS contrast granulation [%] = ', stddev(cont(q))* 100.
print

;-------------------------------------------------
;-- systematic velocity residuals:
;-- 1) due to orbital motion of the satellite
;-- 2) due to spectral curvature along the slit
;-------------------------------------------------

;-------------------------------------------------------------------------
;-- now, we have to look for a suitable fit.
;-- we use the first line (the one which is stronger) as the primary data.
;-------------------------------------------------------------------------
curv =  total(fe[*,*,0],2)/float(ny)

;----------------------------------------------------------
;- replace NANs with an average value, if there is any
;- before we start to fit a curve to oribital velocity.
;----------------------------------------------------------
q = where(finite(curv))
qq = where(~finite(curv), count)
if (count gt 1) then curv(qq) = good_mean(curv)
;plot, curv
q = where(curv ne 0, count)
if (count lt 10) then stop
print, '++++++++++++++++++++++++++++++++++++++++++++++'
print, ' velocity gradient along the scan/time '
print, '++++++++++++++++++++++++++++++++++++++++++++++'
get_systematic_velocity, curv, output

temporal_gradient = output
save, filename=outfile_nuv_temp, temporal_gradient

;-------------------------------------------------
;-- apply correction for the temporal gradient
;-------------------------------------------------
x = reform(fe[*,*,0]) & q1 = where(x eq 0, count1)
if (count1 ge 1) then x(q1)= good_mean(x)   &   fe[*,*,0] = x
x = reform(fe[*,*,1]) & q2 = where(x eq 0, count2)
if (count2 ge 1) then x(q2)= good_mean(x)   &   fe[*,*,1] = x

for j = 0, ny-1 do fe[*,j,0] = fe[*,j,0] - temporal_gradient
for j = 0, ny-1 do fe[*,j,1] = fe[*,j,1] - temporal_gradient

;-------------------------------------------------
;-- normalize Fe velocities to average QS
;-------------------------------------------------
for i = 0, 1 do begin
  x = reform(fe[*,*,i])
  xm = good_mean(x[*,50:150])
  case i of 
    0: x = (x- xm)*dispersion_nuv * 3.0d+8 / 281450.0d
    1: x = (x- xm)*dispersion_nuv * 3.0d+8 / 281470.0d
    ;2: x = (x- xm)*dispersion_nuv * 3.0d+8 / 281500.0d
  endcase
  fe[*,0:ny-1,i] = x
endfor
;-------------------------------------------------
;-- remove the gradient along the slit
;-- a low order polynomial is enough.
;-------------------------------------------------
for i = 0, 1 do begin
  gr = reform(fe[*,*,i])
  gr = total(gr, 1, /double) / float(nx)
  tvsclm, reform(fe[*,*,i]), zm=2, cb=99
  print, '++++++++++++++++++++++++++++++++++++++++++++++'
  print, ' velocity gradient along the slit '
  print, '++++++++++++++++++++++++++++++++++++++++++++++'
  get_systematic_velocity, gr, yfit
  slit_gradient = yfit - mean(yfit)
  for k = 0, nx - 1 do fe[k,*,i] -= slit_gradient
  fe[*,*,i] = fe[*,*,i] - mean(fe[*,*,i]) - 300.   ; -300 m/s convective blueshift (?)
endfor


;----------------------------------------------------------------
; create a dummy file in case O I data is not present
; it will be overwritten, in case one runs the O I data reduction 
;-----------------------------------------------------------------
euv_temporal_gradient = temporal_gradient - temporal_gradient
save, filename=outfile_fuv_temp, euv_temporal_gradient

x = reform(fe[*,*,0])  &  if (count1 ge 1) then x(q1)= 0.  &   fe[*,*,0] = x
x = reform(fe[*,*,1])  &  if (count2 ge 1) then x(q2)= 0.  &   fe[*,*,1] = x

window, 4, xs = nx*3/zfac, ys = ny/zfac, title = 'phot. 280 nm'
tvsclm, cont < 1.9, zm=zfac
for i = 0, 1 do tvsclm,  fe[*,*,i] < 2000. > (-2000.), xp=nx*(i+1)/zfac, yp=0, zm=zfac, cb=99
loadct, 0, /silent

save, filename=outfile_cont, kont, cont, fe, i_cont, temporal_gradient

endif

; to keep the filenames backward compatible
h = check_if_file_exists(outfile_cont)
if (string2num(h[1]) ne 1.) then begin 
   old_filename = inpath + 'iris_cont_analyzed_scan_'+cc+'.sav'
   spawn, 'mv '+old_filename+'  '+outfile_cont
   
   old_filename = inpath + 'iris_NUV_temporal_'+cc+'.sav'
   spawn, 'mv '+old_filename+'  '+outfile_nuv_temp
   
   old_filename = inpath + 'iris_FUV_temporal_'+cc+'.sav'
   spawn, 'mv '+old_filename+'  '+outfile_fuv_temp
endif
restore, outfile_cont



;------------------------------------------------------------------
;++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;------------------------ O I line --------------------------------
;++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;------------------------------------------------------------------
if (do_o1 eq 1) then begin

the_line = 'O'
if (num_files eq 1) then vv = avprof[OI_range[0]:OI_range[1]] else vv = avprof

s = max( gauss_smooth(vv, 3., /edge_truncate), posmax)
master_o1 = posmax

del_o1 = ((139400.0d * 15.0 / speed_of_light)/ dispersion_fuv) ; = max delta lambda in px  
np = n_elements(vv)

tmp = max(vv, master_o1)
best_fit = [1., 4., master_o1,  2.]; continuum, amplitude, line center, line width
if (do_quiet eq 0) then begin
;  window, 10, xs=700, ys=550, title='O I profile'
  window, 8, xs=nx*7/zfac, ys=ny/zfac, title='O I 1356: continuum, velocity, amplitude, line width, amplitude, chi-sq, chi-sq'
endif
ergm = analyze_occ(vv, 0, np-1, master_o1, dispersion_fuv, 1)
best_fit = [ergm.spar1[0], ergm.spar1[2], ergm.spar1[1], ergm.spar1[3] > 2.]
;stop
;-----------------------------------------------------------------
;-- systematic velocity  for the FUV channel
;-- we take median of each slit and fit a polynomial to it.
;-----------------------------------------------------------------
h = check_if_file_exists(outfile_fuv_temp)
;if (string2num(h[1]) eq 0.) then begin
print
print, '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
print, 'We analyze average profile of each slit position to get the systematic velocity '
print, 'from a fit to the O I line and the two C I lines.'
print
print, 'If the O I data does not exist, skip running this step, or you get a crash ! '
print, '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
print
curve1 = fltarr(nx)
curve2 = fltarr(nx)
for i=0, nx-1 do begin
   v = readfits(data, hed, nslice=i, /silent) 
   q = where(~finite(v), count)
   if (count gt 0) then v(q) = -10.
   v = total(v[*, z0:z1], 2)/ float(ny)  & v = median(v,5)
   tmp = gauss_smooth(reform(v[OI_range[0]:OI_range[1]]),3, /edge_truncate)
   if (median(tmp) gt (-2.)) then begin
     erg = analyze_occ(tmp, 0, np-1, master_o1, dispersion_fuv, do_plot)
     curve1[i] = erg.spar1[1]
     curve2[i] = erg.qpar1[0]
   endif  
endfor

curve = median(curve1,3)

qq = where(~finite(curve), count)
if (count gt 1) then curve(qq) = good_mean(curve)
get_systematic_velocity, curve, output

euv_temporal_gradient = output
save, filename=outfile_fuv_temp, euv_temporal_gradient
print, outfile_fuv_temp, ' was created !'
erase
;endif
;-----------------------------------------------------------------

t_start = systime(/seconds)
for i=0, nx-1 do begin

   v = readfits(data, hed, nslice=i, /silent) 
   q = where(~finite(v), count)
   if (count gt 0) then v(q) = -10.
   v = despike_iris_raster(v, master_o1 + OI_range[0], master_o1 + OI_range[0])
 
   for j=0, ny-1 do begin
        d_stp = 50.0 / (dispersion_fuv / 1.298)
        t_stp = 15.0 / (dispersion_fuv / 1.298)
        tmp = gauss_smooth(reform(v[OI_range[0]:OI_range[1], j+z0]),1, /edge_truncate)
        o_i_kont[i,j] = good_mean(tmp[((master_o1 - d_stp)>0):(master_o1 - t_stp)])
        ;-----------------------------------------------------------------------
        if (stddev(tmp) gt 0.5) then begin   
        erg = analyze_occ(tmp, 0, np-1, master_o1, dispersion_fuv, do_plot)
            occ[i, j, *]   = erg.spar1    ; single Gaussian fit
            occ_q[i, j, *] = erg.qpar1    ; triple Gaussain fit

            occ_fit[i, j, *, 0] = erg.sprf
            occ_fit[i, j, *, 1] = erg.sfit 
            occ_fit[i, j, *, 2] = erg.qfit 
        endif
        ;-----------------------------------------------------------------------
   endfor
   ;print, nx-1 -i, FORMAT='(%"\b%d\b\b\b",$)'
   ;plot, occ[i, *, 4], psym=2 & loadct,40,/silent
   ;oplot, occ_q[i, *, 6], psym=1, color=245 & loadct,0,/silent
   ;stop
   if ((i mod 5) eq 0)and(i gt 0)and(do_quiet eq 0) then begin
     for kapa=0,3 do tvsclm, occ[0:i, *, kapa], xp=kapa*nx/zfac, zm=zfac
     tvsclm, occ_q[0:i, *, 1], xp=4*nx/zfac, zm=zfac  ; O I intensity
     tvsclm, occ[0:i, *, 4]<100.,   xp=5*nx/zfac, zm=zfac  ; chi square
     tvsclm, occ_q[0:i, *, 6]<100., xp=6*nx/zfac, zm=zfac  ; chi square
  endif
   if ((i mod 10) eq 0)and(i gt 0) then print, i, nx-1 -i, round( (systime(/seconds) - t_start)/float(i+1)* float(nx-1 -i) / 60.),' min'
endfor

;-----------------------------------------------
;-- now, we have to look for a suitable fit
;-----------------------------------------------
restore, outfile_fuv_temp
;-----------------------------------------------------------
;-- calculation of non-thermal width
;----------------------------------------------------
;- inputs: (1/e) width of the line profile in km/s

instr_fwhm = 0.028556 ; instrumental FWHM in A for short FUV
;----------------------------------------------------
w = reform(occ_q[*,*,2])  > 0.2 < 15. ; 1/e width in pixel
w = (w/2.) * (dispersion_fuv * 1.0d-3 /135.5598) * speed_of_light ; 1/e width in km/s
wnt_o_i1 = iris_nonthermalwidth('O','I', 1355.598, w, instr_fwhm, ti_max=ti_max,Wt_v=wt) ; non-thermal velocity in km/s
wt_o_i1 = wt
ti_o_i1 = ti_max
;-------------------------------------------------------------
w = reform(occ[*,*,3])  > 0.2 < 15. ; 1/e width in pixel
w = (w/2.) * (dispersion_fuv * 1.0d-3 /135.5598) * speed_of_light ; 1/e width in km/s
wnt_o_i2 = iris_nonthermalwidth('O','I', 1355.598, w, instr_fwhm, ti_max=ti_max,Wt_v=wt) ; non-thermal velocity in km/s
wt_o_i2 = wt
ti_o_i2 = ti_max
;-------------------------------------------------------------
w = reform(occ_q[*,*,4])  > 0.2 < 15. ; 1/e width in pixel
w = (w/2.) * (dispersion_fuv * 1.0d-3 /135.5840) * speed_of_light ; 1/e width in km/s
wnt_c_i1 = iris_nonthermalwidth('C','I', 1355.84, w, instr_fwhm, ti_max=ti_max,Wt_v=wt) ; non-thermal velocity in km/s
wt_c_i1 = wt
ti_c_i1 = ti_max
;-------------------------------------------------------------
o_i_vel =  reform(occ[*,*,1])
q1 = where(o_i_vel eq 0., cz)         ; zero pixels
if (at_limb eq 0) then for j=0, ny-1 do o_i_vel[*,j] = o_i_vel[*,j] - euv_temporal_gradient
o_i_vel -= ergm.spar1[1]
if (cz ge 1) then o_i_vel(q1) = 0.
o_i_vel = o_i_vel * dispersion_fuv * speed_of_light / 135559.8d + 1.000 ; km/s
;-----------------------------------------------------------------------------------
o_i_vel2 =  reform(occ_q[*,*,0])
q1 = where(o_i_vel2 eq 0., cz)         ; zero pixels
if (at_limb eq 0) then for j=0, ny-1 do o_i_vel2[*,j] = o_i_vel2[*,j] - euv_temporal_gradient
o_i_vel2 -= ergm.spar1[1]
if (cz ge 1) then o_i_vel2(q1) = 0.
o_i_vel2 = o_i_vel2 * dispersion_fuv * speed_of_light / 135559.8d + 1.000 ; km/s
;-----------------------------------------------------------------------------------
o_i_cont = reform(occ[*, *, 0])


if (do_quiet eq 0) then begin
  window, 13, xs=nx*6/zfac,ys=ny/zfac, title='O I 1356 A + C I 1356 A'
  tvsclm, o_i_vel < 7. > (-7.) , zm=zfac, cb=99; velocity
  tvsclm, o_i_cont,                  xp=nx/zfac,   zm=zfac ; continuum
  tvsclm, alog10(occ_q[*,*,1]>2.1),  xp=2*nx/zfac, zm=zfac ; log amplitude O I
  tvsclm, occ_q[*,*,2]<10. > 1.,     xp=3*nx/zfac, zm=zfac ; line width
  tvsclm, alog10(occ_q[*,*,3]> 1.5), xp=4*nx/zfac, zm=zfac ; amplitude C I 
  tvsclm, occ_q[*,*,3]<15. > 1.,     xp=5*nx/zfac, zm=zfac

  window, 14, title='reduced chi-square statistics, O I + C I fits'
  s=15.
  x= reform(occ[*, *, 4]) & q=where(x ne 0.) & plot_histogram, x(q), xrange=[0, s]
  loadct, 40,/silent
  x= reform(occ_q[*, *, 6]) & q=where(x ne 0.) & plot_histogram, x(q), xrange=[0, s], oplott=245
  loadct,0, /silent
endif 

x= reform(occ[*, *, 4])     &      q=where(x ne 0.)
r = percentiles(x(q), value=[0.05, 0.25, 0.5, 0.75, 0.95])
print, '-----------------------------------------------------'
print, 'Percentiles (singel Gaussian) ', r
x= reform(occ_q[*, *, 6])     &      q=where(x ne 0.)
r = percentiles(x(q), value=[0.05, 0.25, 0.5, 0.75, 0.95])
print, '-----------------------------------------------------'
print, 'Percentiles (triple Gaussian)', r
print, '-----------------------------------------------------'   

a = 0.
save, filename=outfile_oi, wnt_o_i1, wt_o_i1, ti_o_i1, o_i_vel, occ, occ_q, occ_fit, dispersion_fuv, $
      o_i_cont,   wnt_o_i2, wt_o_i2, ti_o_i2, wnt_c_i1, wt_c_i1, ti_c_i1, $
      master_o1, ergm, euv_temporal_gradient

print, 'O I: average run-time per slit= ',round((systime(/seconds) - t_start)) / float(nx), ' sec'
print, systime()

endif






;------------------------------------------------------------------
;++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;------------------------ Cl I line -------------------------------
;++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;------------------------------------------------------------------
if (do_cl eq 1) then begin

the_line = 'Cl'

if (num_files eq 1) then vv = avprof[ClI_range[0]:ClI_range[1]] else vv = avprof
s = max( gauss_smooth(vv, 3., /edge_truncate), posmax)
master_cl = posmax

ca_ilo = (master_cl - 150.) > 0.
ca_ihi = (master_cl + 150.) < (ClI_range[1] - ClI_range[0])


if (do_quiet eq 0) then begin
 ; window, 10, xs=700, ys=550, title='Cl I spectrum'
  window, 8, xs=nx*5/zfac, ys=ny/zfac, title='Cl I 1352 parameters'
endif
ergm = perform_single_gaussian(vv,ca_ilo,ca_ihi, master_cl, dispersion_fuv, 1)
best_fit = [ergm.spar[0], ergm.spar[2], ergm.spar[1], ergm.spar[3]]

t_start = systime(/seconds)
for i=0, nx-1 do begin 
   v = readfits(data, hed, nslice=i, /silent) 
   q = where(~finite(v), count)
   if (count gt 0) then v(q) = -10.
   v = despike_iris_raster(v, master_cl + ClI_range[0], master_cl + ClI_range[0])
 
   for j=0, ny-1 do begin
      tmp = gauss_smooth(reform(v[ClI_range[0]:ClI_range[1], j+z0]),1, /edge_truncate)
        if (max(tmp) gt (-1.))and(stddev(tmp) gt 0.5) then begin
          erg = perform_single_gaussian(tmp,ca_ilo,ca_ihi, master_cl, dispersion_fuv, 0)  

          cl[i, j, *]   = erg.spar
          cl_fit_gauss[i, j, *, 0] = erg.sprf
          cl_fit_gauss[i, j, *, 1] = erg.sfit
        endif
   endfor
   ;print, nx-1 -i, FORMAT='(%"\b%d\b\b\b",$)'
   if ((i mod 10) eq 0)and(i gt 0)and(do_quiet eq 0) then begin
     tvsclm, cl[0:i, *, 0], xp=0, zm=zfac
     tvsclm, cl[0:i, *, 1], xp=1*nx/zfac, zm=zfac
     tvsclm, cl[0:i, *, 2], xp=2*nx/zfac, zm=zfac
     tvsclm, cl[0:i, *, 3], xp=3*nx/zfac, zm=zfac
     tvsclm, cl[0:i, *, 4]<100., xp=4*nx/zfac, zm=zfac ; chi square
     print, i, nx-1 -i, round((systime(/seconds) - t_start)/(float(i+1)) * float(nx-1 -i) / 60.), ' min'
   endif
endfor

restore, outfile_fuv_temp
;----------------------------------------------------
;-- calculation of non-thermal width
;----------------------------------------------------
;- inputs: (1/e) width of the line profile in km/s

;instr_fwhm = 0.0318 ; instrumental FWHM in A for long FUV
instr_fwhm = 0.028556 ; instrumental FWHM in A for short FUV
;----------------------------------------------------
w = reform(cl[*,*,3])  < 20. ; 1/e width in pixel
w = (w/2.) * (dispersion_fuv * 1.0d-3 /135.1657) * speed_of_light ; 1/e width in km/s
wnt_cl = iris_nonthermalwidth('Cl','I', 1351.657, w, instr_fwhm, ti_max=ti_max,Wt_v=wt) ; non-thermal velocity in km/s
wt_cl = wt
ti_cl = ti_max
;-------------------------------------------------------------

cl_i_vel =  reform(cl[*,*,1])
q1 = where(cl_i_vel eq 0., cz)         ; zero pixels
if (at_limb eq 0) then for j=0, ny-1 do cl_i_vel[*,j] = cl_i_vel[*,j] - euv_temporal_gradient
cl_i_vel -= ergm.spar[1]
if (cz gt 0) then cl_i_vel(q1) = 0.
cl_i_vel = cl_i_vel * dispersion_fuv * speed_of_light / 135165.7d + 1.
;-----------------------------------------------------------------------------------
cl_cont = reform(cl[*, *, 0])

if (do_quiet eq 0) then begin
  window, 12, xs=nx*6/zfac,ys=ny*2/zfac, title='Cl I 1352 A: velocity, Log amplitude, FWHM, Log area, Log UV-cont, Log band H3, H2v, H2r'

  tvsclm, cl_i_vel < 15. > (-15.), zm=zfac, cb=99                ; velocity
  tvsclm, cl_cont,               xp=1*nx/zfac, yp=0, zm=zfac     ; continuum
  tvsclm, cl[*,*,3]<15.     ,    xp=2*nx/zfac, yp=0, zm=zfac     ; 1/e width
  tvsclm, alog10(cl[*,*,2]> 1.), xp=3*nx/zfac, yp=0, zm=zfac     ; log amplitude
  tvsclm, wnt_cl,                xp=4*nx/zfac, yp=ny/zfac, zm=zfac ; non-thermal width
  tvsclm, cl[*,*,4] < 10.,       xp=5*nx/zfac, yp=ny/zfac, zm=zfac ; chi^2 map

  window, 14, title='reduced chi-square statistics, Cl I fits'
  x= reform(cl[*, *, 4]) & q=where(x ne 0.) & plot_histogram, x(q), xrange=[0, 10]
endif 

x= reform(cl[*, *, 4])     &      q=where(x ne 0.)
r = percentiles(x(q), value=[0.05, 0.25, 0.5, 0.75, 0.95])
print, '-----------------------------------------------------'
print, 'Percentiles ', r
print, '-----------------------------------------------------'   

a = 0.
save, filename=outfile_cl, wnt_cl, wt_cl, ti_cl, cl, cl_fit_gauss, $
      cl_cont, master_cl, ergm, euv_temporal_gradient

print, 'Cl I: average run-time per slit= ',round((systime(/seconds) - t_start) / float(nx)), ' sec'
print, systime()

endif





;------------------------------------------------------------------
;++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;----------------------- Si IV line -------------------------------
;++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;------------------------------------------------------------------

if (do_si eq 1)or(do_1394) then begin

the_line = 'Si'   
;---------------------------
;- find the line positions
;---------------------------
guess_si_pos, inpath, avprof, dispersion_fuv, SiIV_range1, SiIV_range2
restore, inpath+'si_lines_position.sav'

if (num_files eq 1) then vv = avprof[SiIV_range2[0]:SiIV_range2[1]] else vv = avprof

if (do_quiet ne 1) then begin
;  window, 10, xs=700, ys=550, title='Si IV 1403 spectrum'
  window, 17, xs=nx*6/zfac, ys=ny/zfac, title='Si IV parameters: single Gaussian parameters'
  window, 11, xs=nx*8/zfac, ys=ny/zfac, title ='Si IV parameters: chi-sq in five different fits'
endif
ergm = analyze_si(vv, 0, n2-1, master_si2, dispersion_fuv, 1)
; it contains properties of the average profile and its fit
best_fit = [ergm.spar1[0], ergm.spar1[2], ergm.spar1[1], ergm.spar1[3]]

t_start = systime(/seconds)
for i=0, nx-1 do begin
   v = readfits(data, hed, nslice=i, /silent) 
   q = where(~finite(v), count)
   if (count gt 0) then v(q) = -10.
   v = despike_iris_raster(v, master_si1 + SiIV_range1[0], master_si2 + SiIV_range2[0])
   
   for j=0, ny-1 do begin
      ;print, i,j
      tmp0 = gauss_smooth(reform(v[*, j + z0]), 2., /edge_truncate)
      si_iv_tot[i, j] = good_mean(tmp0)

      if (do_1394 eq 1) then begin
      ;------------------------------------------------------
      ;- Si IV 1394
      ;------------------------------------------------------
        tmp = tmp0[SiIV_range1[0]:SiIV_range1[1]]
        if (max(tmp) ge (10.))and(stddev(tmp) gt 0.5) then begin 
          erg = perform_double_gaussian('Si', tmp, 0, n1-1, master_si1, dispersion_fuv, 0)

          si_iv1[i, j, *]   = erg.spar1
          si_iv1_v[i, j, *]   = erg.dpar1
          si_iv1_bnd[i, j, *] = erg.band1

          si1_fit_gauss[i, j, *, 0] = erg.sprf
          si1_fit_gauss[i, j, *, 1] = erg.sfit 
          si1_fit_gauss[i, j, *, 2] = erg.dfit 
        endif
     endif
     if (do_si eq 1) then begin
     ;------------------------------------------------------
     ;- Si IV 1403
     ;------------------------------------------------------
       tmp = tmp0[SiIV_range2[0]:SiIV_range2[1]]
       if (stddev(tmp) gt (0.5)) then begin
          erg = analyze_si(tmp, 0, n2-1, master_si2, dispersion_fuv, do_plot) ;--- the single/double/penta/hexa-Gaussian fit

          si_iv2[i, j, *]   = erg.spar1        ; single Gaussian fit, 4 pars
          si_iv2_v[i, j, *] = erg.dpar1        ; double Gaussian fit, 7 pars
          si_iv2_q[i, j, *] = erg.qpar1        ; parameters of the penta-Gussian fit, 7 pars
          si_iv2_z[i, j, *] = erg.zpar1        ; parameters of the penta-Gussian fit, 9 pars
          si_iv2_f[i, j, *] = erg.fpar1        ; parameters of the penta-Gussian fit, 12 pars

          si_iv2_bnd[i, j, *] = erg.band1
        
          si2_fit_gauss[i, j, *, 0] = erg.sprf ; observed profile
          si2_fit_gauss[i, j, *, 1] = erg.sfit ; single Gaussian fit, 4 pars
          si2_fit_gauss[i, j, *, 2] = erg.dfit ; double Gaussian fit, 7 pars
          si2_fit_gauss[i, j, *, 3] = erg.qfit ; penta Gaussian fit, 7 pars
          si2_fit_gauss[i, j, *, 4] = erg.zfit ; penta Gaussian fit, 9 pars
          si2_fit_gauss[i, j, *, 5] = erg.ffit ; penta Gaussian fit, 12 pars
       endif
     endif
  endfor
  ;stop 
  ; print, nx-1 -i, FORMAT='(%"\b%d\b\b\b",$)'
  if ((i mod 5) eq 0)and(i gt 0)and(do_quiet eq 0) then begin
      wset, 17
      for ip=0,5 do tvsclm, reform(si_iv2[0:i,*,ip]), xp=nx * ip/zfac, zm=zfac

      wset, 11
      for ip=0,3 do tvsclm, reform(si_iv2_q[0:i,*,ip]), xp=nx*ip/zfac, zm=zfac
      tvsclm, si_iv2_v[0:i, *, 7],         xp=nx*4/zfac, yp=zfac, zm=zfac ; chi square
      tvsclm, si_iv2_q[0:i, *, 7],         xp=nx*5/zfac, yp=zfac, zm=zfac ; chi square
      tvsclm, si_iv2_z[0:i, *, 10],        xp=nx*6/zfac, yp=zfac, zm=zfac ; chi square
      tvsclm, si_iv2_f[0:i, *, 13],        xp=nx*7/zfac, yp=zfac, zm=zfac ; chi square
      ;stop
   endif
   if ((i mod 10) eq 0)and(i gt 0) then print, i, nx-1 -i, round((systime(/seconds) - t_start)/float(i+1) * float(nx-1 -i) / 60.),' min'

endfor

;---------------------------------------------------------------------------------
print, 'initial saving ....'
save, filename=outfile_si, fe, cont, si_iv2, si_iv2_v, si_iv2_q, si_iv2_z, si_iv2_f, si_iv2_bnd, si2_fit_gauss, $
      si_iv1, si_iv1_v, si_iv1_bnd, si1_fit_gauss, euv_temporal_gradient, si_kont

;-------------------------------------------------------------
;-- removing the systematic line-center shifts for all lines
;-------------------------------------------------------------
restore, outfile_fuv_temp

;----------------------------------------------------
;-- calculation of non-thermal line width
;----------------------------------------------------
;- inputs: (1/e) width of the line profile in km/s

;instr_fwhm = 0.0318 ; instrumental FWHM in A for long FUV
instr_fwhm = 0.028556 ; instrumental FWHM in A for short FUV
;----------------------------------------------------
if (do_1394 eq 1) then begin
;-------------------------------------------------------------
   w = reform(si_iv1[*,*,3])  > 2. < 60.            ; 1/e width in pixel                                        
   w = (w/2.) * (dispersion_fuv * 1.0d-3 /139.378) * speed_of_light    ; 1/e width in km/s
   wnt_si_iv1 = iris_nonthermalwidth('Si','IV', 1393.78, w, instr_fwhm, ti_max=ti_max,Wt_v=wt) ; non-thermal velocity in km/s
   wt_si_iv1 = wt
   ti_si_iv1 = ti_max
;----------------------------------------------------------------------------
   si_iv1_vel = reform(si_iv1[*,*,1])
   if (at_limb eq 0) then for j=0, ny-1 do si_iv1_vel[*,j] = si_iv1_vel[*,j] - euv_temporal_gradient
   si_iv1_vel = reform(si_iv1_vel) < (master_si1 + del_si) > (master_si1 - del_si)
   si_iv1_vel -= master_si1
   si_iv1_vel = si_iv1_vel * dispersion_fuv * speed_of_light / 139378.0d + 6.500
;----------------------------------------------------------------------------  
endif

;-------------------------------------------------------------
if (do_si eq 1) then begin
;-------------------------------------------------------------
   w = reform(si_iv2[*,*,3])  > 2. < 60.          ; 1/e width in pixel
   w = (w/2.) * (dispersion_fuv * 1.0d-3 /140.277) * speed_of_light    ; 1/e width in km/s
   wnt_si_iv2 = iris_nonthermalwidth('Si','IV', 1402.77, w, instr_fwhm, ti_max=ti_max,Wt_v=wt) ; non-thermal velocity in km/s
   wt_si_iv2 = wt
   ti_si_iv2 = ti_max
;--------------------------------------------------------------
   w = reform(si_iv2_f[*,*,2])  > 2. < 60.            ; first Gaussain component
   w = (w/2.) * (dispersion_fuv * 1.0d-3 /140.277) * speed_of_light     ; 1/e width in km/s
   wnt_si_iv2_f1 = iris_nonthermalwidth('Si','IV', 1402.77, w, instr_fwhm, ti_max=ti_max,Wt_v=wt) ; non-thermal velocity in km/s
   wt_si_iv2_f1 = wt
   ti_si_iv2_f1 = ti_max
;--------------------------------------------------------------
   w = reform(si_iv2_f[*,*,5])  > 2. < 60.               ; second Gaussain component
   w = (w/2.) * (dispersion_fuv * 1.0d-3 /140.277) * speed_of_light         ; 1/e width in km/s
   wnt_si_iv2_f2 = iris_nonthermalwidth('Si','IV', 1402.77, w, instr_fwhm, ti_max=ti_max,Wt_v=wt) ; non-thermal velocity in km/s
   wt_si_iv2_f2 = wt
   ti_si_iv2_f2 = ti_max
;--------------------------------------------------------------
   w = reform(si_iv2_f[*,*,8])  > 2. < 60.                                                  ; O IV 1401
   w = (w/2.) * (dispersion_fuv * 1.0d-3 /140.1156) * speed_of_light                                   ; 1/e width in km/s
   wnt_o_iv_f = iris_nonthermalwidth('O','IV', 1401.156, w, instr_fwhm, ti_max=ti_max,Wt_v=wt) ; non-thermal velocity in km/s
   wt_o_iv_f = wt
   ti_o_iv_f = ti_max
endif
;---------------------------------------------------------------
;-- calculation electron density using Si IV lines
;---------------------------------------------------------------
; Note 1: Si IV is NOT the best diagnostic for electron density !!
;       use OIV instead !
; Note 2: when the lines approach toward being optically thick, 
;         then the result is basically garbage !
;---------------------------------------------------------------
; input: integral line intensity (background subtracted)
;
; -- integral of a Gaussian is sigma x sqrt(pi) * amplitude
;
if (do_si eq 1)and(do_1394 eq 1) then begin
   int1 = reform(si_iv1[*, *, 2]) * reform(si_iv1[*, *, 3])
   int2 = reform(si_iv2[*, *, 2]) * reform(si_iv2[*, *, 3])
   err1 = 0.5 * sqrt(int1)
   err2 = 0.5 * sqrt(int2)
   denf = 'si4_1392to1404_den.sav'
   ned_si4 = fltarr(nx, ny, 3) 
   si4_ne = fltarr(ny)      & si4_ne1 = fltarr(ny)      & si4_ne2 = fltarr(ny)
   for i=0,nx-1 do begin
      iris_ne, reform(int1[i,*]), reform(err1[i,*]), reform(int2[i,*]), reform(err2[i,*]), denf, rat, rat_err, si4_ne, si4_ne1, si4_ne2
      ned_si4[i,*, 0] = si4_ne  ; n_e = electron density 
      ned_si4[i,*, 1] = si4_ne1 ; lower limit of n_e
      ned_si4[i,*, 2] = si4_ne2 ; upper limit of n_e
   endfor
endif

;---------------------------------------------------------------
;-- calculation of electron density using O IV lines
;---------------------------------------------------------------
;- inputs: intensity or integral intensity (background subtracted)
; -- note that since integral of a Gaussian is sigma x sqrt(pi) *
; -- amplitude, only the amplitude ratio matters in this calculations
; -- since sigma of the two lines are the same.
;---------------------------------------------------------------
if (do_si eq 1) then begin
   int1 = reform(si_iv2_f[*, *, 7]) ; O IV 1400; the fit amplitude
   int2 = reform(si_iv2_f[*, *, 9]) ; O IV 1401; the fit amplitude
   err1 = sqrt(int1) * 0.5 > 0.1
   err2 = sqrt(int2) * 0.5 > 0.1
   denf = 'o4_1399to1402_den.sav'
   ned_o4 = fltarr(nx, ny, 3)
   o4_ne = fltarr(ny)      & o4_ne1 = fltarr(ny)      & o4_ne2 = fltarr(ny)
   for i=0,nx-1 do begin
  ;iris_ne, int1[i,*], err1[i,*], int2[i,*], err2[i,*], denf, rat, rat_err, o4_ne, o4_ne1, o4_ne2
      iris_ne_oiv, int1[i,*], err1[i,*], int2[i,*], err2[i,*], rat, rat_err, o4_ne, o4_ne1, o4_ne2
      ned_o4[i,*, 0] = o4_ne
      ned_o4[i,*, 1] = o4_ne1   ; lower limit
      ned_o4[i,*, 2] = o4_ne2   ; upper limit
   endfor


si_iv2_vel = reform(si_iv2[*,*,1])
q1 = where(si_iv2_vel eq 0., cz) 
if (at_limb eq 0) then for j=0, ny-1 do si_iv2_vel[*,j] = si_iv2_vel[*,j] - euv_temporal_gradient   ; ---
si_iv2_vel = reform(si_iv2_vel) < (ergm.fpar1[0] + del_si) > (ergm.fpar1[0] - del_si)
si_iv2_vel -= ergm.fpar1[0]
si_iv2_vel = si_iv2_vel * dispersion_fuv * speed_of_light / 140276.9d + 4.80d
if (cz gt 0) then si_iv2_vel(q1) = 0.
;-------------------------------------------------------------------------------
si_iv2q_vel = reform(si_iv2_q[*,*,0])
q1 = where(si_iv2q_vel eq 0., cz) 
if (at_limb eq 0) then for j=0, ny-1 do si_iv2q_vel[*,j] = si_iv2q_vel[*,j] - euv_temporal_gradient
si_iv2q_vel = reform(si_iv2q_vel) < (ergm.fpar1[0] + del_si) > (ergm.fpar1[0] - del_si)
si_iv2q_vel -= ergm.qpar1[0]
si_iv2q_vel = si_iv2q_vel * dispersion_fuv * speed_of_light / 140276.9d + 4.80d
if (cz gt 0) then si_iv2q_vel(q1) = 0.
;-------------------------------------------------------------------------------
si_iv2f_vel = reform(si_iv2_f[*,*,0])
q1 = where(si_iv2f_vel eq 0., cz) 
if (at_limb eq 0) then for j=0, ny-1 do si_iv2f_vel[*,j] = si_iv2f_vel[*,j] - euv_temporal_gradient
si_iv2f_vel = reform(si_iv2f_vel) < (ergm.fpar1[0] + del_si) > (ergm.fpar1[0] - del_si)
si_iv2f_vel -= ergm.fpar1[0]
si_iv2f_vel = si_iv2f_vel * dispersion_fuv * speed_of_light / 140276.9d + 4.80d
if (cz gt 0) then si_iv2f_vel(q1) = 0.
;-------------------------------------------------------------------------------
o_iv_vel = reform(si_iv2_f[*,*,6])
q1 = where(o_iv_vel eq 0., cz) 
if (at_limb eq 0) then for j=0, ny-1 do o_iv_vel[*,j] = o_iv_vel[*,j] - euv_temporal_gradient
o_iv_vel = reform(o_iv_vel) < (ergm.fpar1[6] + del_si) > (ergm.fpar1[6] - del_si)
o_iv_vel -= ergm.fpar1[6]
o_iv_vel = o_iv_vel * dispersion_fuv * speed_of_light / 140115.6d + 12.0d
if (cz gt 0) then o_iv_vel(q1) = 0.
;-------------------------------------------------------------------------------

endif

if (do_quiet eq 0) then begin
if (nx gt 400) then begin
 window, 1, xs=nx*5/zfac,ys=ny*2/zfac, title='Si IV 1394/1403 A: velocity, Log amplitude, FWHM, Log area, Log UV-cont, Log band H3, Wnt, ..'

 tvsclm, si_iv2[*,*,0],                 xp=0*nx/zfac, yp=0, zm=zfac         ; continuum
  tvsclm, si_iv2_vel < 2d1 >(-2d1),xp=1*nx/zfac, yp=0, zm=zfac, cb=99; velocity
  tvsclm, si_iv2[*,*,2],                xp=2*nx/zfac, yp=0, zm=zfac ; log amplitude
  tvsclm, si_iv2[*,*,1],                xp=3*nx/zfac, yp=0, zm=zfac ; line center
  tvsclm, si_iv2[*,*,3]> 0.,            xp=4*nx/zfac, yp=0, zm=zfac ; line width
  
  tvsclm, wnt_si_iv2,        xp=0*nx/zfac, yp=ny/zfac, zm=zfac
  tvsclm, si_iv2[*,*,4],     xp=1*nx/zfac, yp=ny/zfac, zm=zfac ; chi-square, single Gaussian
  tvsclm, si_iv2_v[*,*,7],   xp=2*nx/zfac, yp=ny/zfac, zm=zfac ; chi-square, double Gaussian
  tvsclm, si_iv2_q[*,*,7],   xp=3*nx/zfac, yp=ny/zfac, zm=zfac ; chi-square, penta Gaussian
  tvsclm, si_iv2_f[*,*,12],  xp=4*nx/zfac, yp=ny/zfac, zm=zfac ; chi-square, hexaa Gaussian
   
endif else begin

  window, 1, xs=nx*15/zfac,ys=ny/zfac, title='Si IV 1394/1403 A: cont, vel, Log amp, FWHM, Log area, Log UV-cont, Log band H3, Wnt, ..'

  tvsclm, si_iv2[*,*,0],             xp=0*nx/zfac, yp=0, zm=zfac ; continuum
  tvsclm, si_iv2_vel < 2d1 >(-2d1),  xp=1*nx/zfac, yp=0, zm=zfac, cb=99; velocity Si IV
  tvsclm, o_iv_vel < 2d1 >(-2d1),    xp=2*nx/zfac, yp=0, zm=zfac, cb=99; velocity O IV
  tvsclm, si_iv2[*,*,2],             xp=3*nx/zfac, yp=0, zm=zfac   ; log amplitude
  tvsclm, si_iv2_q[*,*,3],           xp=4*nx/zfac, yp=0, zm=zfac ; line center
  tvsclm, wnt_si_iv2,                xp=5*nx/zfac, yp=0, zm=zfac ; line width
  tvsclm, wnt_si_iv2_f1,             xp=6*nx/zfac, yp=0, zm=zfac
  tvsclm, wnt_si_iv2_f2,             xp=7*nx/zfac, yp=0, zm=zfac
  tvsclm, wnt_o_iv_f,                xp=8*nx/zfac, yp=0, zm=zfac
  

  tvsclm, si_iv2[*,*,4],     xp=9*nx/zfac, yp=0, zm=zfac  ; chi-square, single Gaussian
  tvsclm, si_iv2_v[*,*,7],   xp=10*nx/zfac, yp=0, zm=zfac  ; chi-square, double Gaussian
  tvsclm, si_iv2_q[*,*,7],   xp=11*nx/zfac, yp=0, zm=zfac  ; chi-square, penta Gaussian
  tvsclm, si_iv2_z[*,*,7],   xp=12*nx/zfac, yp=0, zm=zfac  ; chi-square, penta Gaussian
  tvsclm, si_iv2_f[*,*,12],  xp=13*nx/zfac, yp=0, zm=zfac ; chi-square, hexaa Gaussian
endelse 

window, 14, title='reduced chi-square statistics, Si IV + O IV fits'
s = 5
  x= reform(si_iv2[*, *, 4]) & q=where(x ne 0.) & plot_histogram, x(q), xrange=[0, s]
  loadct,40,/silent
  x= reform(si_iv2_v[*, *, 7]) & q=where(x ne 0.) & plot_histogram, x(q), xrange=[0, s], oplott=70
  x= reform(si_iv2_q[*, *, 7]) & q=where(x ne 0.) & plot_histogram, x(q), xrange=[0, s], oplott=170
  x= reform(si_iv2_z[*, *, 10]) & q=where(x ne 0.) & plot_histogram, x(q), xrange=[0, s], oplott=200
  x= reform(si_iv2_f[*, *, 13]) & q=where(x ne 0.) & plot_histogram, x(q), xrange=[0, s], oplott=250
  loadct,0,/silent
  
endif 

x= reform(si_iv2[*, *, 4])     &      q=where(x ne 0.)
r = percentiles(x(q), value=[0.05, 0.25, 0.5, 0.75, 0.95])
print, '-------------------------------------------------------------------------------'
print, 'Percentiles (single Gaussian)', r
x= reform(si_iv2_v[*, *, 7])     &      q=where(x ne 0.)
r = percentiles(x(q), value=[0.05, 0.25, 0.5, 0.75, 0.95])
print, '-------------------------------------------------------------------------------'
print, 'Percentiles (double Gaussian)', r
x= reform(si_iv2_q[*, *, 7])     &      q=where(x ne 0.)
r = percentiles(x(q), value=[0.05, 0.25, 0.5, 0.75, 0.95])
print, '-------------------------------------------------------------------------------'
print, 'Percentiles (penta Gaussian I)', r
x= reform(si_iv2_z[*, *, 10])     &      q=where(x ne 0.)
r = percentiles(x(q), value=[0.05, 0.25, 0.5, 0.75, 0.95])
print, '-------------------------------------------------------------------------------'
print, 'Percentiles (penta Gaussian II)', r
x= reform(si_iv2_f[*, *, 13])     &      q=where(x ne 0.)
r = percentiles(x(q), value=[0.05, 0.25, 0.5, 0.75, 0.95])
print, '-------------------------------------------------------------------------------'
print, 'Percentiles (hexa Gaussian)', r
print, '-------------------------------------------------------------------------------'

a = 0.
print, 'saving ....'
if (do_1394 eq 1) then begin
save, filename=outfile_si, fe, cont, si_iv1, si_iv1_v,  dispersion_fuv,$
      si_iv1_bnd, si_iv2, si_iv2_v, si_iv2_q, si_iv2_z, si_iv2_f, si_iv2_bnd, si_kont, $
      si_iv_tot, si1_fit_gauss, si2_fit_gauss, euv_temporal_gradient, master_si1, master_si2, $
      ned_o4, ned_si4, ergm, si_iv1_vel, si_iv2_vel, si_iv2f_vel, si_iv2q_vel, o_iv_vel,$
      wt_si_iv1, ti_si_iv1, wnt_si_iv1, $
      wt_si_iv2, ti_si_iv2, wnt_si_iv2, $
      wt_si_iv2_f1, ti_si_iv2_f1, wnt_si_iv2_f1, si_iv2q_vel, o_iv_vel, $
      wt_si_iv2_f2, ti_si_iv2_f2, wnt_si_iv2_f2, wnt_o_iv_f
endif else begin
save, filename=outfile_si, fe, cont, $
      si_iv2, si_iv2_v, si_iv2_q, si_iv2_z, si_iv2_f, si_iv2_bnd, dispersion_fuv,$
      si_iv_tot, si2_fit_gauss, euv_temporal_gradient, master_si1, master_si2, $
      ned_o4, ergm, si_iv2_vel, si_iv2f_vel, si_iv2q_vel, o_iv_vel,$
      wt_si_iv2, ti_si_iv2, wnt_si_iv2, $
      wt_si_iv2_f1, ti_si_iv2_f1, wnt_si_iv2_f1, si_iv2q_vel, o_iv_vel, $
      wt_si_iv2_f2, ti_si_iv2_f2, wnt_si_iv2_f2, wnt_o_iv_f
endelse

print, 'average run-time per slit= ',round((systime(/seconds) - t_start) / float(nx)), ' sec'
print, systime()
endif






;------------------------------------------------------------------
;++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;----------------------- C II line --------------------------------
;++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;------------------------------------------------------------------
if (do_cii eq 1) then begin

the_line = 'C'
if (num_files eq 1) then vv = avprof[CII_range[0]:CII_range[1]] else vv = avprof
;---------------------------
;- find the line positions
;---------------------------
guess_cii_pos, inpath, vv, dispersion_fuv
restore, inpath+'cii_lines_position.sav'

if (do_quiet eq 0) then begin
  window, 17, xs=nx*5/zfac, ys=ny/zfac, title='C II: continuum, amplitude, ...'
  window, 10, xs=700, ys=550, title='C II and Ni II profiles'
  ergm = analyze_c2(vv, ca_ilo[1], ca_ihi[1], master_c2, dispersion_fuv, 1, 0) 
  window, 11, xs=nx*6/zfac, ys=ny/zfac, title='C II and Ni II: chi squares, ...'
endif else begin
  ergm = analyze_c2(vv, ca_ilo[1], ca_ihi[1], master_c2, dispersion_fuv, 1, 0)   
endelse
best_fit = [ergm.spar1[0], ergm.spar1[2], ergm.spar1[1], ergm.spar1[3]]

t_start = systime(/seconds)
for i=0, nx-1 do begin
   ;v = reform(a[i,*, CII_range[0]:CII_range[1], cycle])
   ;v = iris_despike_raster(v, master_c1, master_c2)

   v = readfits(data, hed, nslice=i, /silent) 
   q = where(~finite(v), count)
   if (count gt 0) then v(q) = -10.
   v = despike_iris_raster(v, master_c1 + CII_range[0], master_c2 + CII_range[0])
   
   for j=0, ny-1 do begin
        ;print, i, j
        tmp = reform(v[CII_range[0]:CII_range[1], j+z0])
        tmp = gauss_smooth(tmp,1, /edge_truncate)
        c_ii_kont[i,j] = good_mean(tmp[((master_c1 - 80)>0):(master_c1 - 15)])
        ;if (max(tmp) gt (5.))and(stddev(tmp) gt (0.5)) then begin
        if (stddev(tmp) gt (0.5)) then begin
        ;------------------------------------------------------
        ;- C II 133.45 + C II 133.57 + Ni II 133.52
        ;------------------------------------------------------
          erg = analyze_c2(tmp, ca_ilo[1], ca_ihi[1], master_c2, dispersion_fuv, do_plot, fast)
          cii_s[i, j, *]   = erg.spar1 ; one single Gaussian to the two C II 133.57 nm line
          cii_d[i, j, *]   = erg.dpar1 ; two single Gaussian to the two C II lines
          cii_q[i, j, *]   = erg.qpar1 ; quad Gaaussian fit
          cii_p[i, j, *]   = erg.ppar1 ; penta Gaaussian fit
          cii_fit_gauss[i, j, *, 0] = erg.sprf
          cii_fit_gauss[i, j, *, 1] = erg.sfit
          cii_fit_gauss[i, j, *, 2] = erg.dfit
          cii_fit_gauss[i, j, *, 3] = erg.qfit
          cii_fit_gauss[i, j, *, 4] = erg.pfit
        endif
     endfor
   ;print, nx-1 -i, FORMAT='(%"\b%d\b\b\b",$)'
   if ((i mod 5) eq 0)and(i gt 0)and(do_quiet eq 0) then begin
      wset, 17
      tvsclm, cii_s[0:i,*,0], zm=zfac
      tvsclm, cii_d[0:i,*,0], xp=nx/zfac,   yp=0, zm=zfac
      tvsclm, cii_d[0:i,*,1], xp=2*nx/zfac, yp=0, zm=zfac
      tvsclm, cii_d[0:i,*,2], xp=3*nx/zfac, yp=0, zm=zfac
      tvsclm, cii_d[0:i,*,3], xp=4*nx/zfac, yp=0, zm=zfac
      wset, 11
      tvsclm, cii_p[0:i,*,1], zm=zfac                       
      tvsclm, cii_p[0:i,*,2], xp=nx/zfac, yp=0, zm=zfac 
      tvsclm, cii_p[0:i,*,11],xp=2*nx/zfac, yp=0, zm=zfac ; line position Ni II
      tvsclm, cii_d[0:i,*,5], xp=3*nx/zfac, yp=0, zm=zfac ; chi-square, single Gaussian for each C II line
      tvsclm, cii_q[0:i,*,10], xp=4*nx/zfac, yp=0, zm=zfac ; chi-square, double Gaussian for each C II line
      tvsclm, cii_p[0:i,*,13],xp=5*nx/zfac, yp=0, zm=zfac ; chi-square, penta Gaussian
      ;stop
   endif
   if ((i mod 10) eq 0)and(i gt 0) then print, i, nx-1 -i, round((systime(/seconds)- t_start)/(float(i+1))* float(nx-1 -i) / 60.),' min'
endfor

restore, outfile_fuv_temp

cii_cont = reform(cii_s[*,*,0])
;----------------------------------------------------
;-- calculation of non-thermal width
;----------------------------------------------------
;- inputs: (1/e) width of the line profile in km/s

;instr_fwhm = 0.0318 ; instrumental FWHM in A for long FUV  (1380-1406)
instr_fwhm = 0.028556 ; instrumental FWHM in A for short FUV (1331-1358)
;
; Note: the code assumes a fixed temperature so for many 
;       ions, it returns just the observed width.
;----------------------------------------------------
w = reform(cii_q[*,*,7])  > 2. < 60. ; 1/e width in pixel, first component
w = (w/2.) * (dispersion_fuv * 1.0d-3 /133.453) * speed_of_light ; 1/e width in km/s
wnt = iris_nonthermalwidth('C','II', 1334.53, w, instr_fwhm, ti_max=ti_max,Wt_v=wt) ; non-thermal velocity in km/s
wnt_c_ii1 = wnt > 0.1
wt_c_ii1 = wt
ti_c_ii1 = ti_max
;-----------------
w = reform(cii_q[*,*,2])  > 2. < 60. ; 1/e width in pixel, second component
w = (w/2.) * (dispersion_fuv * 1.0d-3 /133.570) * speed_of_light ; 1/e width in km/s
wnt = iris_nonthermalwidth('C','II', 1335.70, w, instr_fwhm, ti_max=ti_max,Wt_v=wt) ; non-thermal velocity in km/s
wnt_c_ii2 = wnt > 0.1
wt_c_ii2 = wt
ti_c_ii2 = ti_max
;----------------------------------

;------------------------------------------------------------------
;- for single Gaussian fit, create a velosity in physical unit
;------------------------------------------------------------------
cii1_vel = reform(cii_d[*,*,3])
q1 = where(cii1_vel eq 0., cz) 
if (at_limb eq 0) then for j=0, ny-1 do cii1_vel[*,j] = cii1_vel[*,j] - euv_temporal_gradient
cii1_vel = reform(cii1_vel) < (master_c1 + del_c1) > (master_c1 - del_c1)
cii1_vel -= master_c1
cii1_vel = cii1_vel * dispersion_fuv * speed_of_light / 133453.2d + 3.0
if (cz gt 0) then cii1_vel(q1) = 0.
;------------------------------------------------------------------
;------------------------------------------------------------------
cii2_vel = reform(cii_d[*,*,0])
q1 = where(cii2_vel eq 0., cz) 
if (at_limb eq 0) then for j=0, ny-1 do cii2_vel[*,j] = cii2_vel[*,j] - euv_temporal_gradient
cii2_vel = reform(cii2_vel) < (master_c2 + del_c2) > (master_c2 - del_c2)
cii2_vel -= master_c2
cii2_vel = cii2_vel * dispersion_fuv * speed_of_light / 133570.8d + 3.0
if (cz gt 0) then cii2_vel(q1) = 0.
;------------------------------------------------------------------

if (do_quiet eq 0) then begin
  window, 2, xs=nx*6/zfac,ys=ny*3/zfac, title='C II 133.4/133.6 nm'

  tvsclm, cii_cont,                  xp=0*nx/zfac, yp=0, zm=zfac    ; continuum
  tvsclm, alog10(cii_d[*,*,1]>.1),   xp=1*nx/zfac, yp=0, zm=zfac ; log amplitude
  tvsclm, alog10(cii_d[*,*,4]>.1),   xp=2*nx/zfac, yp=0, zm=zfac ; log amplitude
  tvsclm, cii_d[*,*,2]<15. > 1.,     xp=3*nx/zfac, yp=0, zm=zfac    ; line width
  tvsclm, cii1_vel < 15. > (-15.0),  xp=4*nx/zfac, yp=0, cb=99, zm=zfac  ; velocity
  tvsclm, cii2_vel<1.3d4 > (-7.0d3), xp=5*nx/zfac, yp=0, zm=zfac, cb=99; velocity
  tvsclm, wnt_c_ii1,                 xp=6*nx/zfac, yp=0, zm=zfac ; non-thermal width, first component
  tvsclm, wnt_c_ii2,                 xp=7*nx/zfac, yp=0, zm=zfac ; non-thermal width, second component

  tvsclm, cii_s[*,*,4]<100.,  xp=0*nx/zfac, yp=ny/zfac ; chi^2 map
  tvsclm, cii_d[*,*,5]<100.,  xp=1*nx/zfac, yp=ny/zfac ; chi^2 map
  tvsclm, cii_q[*,*,10]<100., xp=2*nx/zfac, yp=ny/zfac ; chi^2 map
  tvsclm, cii_p[*,*,13]<100., xp=3*nx/zfac, yp=ny/zfac ; chi^2 map
  tvsclm, cii_p[*,*,10],      xp=4*nx/zfac, yp=ny/zfac ; Ni II line center
  tvsclm, cii_p[*,*,11],      xp=5*nx/zfac, yp=ny/zfac ; Ni II amplitude
  tvsclm, cii_p[*,*,12],      xp=6*nx/zfac, yp=ny/zfac ; Ni II line width

  window, 14, title='reduced chi-square statistics, C II fits'
  s = 5.
  x= reform(cii_s[*, *, 4])  & q=where(x ne 0.) & plot_histogram, x(q), xrange=[0, s]
  loadct,40,/silent
  x= reform(cii_d[*, *, 5])  & q=where(x ne 0.) & plot_histogram, x(q), xrange=[0, s], oplott=70
  x= reform(cii_q[*, *, 10]) & q=where(x ne 0.) & plot_histogram, x(q), xrange=[0, s], oplott=250
  loadct,0,/silent
  set_plot, 'ps'
  device, filename=filepath+'red_chi_pdf_cii.eps', /encapsulated, /color
  x= reform(cii_q[*, *, 10]) & q=where(x ne 0.) & plot_histogram, x(q), xrange=[0, s]
  loadct,40,/silent
  x= reform(cii_d[*, *, 5])  & q=where(x ne 0.) & plot_histogram, x(q), xrange=[0, s], oplott=80

  x= reform(cii_q[*, *, 10]) & q=where(x ne 0.) & plot_histogram, x(q), xrange=[0, s], oplott=250
  loadct,0,/silent
  set_plot, 'x'
  
endif

x= reform(cii_s[*, *, 4])     &      q=where(x ne 0.)
r = percentiles(x(q), value=[0.05, 0.25, 0.5, 0.75, 0.95])
print, '-----------------------------------------------------'
print, 'Percentiles (single Gaussian)', r
x= reform(cii_d[*, *, 5])     &      q=where(x ne 0.)
r = percentiles(x(q), value=[0.05, 0.25, 0.5, 0.75, 0.95])
print, '-----------------------------------------------------'
print, 'Percentiles (double Gaussian)', r
x= reform(cii_q[*, *, 10])     &      q=where(x ne 0.)
r = percentiles(x(q), value=[0.05, 0.25, 0.5, 0.75, 0.95])
print, '-----------------------------------------------------'
print, 'Percentiles (quad Gaussian)', r

a = 0.
save, filename=outfile_cii, c_ii_kont, cii_s, cii_d, cii_q, cii_p, ca_ilo, ca_ihi, del_c1, del_c2, $
      cii_cont, dispersion_fuv, master_c1, master_c2, cii_fit_gauss, ergm, $
      wt_c_ii1, wt_c_ii2, cii1_vel, cii2_vel, dispersion_fuv, $
      wnt_c_ii1, wt_c_ii1, ti_c_ii1, wnt_c_ii2, wt_c_ii2, ti_c_ii2

print, 'average run-time per slit= ',round((systime(/seconds) - t_start) / float(nx)), ' sec'
print, systime()
endif








;------------------------------------------------------------------
;++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;----------------------- Mg II line -------------------------------
;++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;------------------------------------------------------------------
if (do_mg eq 1) then begin

the_line = 'Mg'
mg_range[1] = mg_range[1] < (n_elements(avprof)-1)
if (num_files eq 1) then vv = avprof[Mg_range[0]:Mg_range[1]] else vv = avprof
wing_master =  mean(vv[60:64])

;---------------------------
;- find the line positions
;---------------------------
h = check_if_file_exists(inpath+'mg_lines_position.sav')
if (string2num(h[1]) eq 0.) then guess_mg_pos, inpath, vv, dispersion_nuv, /do_h_line
restore, inpath+'mg_lines_position.sav'
master = [master_k1v1, master_k1v2, master_k1r1, master_k1r2]

if (num_files eq 1) then vv = avprof[Mg_range[0]:Mg_range[1]] else vv = avprof
wing = mean(vv[60:64])

ergm =  analyze_mg(vv, wing, ca_ilo, ca_ihi, master, dispersion_nuv, 20, /plt, do_gauss=1) ;
k_fit_gauss = fltarr(nx, ny, n_elements(ergm.sfit), 4)

if (do_h eq 1) then begin
    master = [master_h1v1, master_h1v2, master_h1r1, master_h1r2]
    ergh =  analyze_mg(vv, wing, ca_ilo, ca_ihi, master, dispersion_nuv, 20, /plt, /do_H_line, do_gauss=1) ;
    h_fit_gauss = fltarr(nx, ny, n_elements(ergh.sfit),  4)
endif

if (do_quiet ne 1) then begin
  window, 13, xs=nx*9/zfac, ys=ny*2/zfac, title='Mg II k line parameters'
  if (do_gauss eq 1) then window, 19, xs=nx*7/zfac, ys=ny/zfac, title='Mg II k Gaussian fits'
endif

t_start = systime(/seconds)
for i=0, nx-1 do begin
   v = readfits(data, hed, nslice=i, /silent) 
   q = where(~finite(v), count)
   if (count gt 0) then v(q) = -10.
   v = v[Mg_range[0]:Mg_range[1], *]
   v = despike_iris_raster(v, master_k3, master_h3)
   
   for j=0, ny-1 do begin
        ;-----------------------------------------------------------------------
        tmp = reform(v[*, j + z0]) & tmp = gauss_smooth(tmp, 1., /edge_truncate)

        ; a symmetric wing range at ~ 1.2 A
        d_stp = 50.0 / (dispersion_nuv / 2.544)
        t_stp = 46.0 / (dispersion_nuv / 2.544)
        wing = 0.5 * (mean(tmp[(master_k3 - d_stp):(master_k3 - t_stp)]) + mean(tmp[(master_k3 + t_stp):(master_k3 + d_stp)]))
        
        if (wing lt 1.)or(~finite(wing)) then wing = wing_master > 1.
        master = [master_k1v1, master_k1v2, master_k1r1, master_k1r2]

        k_check = tmp[master_k1v1:master_k1r2]
        e = check_input_profile(wing, master_k3, master_k1v1, master_k1r2, tmp)
        if (e eq 0) then begin
           erg =  analyze_mg(tmp, wing, ca_ilo, ca_ihi, master, dispersion_nuv, plt=do_plot, do_gauss=do_gauss)

              fe_par[i, j, *]   = erg.velpos
              fe_int[i, j, *]   = erg.fe_int
              mg_wing[i, j]     = wing
              type_k[i, j]      = erg.type
              xmins_k[i, j, *]  = reform(erg.xmins)
              xmaxs_k[i, j, *]  = reform(erg.xmaxs)
              kpar[i, j, *]     = erg.hpar
              em_lines[i,j,*]   = erg.eml
              em_vel[i,j,*]     = erg.ems
              min1v_k[i,j,*]    = erg.xmh1[0,*]
              min1r_k[i,j,*]    = erg.xmh1[1,*]
              mg_bnd[i, j, *]   = erg.band
             
              if (do_gauss eq 1) then begin
                 k_spar[i, j, *]   = erg.spar
                 k_dpar[i, j, *]   = erg.dpar
                 k_tpar[i, j, *]   = erg.tpar

                 k_fit_gauss[i, j, *, 0] = erg.sprf
                 k_fit_gauss[i, j, *, 1] = erg.sfit
                 k_fit_gauss[i, j, *, 2] = erg.dfit
                 k_fit_gauss[i, j, *, 3] = erg.tfit
              endif
              
              ;----------------------------------------------------------------------
              if (do_h eq 1) then begin 
                h_check = tmp[master_h1v1:master_h1r2]
                if (median(h_check) gt 2.)and(max(h_check) gt 2.)and(e eq 0) then begin
                tmp = reform(v[*, j + z0]) & tmp = gauss_smooth(tmp, 1., /edge_truncate)
                master = [master_h1v1, master_h1v2, master_h1r1, master_h1r2]
                erg =  analyze_mg(tmp, wing, ca_ilo, ca_ihi, master, dispersion_nuv, plt=do_plot, /do_H_line, do_gauss=do_gauss)
                fe_par[i, j, *] = erg.velpos
                fe_int[i, j, *] = erg.fe_int

                type_h[i, j]      = erg.type
                xmins_h[i, j, *]  = reform(erg.xmins)
                xmaxs_h[i, j, *]  = reform(erg.xmaxs)
                hpar[i, j, *]   = erg.hpar
                min1v_h[i,j,*] =  erg.xmh1[0,*]
                min1r_h[i,j,*] =  erg.xmh1[1,*]
                mg_bnd[i, j, *]   = erg.band
               
                if (do_gauss eq 1) then begin
                  h_spar[i, j, *]   = erg.spar
                  h_dpar[i, j, *]   = erg.dpar
                  h_tpar[i, j, *]   = erg.tpar

                  h_fit_gauss[i, j, *, 0] = erg.sprf
                  h_fit_gauss[i, j, *, 1] = erg.sfit
                  h_fit_gauss[i, j, *, 2] = erg.dfit
                  h_fit_gauss[i, j, *, 3] = erg.tfit
               endif
           endif
           endif
        ;--------------------------------------------------------------------
        endif
endfor
;print, nx-1 -i, format='(%"\b%d\b\b\b",$)'
  if ((i mod 5) eq 0)and(i gt 0)and(do_quiet eq 0) then begin
     if (do_gauss eq 1) then begin
        wset, 19
        for kapa=0,4 do tvsclm, k_spar[0:i,*,kapa],zm=zfac, xp=kapa*nx/zfac
        tvsclm, k_dpar[0:i,*,7] < 100., zm=zfac, xp=5*nx/zfac
        tvsclm, k_tpar[0:i,*,9] < 20., zm=zfac, xp=6*nx/zfac
     endif
     wset, 13
     for kapa=0,3 do tvsclm, kpar[0:i,*,kapa], zm=zfac, xp=kapa*nx/zfac
     tvsclm, xmaxs_k[0:i,*,0], zm=zfac, xp=4*nx/zfac, cb=99
     tvsclm, alog(reform(xmaxs_k[0:i,*,2])>.1), zm=zfac, xp=5*nx/zfac
     tvsclm, alog(reform(xmaxs_k[0:i,*,3])>.1), zm=zfac, xp=6*nx/zfac
     tvsclm, type_k, zm=zfac, xp=7*nx/zfac
     tvsclm, em_lines[0:i,*,7], zm=zfac, xp=8*nx/zfac
     
     for kapa=6,9 do tvsclm, mg_bnd[0:i,*,kapa]*mg_wing[0:i,*], zm=zfac, xp=(kapa-6)*nx/zfac, yp=ny/zfac
     tvsclm, (mg_bnd[0:i,*,7]/mg_bnd[0:i,*,8]) < 2. > .2, zm=zfac, xp=4*nx/zfac, yp=ny/zfac, cb=99
     tvsclm, mg_bnd[0:i,*,15], zm=zfac, xp=5*nx/zfac, yp=ny/zfac, cb=99
     tvsclm, mg_bnd[0:i,*,17], zm=zfac, xp=6*nx/zfac, yp=ny/zfac, cb=99
     tvsclm, em_lines[0:i,*,1], zm=zfac, xp=7*nx/zfac, yp=ny/zfac
     tvsclm, em_lines[0:i,*,3], zm=zfac, xp=8*nx/zfac, yp=ny/zfac
  endif
  if ((i mod 10) eq 0)and(i gt 0) then print, i, nx-1 -i, round( (systime(/seconds)- t_start)/float(i+1)* float(nx-1 -i) / 60.), ' min'

;  print, i, nx-1 -i, round((systime(/seconds)- t_start)/float(i+1)* float(nx-1 -i) / 60.),'     min'
;  if ((i mod 5) eq 0) then begin
;     print, '-------------------------------------------------------------------------'
;     print, good_mean(k_dpar[0:i, *,6]), stddev(k_dpar[0:i, *,6]), max(k_dpar[0:i, *,6])
;     print, good_mean(k_tpar[0:i, *,9]), stddev(k_tpar[0:i, *,9]), max(k_tpar[0:i, *,9])
;     print, '-------------------------------------------------------------------------'
;  endif
 
endfor
save, filename=filepath+'em_lines.sav', em_lines
;-------------------------------------------------
;-- systematic velocity residuals
;-------------------------------------------------
restore, outfile_nuv_temp
temporal_gradient_mg = temporal_gradient
;--------------------------------
;-- remove wavelength slope
;--------------------------------
if (at_limb eq 0) then begin
   for kapa=0, 1 do for j=0, ny-1 do fe_par[*,j,kapa] = fe_par[*,j,kapa] - temporal_gradient_mg
endif

for eta=0,3 do begin
   vm = reform(em_vel[*,*,eta])
   if (at_limb eq 0) then for j=0, ny-1 do vm[*,j] -= temporal_gradient_mg
   mean_vel = good_mean(vm)
   em_vel[*,*,eta] = (vm - mean_vel) * dispersion_nuv * speed_of_light / 279500.0d   ; velocity in km/s
endfor   

if (at_limb eq 0) then begin
   for kapa=0, 3 do for j=0, ny-1 do xmins_k[*,j,kapa*2] = xmins_k[*,j,kapa*2] - temporal_gradient_mg
   for kapa=0, 2 do for j=0, ny-1 do xmaxs_k[*,j,kapa*2] = xmaxs_k[*,j,kapa*2] - temporal_gradient_mg
   for j=0, ny-1 do mg_bnd[*,j, 17] = mg_bnd[*,j, 17] - temporal_gradient_mg
   for j=0, ny-1 do mg_bnd[*,j, 15] = mg_bnd[*,j, 15] - temporal_gradient_mg
   for j=0, ny-1 do k_spar[*,j,1] = k_spar[*,j,1] - temporal_gradient_mg
endfor
;----------------------------------------
;-- normalize emission peak velocities 
;----------------------------------------
k1_width = (reform(min1r_k[*,*,0]) - reform(min1v_k[*,*,0])) > 20. < 80.

x = reform(xmaxs_k[*,*, 0])
y = reform(xmins_k[*,*, 1]) 
q = where(type_k eq 5, count)
if (count ge 1) then y(q) = x(q)
y = y -  good_mean(y) ;master_k3
mgk_h3_vel = y * dispersion_nuv * speed_of_light / 279500.0d
if (abs(good_mean(mgk_h3_vel)) gt 1.) then mgk_h3_vel = mgk_h3_vel - good_mean(mgk_h3_vel(where(type_k eq 1)))
;---------------------------------------------------------------
;---------------------------------------------------------------
x = reform(xmins_k[*,*, 0])  ; first minimum > K1v position
x = x - good_mean(x) ;master_k1v
x = x < del_mg > (- del_mg)
mgk_h1v_vel = x * dispersion_nuv * speed_of_light / 279500.0d
if (abs(good_mean(mgk_h1v_vel)) gt 1.) then mgk_h1v_vel = mgk_h1v_vel - good_mean(mgk_h1v_vel(where(type_k eq 1)))
;---------------------------------------------------------------
x = reform(xmins_k[*,*, 2])  ; last minimum > K1r position
y = reform(xmins_k[*,*, 1])  ; in case of umbral profile, this is the last minimum
q = where(type_k eq 5, count)
if (count ge 1) then x(q) = y(q)
x = x -  good_mean(x) ;master_k1r
x = x < del_mg > (- del_mg)
mgk_h1r_vel = x * dispersion_nuv * speed_of_light / 279500.0d
if (abs(good_mean(mgk_h1r_vel)) gt 1.) then mgk_h1r_vel = mgk_h1r_vel - good_mean(mgk_h1r_vel(where(type_k eq 1)))
;---------------------------------------------------------------
x = reform(xmaxs_k[*,*, 1])  ; K2r
q=where(type_k eq 1) & x = x -  good_mean(x) ;master_k2r  ; if normal profile
xm = good_mean(x(q))
q=where(type_k eq 5, count)   &   if (count ge 1) then  x(q) = xm
;x = x < (del_mg+4.) > (- del_mg )
mgk_h2r_vel = x * dispersion_nuv * speed_of_light / 279500.0d
if (abs(good_mean(mgk_h2r_vel)) gt 1.) then mgk_h2r_vel = mgk_h2r_vel - good_mean(mgk_h2r_vel(where(type_k eq 1)))
;---------------------------------------------------------------
x = reform(xmaxs_k[*,*, 0])  ; K2v
q=where(type_k eq 1) & x = x -  good_mean(x) ;master_k2v  ; if normal profile
xm = good_mean(x(q))
q=where(type_k eq 5, count)   &    if (count ge 1) then x(q) = xm
;x = x < (del_mg) > (- del_mg -4.)
mgk_h2v_vel = x * dispersion_nuv * speed_of_light / 279500.0d
if (abs(good_mean(mgk_h2v_vel)) gt 1.) then mgk_h2v_vel = mgk_h2v_vel - good_mean(mgk_h2v_vel(where(type_k eq 1)))
;---------------------------------------------------------------
;***************************************************************
;----------------------------------------
;-- normalize COG velocities 
;----------------------------------------
cogk_1 = reform(mg_bnd[*,*,15]) < (master_k3 + del_mg) > (master_k3 - del_mg)  &   cogk_1 -= master_k3
cogk_5 = reform(mg_bnd[*,*,17]) < (master_k3 + del_mg) > (master_k3 - del_mg)  &   cogk_5 -= master_k3

cogk_1 *= dispersion_nuv * speed_of_light / 279500.0d   & cogk_1 -= good_mean(cogk_1)     ;&    cogk_1 *= -1.0
cogk_5 *= dispersion_nuv * speed_of_light / 279500.0d   & cogk_5 -= good_mean(cogk_5)     ;&    cogk_5 *= -1.0

gauss_vel_k = reform(k_spar[*,*,1])
mean_v = good_mean(gauss_vel_k)  &    gauss_vel_k  -= mean_v
;gauss_vel_k -= mean(median(gauss_vel_k, 7))
gauss_vel_k *= dispersion_nuv * speed_of_light / 279500.0d
;----------------------------------------------------
dpar_v_k = reform(k_dpar[*,*,1])
mean_v = good_mean(dpar_v_k)  &    dpar_v_k  -= mean_v
dpar_v_k *= dispersion_nuv * speed_of_light / 279500.0d

dpar_r_k = reform(k_dpar[*,*,4])
mean_v = good_mean(dpar_r_k)  &    dpar_r_k  -= mean_v
dpar_r_k *= dispersion_nuv * speed_of_light / 279500.0d
;----------------------------------------------------
tpar_v_k = reform(k_tpar[*,*,4])
mean_v = good_mean(tpar_v_k)  &    tpar_v_k  -= mean_v
tpar_v_k *= dispersion_nuv * speed_of_light / 279500.0d

tpar_r_k = reform(k_tpar[*,*,7])
mean_v = good_mean(tpar_r_k)  &    tpar_r_k  -= mean_v
tpar_r_k *= dispersion_nuv * speed_of_light / 279500.0d
;----------------------------------------------------

;-------------------------------------------------
;-- normalize Fe velocities to average QS
;-------------------------------------------------
for i=0, 1 do begin
  s = percentiles(cont, value=[0.5, .8, .85])
  q=where((cont gt (s[0]<0.5))and(cont lt (s[2] < 1.5)))

  x= reform(fe_par[*,*,i])
  xm = good_mean(x(q))
  case i of 
    0: x = (x- xm)*dispersion_nuv * 3.0d+8 / 279768.0d
    1: x = (x- xm)*dispersion_nuv * 3.0d+8 / 279886.0d
  endcase
  fe_par[*, *,i] = x - 200. ; 200 m/s convective Blueshift, should be perhaps larger.
endfor

;------------------------------------------
k2_width = reform(kpar[*,*,1]) > 2.
k_wb_width = reform(kpar[*,*,2]) > 10.
k_wb_int = reform(kpar[*,*,3])
kvr = reform(kpar[*,*,4])  
sw = where(~finite(kvr), count) &  if (count gt 1) then kvr(sw) = 1.0 

k1v_b = reform(mg_bnd[*,*,13])
k1r_b = reform(mg_bnd[*,*,0])
k2v_b = reform(mg_bnd[*,*,7])
k2r_b = reform(mg_bnd[*,*,8])
k3_b = reform(mg_bnd[*,*,6])
kindex1 =  reform(mg_bnd[*,*,9])
kindex5 =  reform(mg_bnd[*,*,11])
;------------------------------------------
;------------------------------------------
k2v = cogk_1 - cogk_1  & t1 =  reform(xmaxs_k[*,*,1])
k2r = cogk_1 - cogk_1  & t2 =  reform(xmaxs_k[*,*,3])
k3 = cogk_1 - cogk_1   & t3 =  reform(xmins_k[*,*,3])
k1 = cogk_1 - cogk_1   & t4 =  reform(xmins_k[*,*,5])

q = where(type_k eq 1, count) &  if (count gt 0) then k2v(q) = t1(q)
q = where(type_k eq 1, count) &  if (count gt 0) then k2r(q) = t2(q)
q = where(type_k eq 1, count) &  if (count gt 0) then k3(q) = t3(q)
q = where(type_k eq 1, count) &  if (count gt 0) then k1(q) = t4(q)
q = where(type_k eq 5, count) &  if (count gt 0) then k1(q) = t3(q)
q=where(type_k eq 1) & qq=where(type_k ne 1)
k2v(qq) =min(k2v(q))
k2r(qq) =min(k2r(q))
k3(qq) =min(k3(q))

mgk_core = cogk_1 - cogk_1
m1 = reform(xmaxs_k[*,*,1])
m2 = reform(xmaxs_k[*,*,3])
q = where(type_k eq 1, count) &  if (count gt 0) then mgk_core(q) = m2(q)
q = where(type_k eq 5, count) &  if (count gt 0) then mgk_core(q) = m1(q)
q = where(type_k eq 3, count) &  if (count gt 0) then mgk_core(q) = m1(q)
q = where(type_k eq 2, count) &  if (count gt 0) then mgk_core(q) = m2(q)

;------------------------------------------
;------------------------------------------
if (do_h eq 1) then begin
  ;--------------------------------
  ;-- remove wavelength slope
  ;--------------------------------
  for kapa=0, 3 do for j=0, ny-1 do xmins_h[*,j,kapa*2] = xmins_h[*,j,kapa*2] - temporal_gradient_mg
  for kapa=0, 2 do for j=0, ny-1 do xmaxs_h[*,j,kapa*2] = xmaxs_h[*,j,kapa*2] - temporal_gradient_mg
  for j=0, ny-1 do mg_bnd[*,j, 26] = mg_bnd[*,j, 26] - temporal_gradient_mg
  for j=0, ny-1 do mg_bnd[*,j, 28] = mg_bnd[*,j, 28] - temporal_gradient_mg
  for j=0, ny-1 do h_spar[*,j,0] = h_spar[*,j,0] - temporal_gradient_mg

  ;----------------------------------------
  ;-- normalize emission peak velocities 
  ;----------------------------------------
  h1_width = (reform(min1r_h[*,*,0]) - reform(min1v_h[*,*,0])) > 5. < 80.
  ;---------------------------------------------------------------
  x = reform(xmaxs_h[*,*, 0])
  y = reform(xmins_h[*,*, 1]) 
  q = where(type_h eq 5, count)
  if (count ge 1) then y(q) = x(q)
  y = y -  good_mean(y) ;master_h3
  mgh_h3_vel = y * dispersion_nuv * speed_of_light / 280200.0d
  if (abs(good_mean(mgh_h3_vel)) gt 1.) then mgh_h3_vel = mgh_h3_vel - good_mean(mgh_h3_vel(where(type_h eq 1)))
  ;---------------------------------------------------------------
  ;---------------------------------------------------------------
  x = reform(xmins_h[*,*, 0])  ; first minimum > H1v position
  x = x -  good_mean(x) ;master_h1v
  x = x < del_mg > (- del_mg)
  mgh_h1v_vel = x * dispersion_nuv * speed_of_light / 280200.0d
  if (abs(good_mean(mgh_h1v_vel)) gt 1.) then mgh_h1v_vel = mgh_h1v_vel - good_mean(mgh_h1v_vel(where(type_h eq 1)))
  ;---------------------------------------------------------------
  x = reform(xmins_h[*,*, 2])  ; last minimum > H1r position
  y = reform(xmins_h[*,*, 1])  ; in case of umbral profile, this is the last minimum
  q = where(type_h eq 5, count)
  if (count ge 1) then x(q) = y(q)
  x = x -  good_mean(x) ;master_h1r
  x = x < del_mg > (- del_mg)
  mgh_h1r_vel = x * dispersion_nuv * speed_of_light / 280200.0d
  if (abs(good_mean(mgh_h1r_vel)) gt 1.) then mgh_h1r_vel = mgh_h1r_vel - good_mean(mgh_h1r_vel(where(type_h eq 1)))
  ;---------------------------------------------------------------
  x = reform(xmaxs_h[*,*, 1])  ; H2r
  q=where(type_h eq 1, count) & x = x -  good_mean(x) ;master_h2r  ; if normal profile
  if (count ge 1) then xm = good_mean(x(q))  else xm = 0.
  q=where(type_h eq 5, count)   &    if (count ge 1) then x(q) = xm
  x = x < (del_mg+4.) > (- del_mg)
  mgh_h2r_vel = x * dispersion_nuv * speed_of_light / 280200.0d
  if (abs(good_mean(mgh_h2r_vel)) gt 1.) then mgh_h2r_vel = mgh_h2r_vel - good_mean(mgh_h2r_vel(where(type_h eq 1)))
  ;---------------------------------------------------------------
  x = reform(xmaxs_h[*,*, 0])  ; H2v
  q=where(type_h eq 1, count) & x = x -  good_mean(x) ;master_h2v  ; if normal profile
  if (count ge 1) then xm = good_mean(x(q)) else xm=0.
  q=where(type_h eq 5, count)   &  if (count ge 1) then   x(q) = xm
  x = x < (del_mg) > (- del_mg -4.)
  mgh_h2v_vel = x * dispersion_nuv * speed_of_light / 280200.0d
  if (abs(good_mean(mgh_h2v_vel)) gt 1.) then mgh_h2v_vel = mgh_h2v_vel - good_mean(mgh_h2v_vel(where(type_h eq 1)))
  ;----------------------------------------
  ;-- normalize COG velocities 
  ;----------------------------------------
  cogh_1 = reform(mg_bnd[*,*,26]) < (master_h3 + del_mg) > (master_h3 - del_mg)  &   cogh_1 -= master_h3
  cogh_5 = reform(mg_bnd[*,*,28]) < (master_h3 + del_mg) > (master_h3 - del_mg)  &   cogh_5 -= master_h3

  cogh_1 *= dispersion_nuv * speed_of_light / 280200.0d   & cogh_1 -= good_mean(cogh_1) 
  cogh_5 *= dispersion_nuv * speed_of_light / 280200.0d   & cogh_5 -= good_mean(cogh_5)

  gauss_vel_h = reform(h_spar[*,*,0])
  gauss_vel_h -= good_mean(gauss_vel_h)
  gauss_vel_h *= dispersion_nuv * speed_of_light / 280200.0d

  h2_width = reform(hpar[*,*,1]) > 2.
  h_wb_width = reform(hpar[*,*,2]) > 10.
  h_wb_int = reform(hpar[*,*,3])
  hvr = reform(hpar[*,*,4])     ;& hvr(where(~finite(hvr))) = 1.0 
  sw = where(~finite(hvr), count) &  if (count gt 1) then hvr(sw) = 1.0 

  h1v_b = reform(mg_bnd[*,*,25])
  h1r_b = reform(mg_bnd[*,*,25])
  h2v_b = reform(mg_bnd[*,*,20])
  h2r_b = reform(mg_bnd[*,*,21])
  h3_b  = reform(mg_bnd[*,*,19])
  hindex1 =  reform(mg_bnd[*,*,22])
  hindex5 =  reform(mg_bnd[*,*,24])

  h2v = cogh_1 - cogh_1  & t1 =  reform(xmaxs_h[*,*,1])
  h2r = cogh_1 - cogh_1  & t2 =  reform(xmaxs_h[*,*,3])
  h3 = cogh_1 - cogh_1   & t3 =  reform(xmins_h[*,*,3])
  h1 = cogh_1 - cogh_1   & t4 =  reform(xmins_h[*,*,5])

  q = where(type_h eq 1, count) &  if (count gt 0) then h2v(q) = t1(q)
  q = where(type_h eq 1, count) &  if (count gt 0) then h2r(q) = t2(q)
  q = where(type_h eq 1, count) &  if (count gt 0) then h3(q) = t3(q)
  q = where(type_h eq 1, count) &  if (count gt 0) then h1(q) = t4(q)
  q = where(type_h eq 5, count) &  if (count gt 0) then h1(q) = t3(q)

  q=where(type_h eq 1, count_q) & qq=where(type_h ne 1, count_qq)
  if (count_q gt 0)and(count_qq gt 0) then begin
     h2v(qq) =min(h2v(q))
     h2r(qq) =min(h2r(q))
     h3(qq) =min(h3(q))
  endif

  mgh_core = cogh_1 - cogh_1
  m1 = reform(xmaxs_h[*,*,1])
  m2 = reform(xmaxs_h[*,*,3])
  q = where(type_h eq 1, count) &  if (count gt 0) then mgh_core(q) = m2(q)
  q = where(type_h eq 5, count) &  if (count gt 0) then mgh_core(q) = m1(q)
  q = where(type_h eq 3, count) &  if (count gt 0) then mgh_core(q) = m1(q)
  q = where(type_h eq 2, count) &  if (count gt 0) then mgh_core(q) = m2(q)

endif
;------------------------------------------
;------------------------------------------
if (do_quiet eq 0) then begin
  window, 0,xs=nx*5/zfac,ys=ny/zfac,title='Mg II k velocities: belnds 1, 2, H2v, H2r, H23, COG1, COG5'
  tvsclm, fe_par[*,*,0] < 2000. > (-2000.), xp=0*nx/zfac, zm=zfac, cb=99
  tvsclm, mgk_h2v_vel, ilim=[-10.0, 10.0], xp=1*nx/zfac, zm=zfac, cb=99
  tvsclm, mgk_h2r_vel, ilim=[-10.0, 10.0], xp=2*nx/zfac, zm=zfac, cb=99
  tvsclm, mgk_h3_vel, ilim=[-15.0, 15.0], xp=3*nx/zfac, zm=zfac, cb=99
  tvsclm, cogk_1,     ilim=[-7., 7.], xp=4*nx/zfac, zm=zfac, cb=99
;------------------------------------------
  window, 7,xs=nx*6/zfac,ys=ny*2/zfac,title='Mg II intensities: K3, K2v, K2 width, K-index, K3/V, V/R, V/R, WB emission widths, WB-int, K3, K1, Mg II k core'
  tvsclm, k3_b * mg_wing, zm=zfac
  tvsclm, k2v_b * mg_wing,     xp=1*nx/zfac, zm=zfac
  tvsclm, k2_width,     xp=2*nx/zfac, zm=zfac
  tvsclm, kindex5 * mg_wing,   xp=3*nx/zfac, zm=zfac
  tvsclm, mg_bnd[*,*,6]/(mg_bnd[*,*,7]>0.1), xp=4*nx/zfac, zm=zfac
  tvsclm, mg_bnd[*,*,7]/(mg_bnd[*,*,8]>0.1), xp=5*nx/zfac, zm=zfac, cb=99

  tvsclm, kvr < 3. > .3,  xp=0*nx/zfac, yp=ny/zfac, zm=zfac, cb=99
  tvsclm, k_wb_width,     xp=1*nx/zfac, yp=ny/zfac, zm=zfac
  tvsclm, k_wb_int,       xp=2*nx/zfac, yp=ny/zfac, zm=zfac
  tvsclm, k3 < 3.,        xp=3*nx/zfac, yp=ny/zfac, zm=zfac
  tvsclm, k1,             xp=4*nx/zfac, yp=ny/zfac, zm=zfac
  tvsclm, alog10(mgk_core > 1.), xp=5*nx/zfac, yp=ny/zfac, zm=zfac
;------------------------------------------
;------------------------------------------
  if (do_h eq 1) then begin
     window, 1,xs=nx*5/zfac,ys=ny/zfac,title='Mg II h velocities: H2v, H2r, H23, COG5, COG1'
     tvsclm, mgh_h2v_vel, ilim=[-10.0, 10.0], xp=0*nx/zfac, zm=zfac, cb=99     
     tvsclm, mgh_h2r_vel, ilim=[-10.0, 10.0], xp=1*nx/zfac, zm=zfac, cb=99
     tvsclm, mgh_h3_vel,  ilim=[-15.0, 15.0],  xp=2*nx/zfac, zm=zfac, cb=99
     tvsclm, cogh_5,  ilim=[-6., 6.],       xp=3*nx/zfac, zm=zfac, cb=99
     tvsclm, cogh_1,  ilim=[-6., 6.],       xp=4*nx/zfac, zm=zfac, cb=99

 window, 8,xs=nx*6/zfac,ys=ny*2/zfac,title='Mg II intensities: H3, H2v, H2 width, H-index, H3/V, V/R, V/R, WB emission widths, WB-int, H3, H1, Mg II h core'
     tvsclm, h3_b * mg_wing, zm=zfac
     tvsclm, h2v_b * mg_wing,    xp=1*nx/zfac, zm=zfac
     tvsclm, h2_width,    xp=2*nx/zfac, zm=zfac
     tvsclm, hindex5 * mg_wing,  xp=3*nx/zfac, zm=zfac
     tvsclm, mg_bnd[*,*,20]/(mg_bnd[*,*,21]>0.1), xp=4*nx/zfac, zm=zfac
     tvsclm, mg_bnd[*,*,20]/(mg_bnd[*,*,19]>0.1), xp=5*nx/zfac, zm=zfac

     tvsclm, hvr < 3. > .3,  xp=0*nx/zfac, yp=ny/zfac, zm=zfac, cb=99
     tvsclm, h_wb_width,     xp=1*nx/zfac, yp=ny/zfac, zm=zfac
     tvsclm, h_wb_int,       xp=2*nx/zfac, yp=ny/zfac, zm=zfac
     tvsclm, h3 < 3.,        xp=3*nx/zfac, yp=ny/zfac, zm=zfac
     tvsclm, h1,             xp=4*nx/zfac, yp=ny/zfac, zm=zfac
     tvsclm, alog10(mgh_core > 1.), xp=5*nx/zfac, yp=ny/zfac, zm=zfac
     ;------------------------------------------
     print, 'profile type for the h line (0: reversal free, 1: normal, 5: umbral)'
     for i=0,5 do begin & q=where(type_h eq i, count) &  print, i, count & endfor
     ;------------------------------------------
endif
  if (do_gauss ne 0) then begin
    window, 14, title='reduced chi-square statistics, Mg II fits'
    s = 9.
    x= reform(k_spar[*, *, 4]) & q=where(x ne 0.) & plot_histogram, x(q), xrange=[0, s]
    loadct, 40,/silent
    x= reform(k_dpar[*, *, 7]) & q=where(x ne 0.) & plot_histogram, x(q), xrange=[0, s], oplott=70
    x= reform(k_tpar[*, *, 9]) & q=where(x ne 0.) & plot_histogram, x(q), xrange=[0, s], oplott=250
    loadct, 0,/silent
  endif  
endif ;else begin
   if (do_gauss ne 0) then begin
     x= reform(k_spar[*, *, 4])     &      q=where(x ne 0.)
     r = percentiles(x(q), value=[0.05, 0.25, 0.5, 0.75, 0.95])
     print, '-----------------------------------------------------'
     print, 'Percentiles (single Gaussian)', r
     x= reform(k_dpar[*, *, 7])     &      q=where(x ne 0.)
     r = percentiles(x(q), value=[0.05, 0.25, 0.5, 0.75, 0.95])
     print, '-----------------------------------------------------'
     print, 'Percentiles (double Gaussian)', r
     x= reform(k_tpar[*, *, 9])     &      q=where(x ne 0.)
     r = percentiles(x(q), value=[0.05, 0.25, 0.5, 0.75, 0.95])
     print, '-----------------------------------------------------'
     print, 'Percentiles (triple Gaussian)', r
     print, '-----------------------------------------------------'
  endif  
;endelse

print, 'profile type for the k line (0: reversal free, 1: normal, 5: umbral)'
for i=0,5 do begin & q=where(type_k eq i, count) &  print, i, count & endfor

mg_overview = (k3_b * mg_wing)>.1  
r = percentiles(mg_overview, value=[0.05, 0.25, 0.5, 0.75, 0.95])   
write_jpeg, name_term+'k3_log.jpg', bytscl(alog(mg_overview > r[0] < r[4])), quality=100
mg_overview = (k3_b * mg_wing)>.1
r = percentiles(mg_overview, value=[0.05, 0.25, 0.5, 0.75, 0.95])
write_jpeg, name_term+'k3.jpg', bytscl(mg_overview > r[0] < r[4]), quality=100
mg_overview = (k1r_b * mg_wing > 10.)
r = percentiles(mg_overview, value=[0.05, 0.25, 0.5, 0.75, 0.95])
write_jpeg, name_term+'k1r_log.jpg', bytscl(alog(mg_overview > r[0] < r[4])), quality=100
mg_overview = mg_bnd[*,*,7]/(mg_bnd[*,*,6]>0.1) < 8. > .2
r = percentiles(mg_overview, value=[0.05, 0.25, 0.5, 0.75, 0.95])
write_jpeg, name_term+'em_strength.jpg', bytscl(mg_overview > r[0] < r[4]), quality=100

mg_overview = reform(mg_bnd[*,*,15]) & q=where(mg_overview gt 1.)
r = percentiles(mg_overview(q), value=[0.05, 0.25, 0.5, 0.75, 0.95])
write_jpeg, name_term+'k_cog.jpg', bytscl(mg_overview > r[0] < r[4]), quality=100
   
mg_overview = ((mg_bnd[*,*,0]/(mg_bnd[*,*,5])>.1) > .3)
r = percentiles(mg_overview, value=[0.05, 0.25, 0.5, 0.75, 0.95])
write_jpeg, name_term+'k1r_2_wing.jpg', bytscl(alog(mg_overview > r[0] < r[4])), quality=100
if (do_h eq 1) then begin 
    mg_overview = (k3_b / (h3_b > 0.1)) < 1.9 > .9
    r = percentiles(mg_overview, value=[0.05, 0.25, 0.5, 0.75, 0.95])
    write_jpeg, name_term+'k3_2_h3.jpg', bytscl(mg_overview > r[0] < r[4]), quality=100
 endif


if (do_gauss eq 1) then begin
  if (do_h eq 1) then begin;----------- for h+k profiles + /do_gauss
    save, filename=outfile_mg, fe_par, fe_int, xmaxs_k, xmins_k,xmaxs_h, xmins_h, mg_bnd, mg_wing, $
    h1_width, h2_width, h_wb_width, h_wb_int, hvr, h3_b,h2v_b,h2r_b, mgh_h2v_vel, mgh_h2r_vel, mgh_h3_vel, $
    k1_width, k2_width, k_wb_width, k_wb_int, kvr, k3_b,k2v_b,k2r_b, mgk_h2v_vel, mgk_h2r_vel, mgk_h3_vel, $
    mgk_h1v_vel, mgk_h1r_vel, cogk_1, cogk_5, mgk_core, kpar, k2v, k2r, k3, kindex1, kindex5, k_spar, gauss_vel_k, k1v_b, k1r_b, type_k,$
    mgh_h1v_vel, mgh_h1r_vel, cogh_1, cogh_5, mgh_core, hpar, h2v, h2r, h3, hindex1, hindex5, h_spar, gauss_vel_h, h1v_b, h1r_b, type_h,$
    master_h3, master_k3, dispersion_nuv, master_k2v, master_k2r, master_h2v, master_h2r, temporal_gradient, $
    em_lines, em_vel, k_spar, k_dpar, k_tpar, h_spar, h_dpar, h_tpar, tpar_v_k, tpar_r_k, dpar_v_k, dpar_r_k , $
    k_fit_gauss, h_fit_gauss
  endif else begin;----------- for k lines profiles + /do_gauss
    save, filename=outfile_mg, fe_par, fe_int, xmaxs_k, xmins_k, mg_bnd, mg_wing, $
    k1_width, k2_width, k_wb_width, k_wb_int, kvr, k3_b,k2v_b,k2r_b, $
    mgk_h1v_vel, mgk_h1r_vel, mgk_h2v_vel, mgk_h2r_vel, mgk_h3_vel, $
    cogk_1, cogk_5, mgk_core, kpar, k2v, k2r, k3, kindex1, kindex5, k_spar, gauss_vel_k,$
    master_k3, dispersion_nuv, master_k2v, master_k2r,  temporal_gradient, k1v_b, k1r_b, type_k,$
    em_lines, em_vel, k_spar, k_dpar, k_tpar, tpar_v_k, tpar_r_k, dpar_v_k, dpar_r_k, k_fit_gauss
  endelse
  endif else begin
    if (do_h eq 1) then begin;----------- for h+k profiles
    save, filename=outfile_mg, fe_par, fe_int, xmaxs_k, xmins_k,xmaxs_h, xmins_h, mg_bnd, mg_wing, $
    h1_width, h2_width, h_wb_width, h_wb_int, hvr, h3_b,h2v_b,h2r_b, mgh_h2v_vel, mgh_h2r_vel, mgh_h3_vel, $
    k1_width, k2_width, k_wb_width, k_wb_int, kvr, k3_b,k2v_b,k2r_b, mgk_h2v_vel, mgk_h2r_vel, mgk_h3_vel, $
    mgk_h1v_vel, mgk_h1r_vel, cogk_1, cogk_5, mgk_core, kpar, k2v, k2r, k3, kindex1, kindex5, k_spar, gauss_vel_k, k1v_b, k1r_b, type_k,$
    mgh_h1v_vel, mgh_h1r_vel, cogh_1, cogh_5, mgh_core, hpar, h2v, h2r, h3, hindex1, hindex5, h_spar, gauss_vel_h, h1v_b, h1r_b, type_h,$
    master_h3, master_k3, dispersion_nuv, master_k2v, master_k2r, master_h2v, master_h2r, temporal_gradient, $
    em_lines, em_vel
  endif else begin;----------- for k lines profiles
    save, filename=outfile_mg, fe_par, fe_int, xmaxs_k, xmins_k, mg_bnd, mg_wing, $
    k1_width, k2_width, k_wb_width, k_wb_int, kvr, k3_b,k2v_b,k2r_b, $
    mgk_h1v_vel, mgk_h1r_vel, mgk_h2v_vel, mgk_h2r_vel, mgk_h3_vel, $
    cogk_1, cogk_5, mgk_core, kpar, k2v, k2r, k3, kindex1, kindex5, k_spar, gauss_vel_k,$
    master_k3, dispersion_nuv, master_k2v, master_k2r,  temporal_gradient, k1v_b, k1r_b, type_k,$
    em_lines, em_vel
  endelse
endelse   

print, 'Mg: average run-time per slit = ',round((systime(/seconds) - t_start) / float(nx)), ' sec'
print, systime()

endif


endfor




end
