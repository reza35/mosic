pro iris_dispersion, polis_path

;+
;===============================================================  
; procedure : iris_dispersion.pro
;  
; purpose: to double check the dispersion in IRIS data.
; It returns estimated dispersion as well as the dispersion in the file header.
; It is also the place to select spectral ranges interactively.
; The program first makes a quick check of the dispersion and
; then the user interactively selects the spectral range
;  
; 2015/12/27 : created
; 2016/04/13 : updated documentation
; 2016/05/08 : facelift !  
; 2016/05/24 : bug fix for short Mg range
; 2016/10/07 : works when only the Mg II k line was observed.
;
;  
; R.Rezaei @ IAC                         e-mail:  rrezaei@iac.es      
;===============================================================
;-

print, polis_path
files_d = read_dir(polis_path , filter='aver_prof*.sav')
num = n_elements(files_d.files)
if (num gt 1) then begin
   files_d = read_dir(polis_path , filter='aver_prof*CII*.sav')
endif   
restore, polis_path + files_d.files[0]

f_d = read_dir(polis_path , filter='iris_*im.fits')
new_file = strsplit(f_d.files[0], '.,_',/extract)
outfile= polis_path + new_file[0]+'_'+new_file[1]+'_'+new_file[2]+'_'+new_file[3]+'_'+'avprof.sav'

print
print, outfile
print, '------------------'

plot, avprof
print
print, '-----------------------------------------------------------------'
print, 'To estimate the spectral dispersion in FUV FUV data ....'
print
n3 = fxpar(hed, 'NAXIS3') &  n4 = fxpar(hed, 'NAXIS4') 
if (n3 gt n4)and(n3 gt 1000.) then disp= fxpar(hed, 'CDELT3') else disp=fxpar(hed, 'CDELT4')

vv = avprof
print, 'select a region around C II lines (FUV dispersion)'

read_plot_click, 2, vv, vv, pos


vv = avprof[pos[0]:pos[1]]
print, 'Now click on the line core'
read_plot_click, 2, vv, vv, pos

lpff, vv[pos[0]-10:pos[0]+10], lp  &  lp += (pos[0]-10)  & print, pos[0], lp
pos[0] = lp  & vline, lp
lpff, vv[pos[1]-10:pos[1]+10], lp  &  lp += (pos[1]-10)  & print, pos[1], lp
pos[1] = lp  & vline, lp

disp_fuv = (1335.708 - 1334.532)/(float(pos[1] - pos[0]))  &  disp3 = disp_fuv
if (disp_fuv lt 0.014) then disp_fuv = 0.012980
if (disp_fuv gt 0.014)and(disp_fuv lt 0.028) then disp_fuv = 0.02596
if (disp_fuv gt 0.028)and(disp_fuv lt 0.06) then disp_fuv = 0.05192
if (disp_fuv gt 0.06) then begin
   print, 'The FUV dispersion is larger than 6 mA which is unusual !'
   print, 'Please check the average profile.'
   stop
endif   
;----------------------------------------------------------------------
print, '-----------------------------------------------------------------'
print, 'To estimate the spectral dispersion in NUV data ....'
if (num gt 1) then begin
   files_d = read_dir(polis_path , filter='aver_prof*MgII*.sav')
   infile= polis_path + files_d.files[0]
   restore, infile
endif

vv = avprof
plot, vv, /xstyle, psym=-1
s = 'd'
print, 'If both Mg h and k are observed, enter h, otherwise enetr k'
read, s
if (s eq 'h') then begin
  print, 'select a region around both k and h lines'
  read_plot_click, 2, vv, vv, pos
  vv = avprof[pos[0]:pos[1]]
  print, 'click on the two emission cores of k3 and h3'
  read_plot_click, 2, vv, vv, pos

  lpff, vv[pos[0]-2:pos[0]+2], lp  &  lp += (pos[0]-2)  & print, pos[0], lp
  pos[0] = lp  & vline, lp
  lpff, vv[pos[1]-2:pos[1]+2], lp  &  lp += (pos[1]-2)  & print, pos[1], lp
  pos[1] = lp  & vline, lp
  disp_nuv = (2802.704 - 2795.528)/(float(pos[1] - pos[0]))
endif else begin
  print, 'select a region around the Mg k line'
  read_plot_click, 2, vv, vv, pos
  vv = avprof[pos[0]:pos[1]]
  print, 'click on the two emission peaks of k2v and k2r'
  read_plot_click, 2, vv, vv, pos

  lpff, vv[pos[0]-2:pos[0]+1], lp  &  lp += (pos[0]-2)  & print, pos[0], lp
  vline, lp
  
  lpff, vv[pos[1]-1:pos[1]+2], rp  &  rp += (pos[1]-1)  & print, pos[1], rp
  vline, rp
  disp_nuv = (0.3)/(rp - lp)
  if (disp_nuv lt 0.03) then begin
    maximum,  5, vv[pos[0]-5:pos[0]+4], lpx & print, lpx[3]+ (pos[0] - 5.0)
    maximum,  4, vv[pos[1]-4:pos[1]+5], rpx & print, rpx[3]+ (pos[1] - 4.0)
    disp_nuv = (0.3)/(rp - lp)
  endif   
endelse   

disp4 = disp_nuv 
if (disp_nuv lt 0.014) then disp_nuv = 0.012720
if (disp_nuv gt 0.014)and(disp_nuv lt 0.028) then disp_nuv = 0.02544
if (disp_nuv gt 0.028)and(disp_nuv lt 0.06) then disp_nuv = 0.05088
if (disp_nuv gt 0.06)and(disp_nuv lt .11) then disp_nuv = 0.10176
if (disp_nuv gt 0.11) then begin
   print, 'The NUV dispersion is larger than 10 mA/px which is unusual !'
   print, 'Please check the average profile.'
   stop
endif   
;----------------------------------------------------------------------

disp_nuv *= 100.   ; pm/px
disp_fuv *= 100.
disp3 = disp3 * 100.
disp4 = disp4 * 100.

print,'+++++++++++++++++++++++++++++++++++++++'
print, 'spectral dispersion NUV/FUV = ', disp_nuv, '     / ', disp_fuv, ' pm/px'
print, 'my estimated disp   NUV/FUV = ', disp4, '     / ', disp3
print,'+++++++++++++++++++++++++++++++++++++++'
print,'spectral ranges were selected !        '
print,'+++++++++++++++++++++++++++++++++++++++'
wait, 1.
;----------------------------------------------------------------------
print
print, 'now select ranges for the spectral lines ...'

CII_range = [10, 660]  

OI_range = [662, 1460]
ClI_range = [679, 1060]

OIV_range = [1890,2190]
OIV_range1 = [1890, 2190]  
OIV_range2 = [1890, 2190] 

SiIV_range1 = [1490, 1840]
SiIV_range2 = [1890, 2500]
SiIV_range = [1490, 2500] 

Mg_range = [2900, 3510]  
ct_range = [2540, 2770]

;######################################################################################
if (num eq 1) then begin  ; one average profile for the whole data
vv = avprof
u1 = fxpar(hed, 'WSTART1')
v1 = fxpar(hed, 'WWIDTH1') + u1

if (u1 ne 0)and(v1 ne 0) then vv = avprof[u1:v1]
print, '-----------------------------------'
print, 'select ROUGH C II range'
read_plot_click, 2, vv, vv, pos
pos[0] = pos[0] > 0
pos[1] = pos[1] < (n_elements(vv)-1)
vv = avprof[pos[0]:pos[1]]  &    s = pos[0]
print, 'select exact C II range'
read_plot_click, 2, vv, vv, pos
CII_range = [pos[0], pos[1]] + s

vv = avprof
print, '-----------------------------------'
print, 'select ROUGH O I + C I range (lines are blended)'
read_plot_click, 2, vv, vv, pos
vv = avprof[pos[0]:pos[1]]  &    s = pos[0]
print, 'select exact O I + C I range (lines are blended)'
read_plot_click, 2, vv, vv, pos
OI_range = [pos[0], pos[1]] + s
if ((pos[1]-pos[0]) lt 20) then OI_range=[0,0]
vv = avprof
print, '-----------------------------------'
print, 'select ROUGH strong Cl I + ultra weak Fe XII'
read_plot_click, 2, vv, vv, pos
vv = avprof[pos[0]:pos[1]]  &    s = pos[0]
print, 'select exact Cl I range'
read_plot_click, 2, vv, vv, pos
ClI_range = [pos[0], pos[1]] + s
if ((pos[1]-pos[0]) lt 20) then ClI_range=[0,0]
vv = avprof
print, '-----------------------------------'
print, 'select ROUGH Si IV 1394 range'
read_plot_click, 2, vv, vv, pos
vv = avprof[pos[0]:pos[1]]  &    s = pos[0]
print, 'select exact Si IV 1394 range'
read_plot_click, 2, vv, vv, pos
SiIV_range1 = [pos[0], pos[1]] + s
if ((pos[1]-pos[0]) lt 20) then SiIV1_range=[0,0]
vv = avprof
print, '-----------------------------------'
print, 'select ROUGH Si IV 1403 and O IV'
read_plot_click, 2, vv, vv, pos
vv = avprof[pos[0]:pos[1]]  &    s = pos[0]
print, 'select exact Si IV 1403 and O IV range'
read_plot_click, 2, vv, vv, pos
SiIV_range2 = [pos[0], pos[1]] + s
if ((pos[1]-pos[0]) lt 20) then SiIV2_range=[0,0]

if ((SiIV_range1[0] ne 0)and(SiIV_range2[1] ne 0)) then SiIV_range = [SiIV_range1[0], SiIV_range2[1]]

vv = avprof
print, '-----------------------------------'
print, 'select rough Mg range. Mg II k core should be at px=110'
read_plot_click, 2, vv, vv, pos
pos[0] = pos[0] > 0
pos[1] = pos[1] < (n_elements(vv)-1)
vv = avprof[pos[0]:pos[1]]    &    s = pos[0]
print, 'select exact Mg range'
read_plot_click, 2, vv, vv, pos
pos = pos + s
vv = avprof[pos[0]:pos[1]]    &    s = pos[0]
print, '+++++++++++++++++++++++++++++++++++++++++++++++++++++'
print, 'Click exactly on Mg k3 and then  on any point'
print, '+++++++++++++++++++++++++++++++++++++++++++++++++++++'
read_plot_click, 2, vv, vv, p
pos[0] = pos[0] - (110 - p[0])
Mg_range = [pos[0], pos[1]]

vv = avprof
print, '-----------------------------------'
print, 'select rough 280 nm photospheric range'
read_plot_click, 2, vv, vv, pos
vv = avprof[pos[0]:pos[1]]     &    s = pos[0]
print, 'select exact photospheric range'
read_plot_click, 2, vv, vv, pos
ct_range = [pos[0], pos[1]]  + s

save, filename=outfile, hed, avprof, disp_nuv, disp_fuv, disp3, disp4, ct_range, Mg_range, $
      SiIV_range, SiIV_range1, SiIV_range2, ClI_range, OI_range, CII_range

endif else begin

  f_d = read_dir(polis_path , filter='aver_prof*CII*.sav')
  restore, polis_path + f_d.files[0]
  tay = size(avprof)
  CII_range = [0, tay[1]-1]

  f_d = read_dir(polis_path , filter='aver_prof*SiIV1403*.sav')
  h = check_if_file_exists(polis_path + f_d.files[0])
  if (h[0] eq '      1') then begin
     restore, polis_path + f_d.files[0]
     tay = size(avprof)
     SiIV_range2 = [0, tay[1]-1]
  endif   
  SiIV_range = SiIV_range2
  
  f_d = read_dir(polis_path , filter='aver_prof*SiIV1394*.sav')
  if (h[0] eq '      1') then begin
     restore, polis_path + f_d.files[0]
     tay = size(avprof)
     SiIV_range1 = [0, tay[1]-1]  
     SiIV_range =  [0, SiIV_range1[1] +  SiIV_range2[1]-1]
  endif

  f_d = read_dir(polis_path , filter='aver_prof*OI1356*.sav')
  restore, polis_path + f_d.files[0]
  tay = size(avprof)
  OI_range = [0, tay[1]-1]

  f_d = read_dir(polis_path , filter='aver_prof*FeXII1349*.sav')
  h = check_if_file_exists(polis_path + f_d.files[0])
  if (h[0] eq '      1') then begin
     restore, polis_path + f_d.files[0]
     tay = size(avprof)
     ClI_range = [0, tay[1]-1]
  endif
  
  f_d = read_dir(polis_path , filter='aver_prof*MgIIk2796*.sav')
  restore, polis_path + f_d.files[0]
  tay = size(avprof)
  Mg_range = [0, tay[1]-1]

  f_d = read_dir(polis_path , filter='aver_prof*2832*.sav')
  restore, polis_path + f_d.files[0]
  tay = size(avprof)
  ct_range = [0, tay[1]-1]

restore, outfile
save, filename=outfile, hed, avprof, disp_nuv, disp_fuv, disp3, disp4, ct_range, Mg_range, $
      SiIV_range, SiIV_range1, SiIV_range2, ClI_range, OI_range, CII_range


endelse   

iris_qs_display, polis_path

end  
