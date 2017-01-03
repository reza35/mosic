;             NAME : iris_ne_oiv
; 
;          PURPOSE : Derive the electron density from O IV 1399 & 1401 line pair intensity ratio and plot the result
;
; CALLING SEQUENCE : iris_ne_oiv,int1,err1,int2,err2,rat,rat_err,den,den1,den2
; 
;           INPUTS : int1 - intensity of the first line (total line intensity, continuum subtracted), can be a scaler, a vector, or an array of any dimension
;           
;                    err1 - error of the first line intensity
;           
;                    int2 - intensity of the second line
;                    
;                    err2 - error of the second line intensity
; 
;         KEYWORDS : posi - position of the plot
;                    
;                    ytitle -  lable of y-cooridinate of the density-ratio plot, e.g., 'Ratio of Fe XIV 264/274'
;                    
;                    xrange - range of the x-coordinate of the plot
;                     
;                    savedat - save the result
;
;                    savefig - save the plot
;                    
;          OUTPUTS : rat - line ratio
;                    
;                    rat_err - error of the ratio
;                    
;                    den - density
;                    
;                    den1 - lower bound of the density
;                    
;                    den2 - upper bound of the density
;
;          HISTORY : Written by Hui Tian at CfA, Oct 13, 2013
;
;; **********************************************************************

Pro iris_ne_oiv,int1,err1,int2,err2,rat,rat_err,den,den1,den2,ytitle=ytitle,posi=posi,xrange=xrange,savedat=savedat,savefig=savefig


!p.font=-1

if not(keyword_set(ytitle)) then ytitle = 'Ratio of O IV 1399/1401'
if not(keyword_set(posi)) then posi=[0.10,0.10,0.98,0.95]
if not(keyword_set(xrange)) then xrange=[7,13]

;theoretical relationship between ratio and density
rat0=[0.16710843, 0.1670735, 0.16709631, 0.16718450, 0.16736760, 0.16770638, 0.16831132, 0.16937367, 0.17121214,$
      0.17433385, 0.17949909, 0.18777250, 0.20053728, 0.2193597, 0.24536467, 0.2779118, 0.3135369, 0.34697879,$
      0.3740375, 0.39339120, 0.40603911, 0.41382271, 0.4184366, 0.42111109, 0.42264118]
den0= [1.00000e+07, 1.77828e+07, 3.16228e+07, 5.62341e+07, 1.00000e+08, 1.77828e+08, 3.16228e+08, 5.62341e+08, 1.00000e+09,$
      1.77828e+09, 3.16228e+09, 5.62341e+09, 1.00000e+10, 1.77828e+10, 3.16228e+10, 5.62341e+10, 1.00000e+11, 1.77828e+11,$
      3.16228e+11, 5.62341e+11, 1.00000e+12, 1.77828e+12, 3.16228e+12, 5.62341e+12, 1.00000e+13]
d=interpol(alog10(den0),200,/spline)
r=interpol(rat0,200,/spline)


npt=n_elements(int1)

;ratio and error
rat=int1/int2
rat_err=rat*sqrt((err1/int1)^2+(err2/int2)^2) 

;density and error
den=fltarr(npt) & den1=fltarr(npt) & den2=fltarr(npt)
for i=0,npt-1 do begin
rmin=min(abs(r - rat[i]),sub) 
den[i]=d[sub]
rmin=min(abs(r - (rat[i]-rat_err[i])),sub) 
den1[i]=d[sub]
rmin=min(abs(r - (rat[i]+rat_err[i])),sub) 
den2[i]=d[sub]
endfor

if npt gt 100 then goto,endprogram

window,0,xs=800,ys=600

;plot the theoretical curve
plot,d,r,xrange=xrange,xstyle=1,position=posi,xtitle='Log (Density / cm!E-3!N)',ytitle=ytitle,charsize=1.5

;plot the calculated densities and errors
if npt gt 1 then begin
  oplot,den,rat,psym=6,symsize=0.5,thick=8
  err_plot,den,rat-rat_err,rat+rat_err,width=0.01,thick=1
  err_plot,rat,den1,den2,width=0.02,thick=1,XDir=1
endif else begin
  oplot,[den,den],[rat,rat],psym=6,symsize=0.5,thick=8
  err_plot,[den,den],[rat,rat]-[rat_err,rat_err],[rat,rat]+[rat_err,rat_err],width=0.01,thick=1
  err_plot,[rat,rat],[den1,den1],[den2,den2],width=0.02,thick=1,XDir=1
endelse

;save the plot
if keyword_set(savefig) then write_png,'Ne.png', tvrd(/true), r, g, b

endprogram:
;save the density and ratio results
if keyword_set(savedat) then save,filename='Ne.sav',rat,rat_err,den,den1,den2
end
