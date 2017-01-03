;             NAME : iris_te
; 
;          PURPOSE : Derive the electron temperature from line pair intensity ratio and plot the result
;
; CALLING SEQUENCE : iris_te,int1,err1,int2,err2,temp_file,rat,rat_err,temp,temp1,temp2
; 
;           INPUTS : int1 - intensity of the first line (total line intensity, continuum subtracted), can be a scaler, a vector, or an array of any dimension
;           
;                    err1 - error of the first line intensity
;           
;                    int2 - intensity of the second line
;                    
;                    err2 - error of the second line intensity
;                    
;                    temp_file - file of the theoretical relationship between temperature and ratio. Here is an example on how to 
;                    generate this file:
;                    
;                    temperature_ratios,'fe_12',195,1350,5.5,7.0,temp,rat,desc,density=10^8.5
;                    t=interpol(alog10(temp),250,/spline)
;                    r=interpol(rat,250,/spline)
;                    save,filename='fe12_195to1349_temp_n8.5.sav',t,r
; 
;         KEYWORDS : posi - position of the plot
;                    
;                    ytitle -  lable of y-cooridinate of the temperature-ratio plot, e.g., 'Ratio of Fe XII 195.12/1349.40'
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
;                    temp - temperature
;                    
;                    temp1 - lower bound of the temperature
;                    
;                    temp2 - upper bound of the temperature
;
;          HISTORY : Written by Hui Tian at CfA, April 9, 2013
;
;; **********************************************************************

Pro iris_te,int1,err1,int2,err2,temp_file,rat,rat_err,temp,temp1,temp2,ytitle=ytitle,posi=posi,xrange=xrange,savedat=savedat,savefig=savefig


!p.font=-1

if not(keyword_set(ytitle)) then ytitle = 'Ratio'
if not(keyword_set(posi)) then posi=[0.10,0.10,0.98,0.95]
if not(keyword_set(xrange)) then xrange=[4.5,7.0]

;theoretical relationship between ratio and temperature
restore,temp_file,/ver

npt=n_elements(int1)

;ratio and error
rat=int1/int2
rat_err=rat*sqrt((err1/int1)^2+(err2/int2)^2) 

;temperature and error
temp=fltarr(npt) & temp1=fltarr(npt) & temp2=fltarr(npt)
for i=0,npt-1 do begin
rmin=min(abs(r - rat[i]),sub) 
temp[i]=t[sub]
rmin=min(abs(r - (rat[i]-rat_err[i])),sub) 
temp1[i]=t[sub]
rmin=min(abs(r - (rat[i]+rat_err[i])),sub) 
temp2[i]=t[sub]
endfor

if npt gt 100 then goto,endprogram

window,0,xs=800,ys=600

;plot the theoretical curve
plot,t,r,xrange=xrange,xstyle=1,position=posi,xtitle='Log (Temperature / K)',ytitle=ytitle,charsize=1.5

;plot the calculated densities and errors
if npt gt 1 then begin
  oplot,temp,rat,psym=6,symsize=0.5,thick=8
  err_plot,temp,rat-rat_err,rat+rat_err,width=0.01,thick=1
  err_plot,rat,temp1,temp2,width=0.02,thick=1,XDir=1
endif else begin
  oplot,[temp,temp],[rat,rat],psym=6,symsize=0.5,thick=8
  err_plot,[temp,temp],[rat,rat]-[rat_err,rat_err],[rat,rat]+[rat_err,rat_err],width=0.01,thick=1
  err_plot,[rat,rat],[temp1,temp1],[temp2,temp2],width=0.02,thick=1,XDir=1
endelse

;save the plot
if keyword_set(savefig) then write_png,'Te.png', tvrd(/true), r, g, b

endprogram:
;save the density and ratio results
if keyword_set(savedat) then save,filename='Te.sav',rat,rat_err,temp,temp1,temp2
end
