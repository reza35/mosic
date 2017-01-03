PRO LPFF, arr,pos,WL=wl,PIX=pix,VEL=vel,PLOT=plot, HELP=help
;+
; Procedure:         LPFF
;
; Purpose:           Measures the position of a spectral line (absorption
;                    or emission within the interval specified.
;
; Category:          Spectral Analysis
;
; Calling sequence:  LPFF,interval,pos[,wl=wl,pix=pix,vel=vel,plot=plot]  
;
; Input:             arr: 1-dim array containing the line profile
;
; optional Input:    wl:  Wavelength in nm
;                    pix: pixel size in pm
;
; optional keyword:  plot: ...
;
; Output:            pos: position of line at subpixel accuracy,
;                         measured from left boundary of 'arr'
;                         [pos=1 means: line core is located at arr(1)]
;
; optional output:   vel: Doppler velocity in m/s, measured from
;                         center of interval, Redshift is positive
;
; Restrictions:     The interval may be much  wider than the line 
;                   or contain only the line core. The routine works
;                   for odd and even numbers of elements of 'arr'.
;                   the line core intensity is NOT measured.
;
; Procedure:        The Fourier phase method is used. This method is 
;                   very insensitive to noise and is very fast. 
;-
; History:
; 1997-Jan-?? ws@kis: written
; 1998-???-?? ws@kis: updated.
; 2000-Jul-24 nlte@kis: documentation.
;
; 
if n_params() lt 1 or n_params() gt 2 or keyword_set(help) then begin
   doc_library,'lpff'
   return
endif
; --------------------------------------------------------------------------
on_error,2
;
if n_elements(wl)  le 0 then wl=500.         ; 500 nm wavelength
if n_elements(pix) le 0 then pix=0.01         ; 10 pm pixel size
sz=size(arr)
length=sz(3)
mid=sz(3)/2.
wll=wl*1.0E-9                      ;wl  in m
pixx=pix*1.0E-12                   ;pix in m
c=2.99E8                           ;c in m/s
tpi=2.*3.141592654
dv=c*pixx/wll                     ;v0/pix in m/s
dp=360./length                    ;Grad/pixel
l1=double(arr)
fl1=fft(l1,-1)
lp=-atan(imaginary(fl1(1))/float(fl1(1)))/tpi*360.
pos= lp/dp + mid
vel=float(lp/dp*dv)
return
end
