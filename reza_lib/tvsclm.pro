pro tvsclm, dat, ilim=ilim, xp=xp, yp=yp, cb=cb, zm=zm

;+
;===============================================================
; pro tvsclm.pro
;
; purpose : make a tvscl considering the outlier values and the screen size
;
; dat : 2D array
; ilim: a vector of 2 numbers for the intensity limit, optional.
; xp, yp : location, optional
; cp : colortable, optional
; zm : zoom factor
;  
; Dec 29, 2015 : original
;  
; R.Rezaei @ IAC                         e-mail:  rrezaei@iac.es      
;===============================================================
;-
  on_error, 2
  
  dim = get_screen_size(RESOLUTION=resolution)
  td = size(dat)
  if (td[1] ge dim[0]/3.)or(td[2] ge dim[1]/2.) then zoom = 2 else zoom = 1

  if (n_elements(zm) ne 0) then zoom=zm
  if (n_elements(xp) eq 0) then xp=0
  if (n_elements(yp) eq 0) then yp=0
  if (n_elements(cb) eq 0) then cb=0
  
  if (n_elements(ilim) eq 0) then begin
     q = where(dat ne 0., count)
     if (count gt 1.) then u = percentiles(dat(q), value=[0.02, 0.98]) else return
  endif else begin
     u = ilim
  endelse   
  
  if (zoom eq 1.) then begin
       if (cb le 40) then loadct, cb, /silent else ctab
       tvscl, dat > u[0] < u[1], xp, yp
       loadct, 0, /silent
  endif else begin
       if (cb le 40) then loadct, cb, /silent else ctab
       tvscl, congrid(dat, td[1]/zoom, td[2]/zoom, /interp) > u[0] < u[1], xp, yp
       loadct, 0, /silent
  endelse   

  
end  
