pro vline, x, cc, dashed=dashed
;+
;===============================================================
; procedure: vline.pro
;
; purpose : overplot a vertical line. 
;
; Sep 17, 2006 : created
;
; R.Rezaei @ KIS                         e-mail:  rrezaei@iac.es      
;===============================================================
;-       

if (n_elements(dashed) eq 0) then dashed=0 else dashed=2  
if (n_elements(cc) eq 0) then begin
  oplot,[x,x],[-1.0d4,1.0d4], linestyle=dashed
endif else begin
  loadct, 40, /silent
  oplot,[x,x],[-1.0d4,1.0d4], color=cc, linestyle=dashed
  loadct, 0, /silent
endelse

end
