function check_input_profile, wng, mas_k3, mas_h1a, mas_h1b, inp

;+
;===============================================================
; function :  check_input_profile
;
; purpose: to check if the input Mg II h/k profile is alright
; 
;  
; May 28, 2016 : created
;  
; R.Rezaei @ IAC                         e-mail:  rrezaei@iac.es      
;===============================================================
;-  

  tay = size(inp)
  c1 = 0
  q1 = where(~finite(inp), c1)
  ;if (c1 ge 1) then inp(c1) = 0.
  ;c1 = 0
  c2 = 0
  if (median(inp[mas_h1a:mas_h1b]) lt 1.) then c2 = 1 
  ;q2 = where((inp/wng lt 0.)or(inp/wng gt 1d6), c2)
  ;if (c2 gt 0) then  c2 = 1   
  c3 = 0
  if (wng lt 0.02) then c3 = 1
  ;d = median(inp, 5)
  c4 = 0
  if (inp[mas_k3] lt 2.) then c4=1
  ;q1 = where(abs(inp - d) gt 80., c4)
  ;if (c4 gt 0) then begin
  ;   q2 = mean(inp(q1))
  ;   if (q2 lt 1.) then c4 = 0 else stop
  ;endif   
  ;print, c1, c2, c3
  return, c1+c2+c3+c4
end  
