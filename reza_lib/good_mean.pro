function good_mean, inp

;+
;===============================================================
; function : good_mean.pro
;
; purpose : to calculate a safe mean of the input array without using outliers
;  
; Dec 28, 2015 : created
;
; R.Rezaei @ IAC                         e-mail:  rrezaei@iac.es      
;===============================================================
;-

   vm = reform(inp)
   e = percentiles(vm, value=[0.2, 0.5, 0.8])
   if (e[0] ne e[2]) then begin
     q = where((vm gt e[0])and(vm lt e[2]), count)
     if (count gt 0) then mean_value = mean(vm(q)) else  mean_value = e[1]
   endif else begin
        e = percentiles(vm, value=[0.01, 0.5, 0.99])
        if (e[0] ne e[2]) then begin
           q = where((vm gt e[0])and(vm lt e[2]), count)
           if (count gt 0) then mean_value = mean(vm(q)) else  mean_value = e[1]
        endif else begin
           mean_value = mean(inp)
        endelse
   endelse    
   return, mean_value

end  
