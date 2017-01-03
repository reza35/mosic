function num2string, inp
;+
;===============================================================
; function :  num2string.pro
; 
; purpose : to convert a number to a string
;  
; June 01, 2016 : created  
;  
; R.Rezaei @ IAC                         e-mail:  rrezaei@iac.es      
;===============================================================
;-
  return, strtrim(string(inp),2)
end
