function check_if_file_exists, infile

;+
;===============================================================
; function :  check_if_file_exists.pro
; 
; purpose : to check if a file exists.
;  
; June 01, 2016 : created  
;  
; R.Rezaei @ IAC                         e-mail:  rrezaei@iac.es      
;===============================================================
;-

  return, [infile, string(file_test(infile))]
end
