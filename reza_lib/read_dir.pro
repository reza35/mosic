function read_dir, polis_path, filter=filter
;+
;===============================================================
; function : read_dir.pro
; 
; purpose :  read files in the given directory.
;  
; June 01, 2016 : created  
;  
; R.Rezaei @ IAC                         e-mail:  rrezaei@iac.es      
;===============================================================
;-
  
  on_error, 2
  cd, current = dir1
  cd, polis_path


  if (n_elements(filter) eq 0) then filter = ' '
  list = file_search(filter)

  list_out = {files:list}
  cd, dir1
  return, list_out
end
