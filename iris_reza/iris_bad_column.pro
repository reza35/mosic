pro iris_bad_column, kont, data

;+
;===============================================================
; procedure : iris_bad_column
;  
; purpose: to replace dead columns in some maps with average of left and right columns.
; it is a very simple algorithm and does not fix all of dead pixels !
  
; data : 3D adat cube of the form [lambda, slit, scan]
;
; Note : Obsolete, as we now read the data slit by slit.  
;  
; June 17, 2015: created
;
; R.Rezaei @ IAC                         e-mail:  rrezaei@iac.es      
;===============================================================
;-


vtf = size(data)

b = total(kont, 2)
u = good_mean(b)
for i=1,vtf[3]-2 do begin
   if ((abs(b[i]-u)/u gt 0.9)) then begin
      print, 'dead column = ',i
      data[*,*,i] = (data[*,*,i-1] + data[*,*,i+1])*0.5
      print, '+++++++++++++++++++++++++++++ DEAD COLUMN ++++++++++++++++++++++++++++++++++'
      print
   endif
endfor
; now check the first and last column
i = 0
if ((abs(b[i]-u)/u gt 0.9)) then begin
      print, 'dead column = ',i
      data[*,*,i] = data[*,*,i+1]
      print, '+++++++++++++++++++++++++++++ DEAD COLUMN ++++++++++++++++++++++++++++++++++'
      print
endif
i = vtf[3]-1
if ((abs(b[i]-u)/u gt 0.9)) then begin
      print, 'dead column = ',i
      data[*,*,i] = data[*,*,i-1]
      print, '+++++++++++++++++++++++++++++ DEAD COLUMN ++++++++++++++++++++++++++++++++++'
      print
endif



end
