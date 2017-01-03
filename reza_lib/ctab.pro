pro ctab

;+
;===============================================================
; procedure: ctab.pro
;
; purpose : change the color table to blue-red, suitable for velocity.
;
; Jan 21, 2011 : created
;
; R.Rezaei @ KIS                         e-mail:  rrezaei@iac.es      
;===============================================================
;-       

g = [findgen(128)*2, reverse(findgen(128)*2)]
r = [intarr(128)+255, reverse(findgen(128)*2)]
b = [findgen(128)*2, intarr(128)+255]
tvlct, b, g, r

end
