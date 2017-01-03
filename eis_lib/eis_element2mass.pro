
;+
; NAME:
;
; eis_element2mass
;
; PURPOSE:
;
; Returns the mass (in grams or kilograms) of an element.
;
; CATEGORY:
;
; Data analysis.
;
; CALLING SEQUENCE:
;
; mass = eis_element2mass(element)
; mass = eis_element2mass(element,/kg)
; 
; INPUTS:
;
; The element's symbol (case sensitive).
;
; OPTIONAL INPUTS:
;
; None.
;       
; KEYWORD PARAMETERS:
;
; kg: Mass is output in kilograms rather than grams for use with MKS.
;
; OUTPUTS:
;
; The mass of the element in grams or kilograms.
;
; OPTIONAL OUTPUTS:
;
; None.
;
; COMMON BLOCKS:
;
; None.
;
; SIDE EFFECTS:
;
; None.
;
; RESTRICTIONS:
;
; Only the first 30 elements are available. Returns -1 if the input
; element is not found.
;
; PROCEDURE:
;
; Uses a simple lookup table.
;
; EXAMPLE:
;
; IDL> Fe_mass = element2mass('Fe')
; IDL> print,Fe_mass
;   9.27339e-23
;
; MODIFICATION HISTORY:
;
;   HPW 16-AUG-1996: Based on a routine from John Mariska's SSWHAT
;                    program.
;   HPW 15-MAY-2007: Renamed eis_element2mass.
;   HPW 10-AUG-2007: Added STRUPCASE to make more robust.
;
;-

function eis_element2mass,element,kg=kg
  
amu2g = 1.6605E-24
data  = replicate({datum,element:	'',$
                         z:		 0,$
                         amu:	       0.0},30)

data(0)  = {datum,"H",   1, 1.0080}
data(1)  = {datum,"He",  2, 4.0026}
data(2)  = {datum,"Li",  3, 6.941}
data(3)  = {datum,"Be",  4, 9.0122}
data(4)  = {datum,"B",   5, 10.811}
data(5)  = {datum,"C",   6, 12.0111}
data(6)  = {datum,"N",   7, 14.0067}
data(7)  = {datum,"O",   8, 15.9994}
data(8)  = {datum,"F",   9, 18.9984}
data(9)  = {datum,"Ne", 10, 20.179}
data(10) = {datum,"Na", 11, 22.9898}
data(11) = {datum,"Mg", 12, 24.305}
data(12) = {datum,"Al", 13, 26.9815}
data(13) = {datum,"Si", 14, 28.086}
data(14) = {datum,"P",  15, 30.9738}
data(15) = {datum,"S",  16, 32.06}
data(16) = {datum,"Cl", 17, 35.453}
data(17) = {datum,"Ar", 18, 39.948}
data(18) = {datum,"K",  19, 39.102}
data(19) = {datum,"Ca", 20, 40.08}
data(20) = {datum,"Sc", 21, 44.956}
data(21) = {datum,"Ti", 22, 47.90}
data(22) = {datum,"V",  23, 50.9414}
data(23) = {datum,"Cr", 24, 51.996}
data(24) = {datum,"Mn", 25, 54.9380}
data(25) = {datum,"Fe", 26, 55.847}
data(26) = {datum,"Co", 27, 58.9332}
data(27) = {datum,"Ni", 28, 58.71}
data(28) = {datum,"Cu", 29, 63.546}
data(29) = {datum,"Zn", 30, 65.37}

m = where(STRUPCASE(data.element) eq STRUPCASE(element),count)
if (count gt 0) then begin
  mass = data(m).amu*amu2g 
  if keyword_set(kg) then mass = mass/1000.0
endif else mass = -1

return,mass
end
