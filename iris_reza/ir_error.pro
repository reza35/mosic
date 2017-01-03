function ir_error, DN, dark=dark_k $
                 , fuv=fuv, nuv=nuv, sjinuv=sjinuv, sjifuv=sjifuv
;
; calculate error do a given DN
;
;
; PHOTON NOISE:
;
; See De Pontieu et al (2013), beginning of Sect 7.1
;
;------------------------------------------------------------
; NUV spectra AND NUV SJI:
;------------------------------------------------------------
; 18   electrons / DN
; 1 electrons / photons
;
; -> 18 photons = 1 DN
;                      -> sigma_DN = sqrt(DN/18) = 0.24 sqrt(DN)
;
;------------------------------------------------------------
; FUV spectra:
;------------------------------------------------------------
;
; use conversion: 6   electrons / DN
;                 1.5 electrons / photons
;                 
; -> 4 photons = 1 DN
;                      -> sigma_DN = sqrt(DN/4) = 0.5 sqrt(DN)
;
;------------------------------------------------------------
; NUV spectra AND NUV SJI:
;------------------------------------------------------------
; 18   electrons / DN
; 1.5 electrons / photons
;
; -> 12 photons = 1 DN
;                      -> sigma_DN = sqrt(DN/12) = 0.29 sqrt(DN)
;
;
;
;
;------------------------------------------------------------
; READOUT NOISE
;------------------------------------------------------------
;
; 20 electrons
;
; NUV spectra & SJI:   20/18 = 1.2 DN
; FUV spectra:         20/6  = 3.3 DN
; FUV SJI:             20/18 = 1.2 DN
;
;
;
; in many data it looks like the dark is not subtracted correctly.
; Check the continuum level. The lower end of the histogram of the
; continuum intensity would be the "dark" that has not been subtracted
; correctly. E.g. in the quiet Sun data I had a range of the continuum
; level from 16 DN to 20 DN. Thus I subtracted 16 DN from the count rate
; before calculating the errors. This way zero continuum should correspond
; to 0 DN... Set keyword dark to subtract the base continuum level (in the 
; above mentioned example set dark=16).
;
;
; peter@mps.mpg.de Oct 7, 2013
; added dark keyword, Nov 27 2013
;
;
if n_elements(dark_k) eq 0 then dark=16. else dark=float(dark_k)



error=0.

if keyword_set(fuv)    then error = sqrt( float(abs(DN-dark))/4.  + 3.3^2 )
if keyword_set(sjifuv) then error = sqrt( float(abs(DN-dark))/12. + 1.2^2 )

if keyword_set(nuv)    then error = sqrt( float(abs(DN-dark))/18. + 1.2^2 )
if keyword_set(sjinuv) then error = sqrt( float(abs(DN-dark))/18. + 1.2^2 )

return, error

; return, sqrt( 2*float(abs(DN)) + 3.3^2 )
end

