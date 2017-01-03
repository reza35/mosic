pro plot_histogram, data, nbins = nbins, xrange = xrange, oplott=oplott
;+
;===============================================================
; procedure :  plot_histogram.pro
; 
; purpose :  it is a simplistic 1D histogram program which plots the PDF of the
; input variable
;  
; June 01, 2016 : created  
;  
; R.Rezaei @ IAC                         e-mail:  rrezaei@iac.es      
;===============================================================
;-
  
if n_elements(data) eq 0 then return
if n_elements(xthik) eq 0 then xthik = 1
if n_elements(ythik) eq 0 then ythik = 1
if n_elements(chtk) eq 0 then chtk=1
if n_elements(nbins) eq 0 then nbins = 50.
if n_elements(xrange) eq 0 then xrange=[min(data), max(data)]
  
hs = histogram(data, nbins=nbins, min=xrange[0], max=xrange[1]) * 1.0
xs = findgen(nbins) * (xrange[1] - xrange[0])/nbins + xrange[0]
hs /=  int_tabulated(xs, hs)
tay = size(xs)
xs = [xs, xs[tay[1]-1]]
hs = [hs, 0.]

if (n_elements(oplott) eq 0) then begin
   plot, xs, hs * 100., psym = 10, ytitle = 'fraction [%]'
endif  else begin
   loadct, 40, /silent
   oplot, xs, hs * 100., psym = 10, color=oplott
   loadct, 0, /silent
endelse


end  
  
