pro read_plot_click, np, ff, ffn, x_lobes

;+
;===============================================================
; procedure: read_plot_click.pro
;
; purpose: reads the click position on a line profile.
;
; np : number of points to read/click
; ff : profile
; ffn : another profile, if desired. It will be overplotted.
; x_lobes : output variable  
;
; Aug 12, 2014, 2007 : created  
;
; R.Rezaei @ KIS                         e-mail:  rrezaei@iac.es      
;===============================================================
;-       

x_lobes = 0.
npr = n_elements(ff)
if n_elements(x_lobes) ne np then begin
   print,' Please select '+strtrim(string(np),1)+' points manually:'
   !p.multi = 0
   plot,ff,xr =[0, npr],/xstyle, psym=-1   
   if (npr gt 0) then  loadct,40,/silent & oplot, ffn, color=245 & loadct,0,/silent
   ilobeindex = 0.
   x_lobes = fltarr(np)
   while ilobeindex ne np do begin
      cursor,x,y,/data,/wait
      if !mouse.button eq 1 then begin
         x_lobes(ilobeindex) = x > 0 < (npr-1)
         ilobeindex = ilobeindex+1
         print, x
         oplot,[x,x],[-10,4d4]
         wait,.5
      endif
      if ilobeindex eq np then begin
         print, 'Selection okay ? (y/n) '
         ans=get_kbrd()
         if ans eq 'n' then begin
            ilobeindex = 0
            plot,ff,xr =[0, npr],/xstyle, psym=-1  
            loadct,40,/silent & oplot, ffn, color=245 & loadct,0,/silent
         endif
      endif
   endwhile
   x_lobes = x_lobes(0:ilobeindex-1)
   x_lobes = x_lobes(sort(x_lobes))
endif

x_lobes = round(x_lobes)
x_lobes = x_lobes(sort(x_lobes))
print, x_lobes

end
