pro stdev_despike, data, box, stdev_mult
on_error, 2

; use this program to despike your data based on the standard deviation
; of the data within a small local window surrounding the current point.

sizedata=size(data)
numrng=sizedata(1)
numaz=sizedata(2)
workdata=data


;read,'Enter x and y size for local window (use odd) ',xwin,ywin
xwin = box
ywin = box
;print,' '
;print,'Point is a spike if its value > window average + multiple*standard deviation'
;read,'Multiple of Standard Deviation to be considered a spike ',stdev_mult
;print,'....working....be patient....'

tot_spike=long(0)
tot_zero=long(0)
for y=ywin/2,numaz-1-(ywin/2) do begin
    ;if ( y mod 100 eq 0) then print, '...working on line ...',y
    for x=xwin/2,numrng-1-(xwin/2) do begin

        sample = workdata[(x-xwin/2):(x+xwin/2),(y-ywin/2):(y+ywin/2)]
        sample[box/2, box/2] = sample[0,0]
        xystdev = stddev(sample)
        bmean = mean(sample)

           if (abs(workdata[x,y]) gt abs(abs(bmean)+stdev_mult*xystdev)) then begin
               data(x,y) = sample[0,0]
               tot_spike = tot_spike+1
           endif

    endfor
endfor

;print, 'Total Spikes and Windows with all Zeros ',tot_spike,tot_zero


return
end
