function timeavg, input, time, begavg=begavg
num = n_elements(input)
output = dblarr(num)

temp = min(abs(time - begavg),grid)

for i=0l+grid, num-1 do output(i) = mean(input(grid:i))
return, output 
end


pro  gap_profile_plot, cases=cases, legend=legend, label=label, xrange=xrange,xtickinterval=xtickinterval,yrange=yrange, start=start $
                    ,ytickinterval=ytickinterval, begavg=begavg

if not keyword_set(begavg) then begavg = 30d0

cases=string(cases)
numcases=n_elements(cases)

lines = file_lines(filepath('gap_profile.dat',root_dir='.',subdir=[strcompress(cases(0),/remove_all)]))
array = dblarr(5, lines)
openr,1,filepath('gap_profile.dat',root_dir='.',subdir=[strcompress(cases(0),/remove_all)])
readf,1,array
close,1
time = array(0,*)
;data =0.5*( array(3,*) + array(1,*))
;data = array(4,*) - array(2,*)
;data = array(3,*)
;data = array(4,*)
data = array(4,*) + array(2,*)
;data = array(1,*)
;data = abs(array(2,*))
data_avg = timeavg(data, time, begavg=begavg)

;ytitle = 'Outer gap depth'
;ytitle = 'Inner gap depth'
;ytitle = 'Outer gap width'
;ytitle = 'Inner gap width'
;ytitle = 'Gap width'
;ytitle = 'Gap depth'
ytitle= 'Gap asymmetry'


set_plot, 'ps'
device, filename='gap_profile_plot.ps' $
,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches

plot,time, data_avg,xmargin=[8,2],ymargin=[3.5,1.5] $
  ,ytitle=ytitle,xtitle='t/orbits' $
  ,charsize=1.5,xminor=4, thick=4, xrange=xrange,xtickinterval=xtickinterval $
  ,yrange=yrange,ytickinterval=ytickinterval

for i=1, numcases-1 do begin

    lines = file_lines(filepath('gap_profile.dat',root_dir='.',subdir=[strcompress(cases(i),/remove_all)]))
    array = dblarr(5, lines)
    openr,1,filepath('gap_profile.dat',root_dir='.',subdir=[strcompress(cases(i),/remove_all)])
    readf,1,array
    close,1
    time = array(0,*)
;    data = 0.5*( array(3,*) + array(1,*))
;    data = array(3,*)
;    data = array(4,*)
;    data = array(4,*) - array(2,*)
    data = array(4,*)+ array(2,*)
;    data = array(1,*)
;    data = abs(array(2,*))
    data_avg = timeavg(data, time, begavg=begavg)

    oplot, time, data_avg, thick=4,linestyle=i
endfor

if keyword_set(legend) then begin
    x0=legend(0)
    x1=legend(1)
    y0=legend(2)
    dy=legend(3)
    for j=0, numcases-1 do begin
        oplot, [x0,x1], [y0,y0]-dy*j, thick=4, linestyle=j
        xyouts, x1, y0-dy*j,textoidl(label(j)),charsize=1.5
    endfor
endif
device,/close

end
