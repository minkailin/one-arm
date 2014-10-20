function timeavg, input, time, begavg=begavg
num = n_elements(input)
output = dblarr(num)

temp = min(abs(time - begavg),grid)

for i=0l+grid, num-1 do output(i) = mean(input(grid:i))
return, output
end

PRO torque2, cases=cases, ntime=ntime, legend=legend, label=label, yrange=yrange, xrange=xrange $
             , xtickinterval=xtickinterval, ytickinterval=ytickinterval, scale=scale, begavg=begavg $
             ,xminor=xminor, rp0=rp0

if not keyword_set(begavg) then begavg = 30.0
if not keyword_set(scale) then scale = 1d0

period = 2.0*!dpi*rp0^1.5

cases=string(cases)
numcases=n_elements(cases)

temp=READ_ASCII(filepath('tqwk0.dat',root_dir='.',subdir=[strcompress(cases(0),/remove_all)]))
temp=double(temp.(0))
num=n_elements(temp(0,*))
temp(9,*) /= period

inner_torque_notaper = temp(1,*)*scale
outer_torque_notaper = temp(2,*)*scale
torque_notaper = (temp(2,*) + temp(1,*))*scale
torque_notaper_avg = timeavg(torque_notaper,temp(9,*),begavg=begavg)
torque_outer_notaper_avg = timeavg(outer_torque_notaper,temp(9,*),begavg=begavg)

inner_torque_taper = temp(3,*)*scale
outer_torque_taper = -temp(4,*)*scale
torque_taper = (temp(3,*) + temp(4,*))*scale
torque_taper_avg = timeavg(torque_taper,temp(9,*),begavg=begavg)

;;;;;;;;;;
;PLOTTING;
;;;;;;;;;;
set_plot, 'ps'
device, filename='torque2_.ps', xsize=8, ysize=4.5, xoffset=0, yoffset=0, /inches,/color,bits_per_pixel=8
plot, temp(9,*), torque_notaper, xmargin=[8,2], ymargin=[3,1], xtitle='t/orbits' $
  , ytitle=textoidl('Torque'), yrange=yrange , ytickinterval=ytickinterval $
  , charsize=1.5,thick=4,linestyle=0, xrange=xrange,xtickinterval=xtickinterval,xminor=xminor

oplot, temp(9,*), torque_taper, thick=4, linestyle=
;oplot, temp(9,*), inner_torque_taper, thick=4,linestyle=1
;oplot, temp(9,*), -outer_torque_taper, thick=4,linestyle=2


for k=1, numcases-1 do begin
    temp=READ_ASCII(filepath('tqwk0.dat',root_dir='.',subdir=[strcompress(cases(k),/remove_all)]))
    temp=double(temp.(0))    
    temp(9,*) /= period

    inner_torque_notaper = temp(1,*)*scale
    outer_torque_notaper = temp(2,*)*scale
    torque_notaper = (temp(2,*) + temp(1,*))*scale
    torque_notaper_avg = timeavg(torque_notaper,temp(9,*),begavg=begavg)
    torque_outer_notaper_avg = timeavg(outer_torque_notaper,temp(9,*),begavg=begavg)

    inner_torque_taper = temp(3,*)*scale
    outer_torque_taper = -temp(4,*)*scale
    torque_taper = (temp(3,*) + temp(4,*))*scale
    torque_taper_avg = timeavg(torque_taper,temp(9,*),begavg=begavg)
    
;    oplot, temp(9,*), torque_notaper, thick=4,linestyle=k
;    oplot, temp(9,*), torque_taper_avg, thick=1,linestyle=k
    oplot, temp(9,*), torque_notaper , thick=4,linestyle=k
endfor

if keyword_set(legend) then begin
    x0=legend(0)
    x1=legend(1)
    y1=legend(2)
    dy=legend(3)
    for k=0, n_elements(label)-1 do begin
        oplot,[x0,x1],[y1,y1]-k*dy,thick=4,linestyle=k
        xyouts,x1,y1-k*dy,textoidl(label(k)),charsize=1.5
    endfor
endif
device,/close









END

