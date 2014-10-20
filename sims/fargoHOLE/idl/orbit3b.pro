PRO orbit3b, cases=cases, ntime=ntime, legend=legend, label=label, xrange=xrange, yrange=yrange $
            ,ytickinterval=ytickinterval, xtickinterval=xtickinterval, mp=mp

bigM = mp + 1d0
cases=string(cases)
numcases=n_elements(cases)
lines = file_lines(filepath('bigplanet0.dat',root_dir='.',subdir=[strcompress(cases(0),/remove_all)]))
array=dblarr(11, lines)
openr,1,filepath('bigplanet0.dat',root_dir='.',subdir=[strcompress(cases(0),/remove_all)])
readf,1,array
close,1

time = array(7,*)
x    = array(1,*)
y    = array(2,*)
vx   = array(3,*) 
vy   = array(4,*)
plrad= sqrt(x*x + y*y)
vel = sqrt(vx*vx + vy*vy)
time /= 2d0*!dpi*plrad(0)^1.5

h3 = x*vy - y*vx
angmom = h3
semimajor = bigM*plrad/(2d0*bigM - plrad*vel^2)
eccen     = sqrt(1d0 - angmom^2/(bigM*semimajor))


set_plot, 'ps'
device, filename='orbit3b_a_e.ps', xsize=8, ysize=4.5, xoffset=0, yoffset=0, /inches,/color,bits_per_pixel=8
plot, time,  eccen, xmargin=[6,8], ymargin=[3,1], xtitle='t/orbits', ytitle='e',/nodata $
,charsize=1.5,thick=4,linestyle=0,ystyle=4,xrange=xrange, xtickinterval=xtickinterval
axis,yaxis=0,charsize=1.5,ytitle=textoidl('a'),/save,yrange=yrange;[min(semimajor),max(semimajor)]
oplot, time,semimajor, linestyle=0, thick=4
axis,yaxis=1,charsize=1.5,ytitle='e (dotted)',/save,yrange=[0.0,0.1];[min(eccen),max(eccen)]
oplot,time, eccen, linestyle=1, thick=4
device,/close






;;;;;;;;;;
;PLOTTING;
;;;;;;;;;;
set_plot, 'ps'
device, filename='orbit3b_rp.ps', xsize=8, ysize=4.5, xoffset=0, yoffset=0, /inches,/color,bits_per_pixel=8
plot, time, plrad, xmargin=[6,2], ymargin=[3,2], xtitle='t/orbits' $
, ytitle=textoidl('r_p(t)'), yrange=yrange, ytickinterval=ytickinterval $
, charsize=1.5,yminor=10,thick=4,linestyle=0, xrange=xrange, xtickinterval=xtickinterval, xminor=5

for k=1, numcases-1 do begin

lines = file_lines(filepath('bigplanet0.dat',root_dir='.',subdir=[strcompress(cases(k),/remove_all)]))
array=dblarr(11, lines)
openr,1,filepath('bigplanet0.dat',root_dir='.',subdir=[strcompress(cases(k),/remove_all)])
readf,1,array
close,1

time = array(7,*)
x    = array(1,*)
y    = array(2,*)
plrad= sqrt(x*x + y*y)
time /= 2d0*!dpi*plrad(0)^1.5

oplot, time, plrad,thick=4,linestyle=k
endfor

if keyword_set(legend) then begin
x0=legend(0)
x1=legend(1)
y1=legend(2)
dy=legend(3)
for k=0, numcases-1 do begin
oplot,[x0,x1],[y1,y1]-k*dy,thick=4,linestyle=k
xyouts,x1,y1-k*dy,textoidl(label(k)),charsize=1.5
endfor
endif
 device,/close
END
