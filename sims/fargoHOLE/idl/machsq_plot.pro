pro machsq_plot, cases=cases, legend=legend, label=label, xrange=xrange,xtickinterval=xtickinterval $
                 ,yrange=yrange, scale=scale, xminor=xminor, ytickinterval=ytickinterval

cases=string(cases)
numcases=n_elements(cases)

temp=READ_ASCII(filepath('vortensity_jump.dat',root_dir='.',subdir=[strcompress(cases(0),/remove_all)]))
temp=double(temp.(0))

set_plot, 'ps'
device, filename='machsq_plot.ps' $
,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches

plot,temp(0,*),temp(2,*),xmargin=[8,2],ymargin=[3,2] $
  ,ytitle=textoidl('vortensity jump'),xtitle='r-r_p' $
  ,charsize=1.5, thick=4, xrange=xrange,xtickinterval=xtickinterval $
  ,yrange=yrange,xminor=xminor,ytickinterval=ytickinterval,title='100 orbits'

 ;oplot ,temp(0,*),temp(3,*), thick=2

for i=1, numcases-1 do begin
    temp=READ_ASCII(filepath('vortensity_jump.dat',root_dir='.',subdir=[strcompress(cases(i),/remove_all)]))
    temp=double(temp.(0))
    ;temp(1,1:modes-1) *= 1d2
    oplot,temp(0,*), temp(2,*),thick=4,linestyle=i
    ;oplot,temp(0,*),temp(3,*), thick=2, linestyle=i
endfor

if keyword_set(legend) then begin
x0=legend(0)
x1=legend(1)
y0=legend(2)
dy=legend(3)
for j=0, numcases-1 do begin
    oplot, [x0,x1], [y0,y0]-dy*j, thick=4, linestyle=j
    xyouts, x1, y0-dy*j,label(j),charsize=1.5
endfor
endif

device,/close

end
