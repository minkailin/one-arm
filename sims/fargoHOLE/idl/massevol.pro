PRO massevol, cases=cases, legend=legend, label=label, yrange=yrange, xrange=xrange $
, xtickinterval=xtickinterval, ytickinterval=ytickinterval,xminor=xminor

cases=string(cases)
numcases=n_elements(cases)
p0=2.0*!dpi*5.0^(1.5)

temp=READ_ASCII(filepath('planet0.dat',root_dir='.',subdir=[strcompress(cases(0),/remove_all)]))
temp=double(temp.(0))

set_plot, 'ps'
device, filename='massevol.ps' $
,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches

meff = (1.0+2.0*temp(9,*)/temp(5,*))*temp(5,*)

plot,temp(7,*)/p0,meff,xmargin=[8,2],ymargin=[3,1] $
  ,ytitle='effective planet mass',xtitle='t/orbits' $
  ,charsize=1.5, thick=4, xrange=xrange,xtickinterval=xtickinterval $
  ,yrange=yrange,ytickinterval=ytickinterval,xminor=xminor

for i=1, numcases-1 do begin
    temp=READ_ASCII(filepath('planet0.dat',root_dir='.',subdir=[strcompress(cases(i),/remove_all)]))
    temp=double(temp.(0))
    meff = (1.0+2.0*temp(9,*)/temp(5,*))*temp(5,*)
    oplot,temp(7,*)/p0,meff,thick=4,linestyle=i
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


END
