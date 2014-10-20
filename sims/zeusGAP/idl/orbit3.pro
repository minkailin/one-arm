pro orbit3, cases=cases, r0=r0 $
             , xtickinterval=xtickinterval, ytickinterval=ytickinterval $
             , xrange = xrange, yrange = yrange, legend=legend, label=label

torb = 2d0*!dpi*r0^(3d0/2d0)
ncase = n_elements(cases)

location =strcompress(cases(0),/remove_all)
nlines = file_lines(filepath('planetxy_dt.dat',root_dir='.',subdir=[location]))
tq = dblarr(7,nlines)
openr, 1, filepath('planetxy_dt.dat',root_dir='.',subdir=[location])
readf,1,tq
close,1

xtitle = textoidl('t/orbits')

time = tq(0,*)/torb
x = tq(1,*)
y = tq(2,*)
z = tq(3,*)
r = sqrt(x^2 + y^2 + z^2)

set_plot, 'ps'
device, filename='orbit3_rp.ps' $
  ,/color, bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches
plot, time, r, thick=4, charsize=1.5 $
  ,ymargin=[3.5,0.5],xmargin=[6,2],xrange=xrange,yrange=yrange, xtickinterval=xtickinterval $
  , ytickinterval=ytickinterval, ytitle=textoidl('r_p(t)'), xtitle=xtitle, xminor=5

for k=1, ncase-1 do begin

    location =strcompress(cases(k),/remove_all)
    nlines = file_lines(filepath('planetxy_dt.dat',root_dir='.',subdir=[location]))
    tq = dblarr(7,nlines)
    openr, 1, filepath('planetxy_dt.dat',root_dir='.',subdir=[location])
    readf,1,tq
    close,1
       
    time = tq(0,*)/torb
    x = tq(1,*)
    y = tq(2,*)
    z = tq(3,*)
    r = sqrt(x^2 + y^2 + z^2)
    
    oplot, time, r,thick=4,linestyle=k
endfor

if keyword_set(legend) then begin
    x0=legend(0)
    x1=legend(1)
    y0=legend(2)
    dy=legend(3)
    for j=0, n_elements(label)-1 do begin
        oplot, [x0,x1], [y0,y0]-dy*j, thick=4, linestyle=j
        xyouts, x1, y0-dy*j,textoidl(label(j)),charsize=1.5
    endfor
endif

device,/close

end

