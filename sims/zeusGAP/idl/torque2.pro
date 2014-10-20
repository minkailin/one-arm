pro torque2, loc=loc, r0=r0, scale=scale $
             , xtickinterval=xtickinterval, ytickinterval=ytickinterval $
             , xrange = xrange, yrange = yrange, mp=mp, legend = legend $
             , yyrange=yyrange, twodim=twodim, name=name

torb = 2d0*!dpi*r0^(3d0/2d0)
if not keyword_set(scale) then scale = 1d0

location =strcompress(loc,/remove_all)
nlines = file_lines(filepath('dptorque.dat',root_dir='.',subdir=[location]))
tq = dblarr(8,nlines)
openr, 1, filepath('dptorque.dat',root_dir='.',subdir=[location])
readf,1,tq
close,1

xtitle = textoidl('t/P_0')

time = tq(0,*)/torb
tq_in = tq(1,*)*scale
tq_out= tq(2,*)*scale
tq_tot= tq(3,*)*scale

tq_in_ex  = tq(4,*)*scale
tq_out_ex = tq(5,*)*scale
tq_tot_ex = tq(6,*)*scale

hmass = tq(7,*)/mp

if keyword_set(legend) then begin
    label = ['total','inner','outer']

    x0 = legend(0)
    x1 = legend(1)
    y0 = legend(2)
    dy = legend(3)
endif


set_plot, 'ps'
device, filename=filepath(strcompress('torque2_tq.ps',/remove_all) $
                          ,root_dir='.',subdir=[location]) $
  ,/color, bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches
plot, time, tq_tot, thick=4, charsize=1.5 $
  ,ymargin=[3.5,0.5],xmargin=[8,2],xrange=xrange,yrange=yrange, xtickinterval=xtickinterval $
  , ytickinterval=ytickinterval, ytitle='Torque', xtitle=xtitle
oplot, time, tq_in, thick=4, linestyle=1
oplot, time, tq_out, thick=4, linestyle=2

if keyword_set(legend) then begin
    for j=0, 2 do begin
        oplot, [x0, x1], [y0, y0]- j*dy, thick=4, linestyle=j
        xyouts, x1, y0-j*dy, label(j), charsize=1.5
    endfor
endif
device,/close

set_plot, 'ps'
device, filename=filepath(strcompress('torque2_tqex.ps',/remove_all) $
                          ,root_dir='.',subdir=[location]) $
  ,/color, bits_per_pixel=8,xsize=8, ysize=4. ,xoffset=0,yoffset=0,/inches
plot, time, tq_tot_ex, thick=4 ,ymargin=[3.5,0.5],xmargin=[8,2],xrange=xrange $
  ,yrange=yrange, xtickinterval=xtickinterval, charsize=1.5 $
  , ytickinterval=ytickinterval, ytitle='Torque (tapered)', xtitle=xtitle
oplot, time, tq_in_ex, thick=4, linestyle=1
oplot, time, tq_out_ex, thick=4, linestyle=2

if keyword_set(legend) then begin
    for j=0, 2 do begin
        oplot, [x0, x1], [y0, y0]- j*dy, thick=4, linestyle=j
        xyouts, x1, y0-j*dy, label(j), charsize=1.5
     endfor
    if keyword_set(name) then xyouts, x0, y0+dy, name, charsize=1.5
endif
device,/close

set_plot, 'ps'
device, filename=filepath(strcompress('torque2_mhill.ps',/remove_all) $
                          ,root_dir='.',subdir=[location]) $
  ,/color, bits_per_pixel=8,xsize=8, ysize=4.5 ,xoffset=0,yoffset=0,/inches
plot, time, hmass, thick=4, charsize=1.5  $
  ,ymargin=[3,1],xmargin=[8,2],xrange=xrange,yrange=yyrange, xtickinterval=xtickinterval $
  , ytickinterval=ytickinterval, ytitle=textoidl('m_H/M_p'), xtitle=xtitle
device,/close

stop

end

