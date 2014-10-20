pro orbit, loc=loc, r0=r0 $
             , xtickinterval=xtickinterval, ytickinterval=ytickinterval $
             , xrange = xrange, yrange = yrange, mp=mp, erange = erange $
             , irange=irange

torb = 2d0*!dpi*r0^(3d0/2d0)
bigM = 1d0 + mp


location =strcompress(loc,/remove_all)
nlines = file_lines(filepath('planetxy_dt.dat',root_dir='.',subdir=[location]))
tq = dblarr(7,nlines)
openr, 1, filepath('planetxy_dt.dat',root_dir='.',subdir=[location])
readf,1,tq
close,1

xtitle = 'orbits'

time = tq(0,*)/torb
x = tq(1,*)
y = tq(2,*)
z = tq(3,*)
r = sqrt(x^2 + y^2 + z^2)

vx = tq(4,*)
vy = tq(5,*)
vz = tq(6,*)
vel = sqrt(vx^2 + vy^2 + vz^2)

h1 = y*vz - z*vy
h2 = z*vx - x*vz
h3 = x*vy - y*vx
angmom = sqrt(h1^2 + h2^2 + h3^2)


incl      = atan(sqrt(h1^2 +h2^2)/h3)
semimajor = bigM*r/(2d0*bigM - r*vel^2)
eccen     = sqrt(1d0 - angmom^2/(bigM*semimajor)) 

set_plot, 'ps'
device, filename=filepath(strcompress('orbit_rp.ps',/remove_all) $
                          ,root_dir='.',subdir=[location]) $
  ,/color, bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches
plot, time, r, thick=4, charsize=1.5 $
  ,ymargin=[3,1],xmargin=[6,2],xrange=xrange,yrange=yrange, xtickinterval=xtickinterval $
  , ytickinterval=ytickinterval, ytitle=textoidl('r_p(t)'), xtitle='orbits'
device,/close

set_plot, 'ps'
device, filename=filepath(strcompress('orbit_a.ps',/remove_all) $
                          ,root_dir='.',subdir=[location]) $
  ,/color, bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches
plot, time, semimajor, thick=4, charsize=1.5 $
  ,ymargin=[3,1],xmargin=[6,2],xrange=xrange,yrange=yrange, xtickinterval=xtickinterval $
  , ytickinterval=ytickinterval, ytitle=textoidl('a(t)'), xtitle='orbits'
device,/close


set_plot, 'ps'
device, filename=filepath(strcompress('orbit_e.ps',/remove_all) $
                          ,root_dir='.',subdir=[location]) $
  ,/color, bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches
plot, time, eccen, thick=4, charsize=1.5 $
  ,ymargin=[3,1],xmargin=[8,2],xrange=xrange,yrange=erange, xtickinterval=xtickinterval $
  , ytickinterval=ytickinterval, ytitle=textoidl('e(t)'), xtitle='orbits'
device,/close


set_plot, 'ps'
device, filename=filepath(strcompress('orbit_i.ps',/remove_all) $
                          ,root_dir='.',subdir=[location]) $
  ,/color, bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches
plot, time, incl, thick=4, charsize=1.5 $
  ,ymargin=[3,1],xmargin=[8,2],xrange=xrange,yrange=irange, xtickinterval=xtickinterval $
  , ytickinterval=ytickinterval, ytitle=textoidl('i(t)'), xtitle='orbits'
device,/close


set_plot, 'ps'
device, filename=filepath(strcompress('orbit_ae.ps',/remove_all) $
                          ,root_dir='.',subdir=[location]) $
  , xsize=8, ysize=4.5, xoffset=0, yoffset=0, /inches,/color,bits_per_pixel=8
plot, time,  semimajor, xmargin=[8,8], ymargin=[3,1], xtitle='t/orbits', ytitle='e',/nodata $
  ,charsize=1.5,thick=4,linestyle=0,ystyle=4,xrange=xrange, xtickinterval=xtickinterval
axis,yaxis=0,charsize=1.5,ytitle=textoidl('a (solid)'),/save,yrange=[10,11] ;[min(semimajor),max(semimajor)]
oplot, time, semimajor, linestyle=0, thick=4
axis,yaxis=1,charsize=1.5,ytitle='e (dotted)',/save,yrange=[0.0,0.1] ;[min(eccen),max(eccen)]
oplot,time, eccen, linestyle=1, thick=4
device,/close

device,/close











end

