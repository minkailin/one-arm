pro torque_corot_plot, loc=loc, legend=legend, label=label, xrange=xrange,xtickinterval=xtickinterval $
                       ,yrange=yrange, scale=scale, ytickinterval=ytickinterval, start=start, finish=finish $
                       , out=out, segment=segment, mp=mp, title=title
f0 = (mp/3d0)^(1d0/3d0)

rmin = -4d0
rmax =  4d0

nr = 128
rplot = rmin + (rmax-rmin)*dindgen(nr)/(nr-1d0)
torque_spline = dblarr(nr)

if not keyword_set(scale) then scale = 1d0

;;;;;;;;;;;;;;;
;WHERE IS DATA;
;;;;;;;;;;;;;;;
location=strcompress(loc,/remove_all)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;GET DIMENSIONS AND TIME UNIT;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
dims=(read_ascii(filepath('dims.dat',root_dir='.',subdir=[location]))).(0)
nout=fix(dims(5))
nrad=fix(dims(6))
nsec=fix(dims(7))
if not keyword_set (out) then begin
    info=dblarr(11,nout+1)
endif else info=dblarr(11,out+1)
openr,3,filepath('planet0.dat',root_dir='.',subdir=[location])
readf,3,info
close,3
a0=info(1,0)
dt=info(7,1)
period=2.*!dpi*(a0)^(3./2.)

azi = -1d0 + 2d0*dindgen(nsec)/(nsec - 1d0)

torque1d_rad_input = dblarr(2, nrad)
torque1d_azi_input = dblarr(2, nsec)

torque1d_rad_alltime = dblarr(finish-start+1,nr) 
torque1d_azi_alltime = dblarr(finish-start+1,nsec)

torque1d_rad_avg = dblarr(2, nr)
torque1d_azi_avg = dblarr(2, nsec)


for k=start, finish do begin
    ks=string(k,format='(I03)')

    openr,1,filepath(strcompress('torque1d_rad_'+ks+'.dat',/remove_all),root_dir='.',subdir=[location])
    readf,1,torque1d_rad_input
    close,1
    
;    openw,1,filepath(strcompress('torque1d_rad_nohill'+ks+'.dat',/remove_all),root_dir='.',subdir=[location])
;    for i = 0, nrad-1 do printf,1,rplot(i),torque1d_rad_nohill(i)
;    close,1

    openr,1,filepath(strcompress('torque1d_azi_'+ks+'.dat',/remove_all),root_dir='.',subdir=[location])
    readf,1,torque1d_azi_input
    close,1
    
;    openw,1,filepath(strcompress('torque1d_azi_nohill'+ks+'.dat',/remove_all),root_dir='.',subdir=[location])
;    for j = 0, nsec-1 do printf,1,azi1(j),torque1d_azi(j)
;    close,1

    
    torque_spline = spline(torque1d_rad_input(0,*), torque1d_rad_input(1,*), rplot, /double)
    torque1d_rad_alltime(k-start,0:nr-1) = torque_spline(0:nr-1)


    torque1d_azi_alltime(k-start, 0:nsec-1) = torque1d_azi_input(1,0:nsec-1)


endfor

for i=0, nr -1 do begin
    torque1d_rad_avg(0,i) = mean(torque1d_rad_alltime(segment(0)-start:segment(1)-start,i))
    torque1d_rad_avg(1,i) = mean(torque1d_rad_alltime(segment(1)-start:segment(2)-start,i))
endfor

for j=0, nsec-1 do begin
    torque1d_azi_avg(0,j) = mean(torque1d_azi_alltime(segment(0)-start:segment(1)-start,j))
    torque1d_azi_avg(1,j) = mean(torque1d_azi_alltime(segment(1)-start:segment(2)-start,j))
endfor

set_plot, 'ps'
device, filename=filepath(strcompress('torque_corot_plot_rad.ps',/remove_all),root_dir='.',subdir=[location]) $
  ,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches

plot,rplot,torque1d_rad_avg(0,*)*scale,xmargin=[8,2],ymargin=[3.5,2.0] $
  ,ytitle=textoidl('Torque'),xtitle=textoidl('(r-r_p)/r_h'), title=textoidl(title) $
  ,charsize=1.5,xminor=4, thick=4, xrange=xrange,xtickinterval=xtickinterval $
  ,yrange=yrange, ytickinterval=ytickinterval
oplot, rplot, torque1d_rad_avg(1,*)*scale, thick=4, linestyle=1
oplot, [1,1], [-1,1]*1d10
oplot, [-1,-1], [-1,1]*1d10

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

set_plot, 'ps'
device, filename=filepath(strcompress('torque_corot_plot_azi.ps',/remove_all),root_dir='.',subdir=[location]) $
  ,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches

plot,azi,abs(torque1d_azi_avg(0,*))*scale,xmargin=[6,2],ymargin=[3.5,2.0] $
  ,ytitle=textoidl('|Torque|'),xtitle=textoidl('(\phi-\phi_p)/\pi'), title=textoidl(title) $
  ,charsize=1.5,xminor=4, thick=4, xrange=xrange,xtickinterval=xtickinterval $
  ,yrange=yrange
oplot, azi, abs(torque1d_azi_avg(1,*))*scale, thick=4, linestyle=1
oplot, [f0,f0]/!dpi, [-1,1]*1d10
oplot, [-f0,-f0]/!dpi, [-1,1]*1d10

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
