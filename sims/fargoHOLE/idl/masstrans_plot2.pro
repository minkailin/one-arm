pro masstrans_plot2, cases=cases, legend=legend, label=label, xrange=xrange,xtickinterval=xtickinterval,yrange=yrange, start=start $ 
                    ,ytickinterval=ytickinterval

cases=string(cases)
numcases=n_elements(cases)

lines = file_lines(filepath('masstrans2_time.dat',root_dir='.',subdir=[strcompress(cases(0),/remove_all)]))
temp0 = dblarr(5,lines)
openr,1, filepath('masstrans2_time.dat',root_dir='.',subdir=[strcompress(cases(0),/remove_all)])
readf,1,temp0, format='(5(e22.15,2x))'
close,1

set_plot, 'ps'
device, filename='masstrans_plot2_flux.ps' $
,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches

plot,temp0(0,*),temp0(2,*)*1d2,xmargin=[6,2],ymargin=[3,2] $
  ,ytitle=textoidl('10^2\Sigmau_r/(\Sigma_0c_s)'),xtitle='t/orbits' $
  ,charsize=1.5,xminor=5, thick=4, xrange=xrange,xtickinterval=xtickinterval $
  ,yrange=yrange,ytickinterval=ytickinterval

for i=1, numcases-1 do begin
    lines = file_lines(filepath('masstrans2_time.dat',root_dir='.',subdir=[strcompress(cases(i),/remove_all)]))
    temp = dblarr(5,lines)
    openr,1, filepath('masstrans2_time.dat',root_dir='.',subdir=[strcompress(cases(i),/remove_all)])
    readf,1,temp, format='(5(e22.15,2x))'
    close,1

    oplot,temp(0,*),temp(2,*)*1d2,thick=4,linestyle=i
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


set_plot, 'ps'
device, filename='masstrans_plot2_gap.ps' $
,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches

plot,temp0(0,*),temp0(1,*),xmargin=[8,2],ymargin=[3,2] $
  ,ytitle=textoidl('-\Delta\Sigma/\Sigma'),xtitle='t/orbits' $
  ,charsize=1.5,xminor=5, thick=4, xrange=xrange,xtickinterval=xtickinterval $
  ,yrange=yrange,ytickinterval=ytickinterval

for i=1, numcases-1 do begin
    lines = file_lines(filepath('masstrans2_time.dat',root_dir='.',subdir=[strcompress(cases(i),/remove_all)]))
    temp = dblarr(5,lines)
    openr,1, filepath('masstrans2_time.dat',root_dir='.',subdir=[strcompress(cases(i),/remove_all)])
    readf,1,temp, format='(5(e22.15,2x))'
    close,1

    oplot,temp(0,*),temp(1,*),thick=4,linestyle=i
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


set_plot, 'ps'
device, filename='masstrans_plot2_gapwidth.ps' $
,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches

plot,temp0(0,*),temp0(3,*),xmargin=[8,2],ymargin=[3,2] $
  ,ytitle=textoidl('w/r_h'),xtitle='t/orbits' $
  ,charsize=1.5,xminor=5, thick=4, xrange=xrange,xtickinterval=xtickinterval $
  ,yrange=yrange,ytickinterval=ytickinterval

for i=1, numcases-1 do begin
    lines = file_lines(filepath('masstrans2_time.dat',root_dir='.',subdir=[strcompress(cases(i),/remove_all)]))
    temp = dblarr(5,lines)
    openr,1, filepath('masstrans2_time.dat',root_dir='.',subdir=[strcompress(cases(i),/remove_all)])
    readf,1,temp, format='(5(e22.15,2x))'
    close,1

    oplot,temp(0,*),temp(3,*),thick=4,linestyle=i
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


set_plot, 'ps'
device, filename='masstrans_plot2_discmass.ps' $
,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches

plot,temp0(0,*),temp0(4,*),xmargin=[8,2],ymargin=[3,2] $
  ,ytitle=textoidl('M_d/M_d(0)'),xtitle='t/orbits' $
  ,charsize=1.5,xminor=5, thick=4, xrange=xrange,xtickinterval=xtickinterval $
  ,yrange=yrange,ytickinterval=ytickinterval

for i=1, numcases-1 do begin
    lines = file_lines(filepath('masstrans2_time.dat',root_dir='.',subdir=[strcompress(cases(i),/remove_all)]))
    temp = dblarr(5,lines)
    openr,1, filepath('masstrans2_time.dat',root_dir='.',subdir=[strcompress(cases(i),/remove_all)])
    readf,1,temp, format='(5(e22.15,2x))'
    close,1

    oplot,temp(0,*),temp(4,*),thick=4,linestyle=i
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



if keyword_set(start) then begin
    ks=string(start,format='(I03)')

    filename = strcompress('masstrans2_'+ks(0)+'.dat',/remove_all)    
    lines = file_lines(filepath(filename,root_dir='.',subdir=[strcompress(cases(0),/remove_all)]))
    temp0 = dblarr(3,lines)
    openr,1, filepath(filename,root_dir='.',subdir=[strcompress(cases(0),/remove_all)])
    readf,1,temp0, format='(3(e22.15,2x))'
    close,1
    
    set_plot, 'ps'
    device, filename='masstrans_plot2_flux1d.ps' $
      ,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches
    
    plot,temp0(0,*), temp0(2,*)*1d3, xmargin=[6,2],ymargin=[3,2] $
      ,ytitle=textoidl('10^3\Sigmau_r/(\Sigma_0v_k)'),xtitle='r' $
      ,charsize=1.5, thick=4,xtickinterval=xtickinterval,ytickinterval=ytickinterval $
      ,xrange=xrange, yrange=yrange, linestyle=0
    
    for i=1, numcases-1 do begin
        filename = strcompress('masstrans2_'+ks(i)+'.dat',/remove_all)    
        lines = file_lines(filepath(filename,root_dir='.',subdir=[strcompress(cases(i),/remove_all)]))
        temp = dblarr(3,lines)
        openr,1, filepath(filename,root_dir='.',subdir=[strcompress(cases(i),/remove_all)])
        readf,1,temp, format='(3(e22.15,2x))'
        close,1
    
        oplot,temp(0,*), temp(2,*)*1d3, thick=4,linestyle=i
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

    
    set_plot, 'ps'
    device, filename='masstrans_plot2_gap1d.ps' $
      ,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches
    
    plot,temp0(0,*), temp0(1,*), xmargin=[8,2],ymargin=[3,2] $
      ,ytitle=textoidl('\Delta\Sigma/\Sigma'),xtitle='r' $
      ,charsize=1.5, thick=4,xtickinterval=xtickinterval,ytickinterval=ytickinterval $
      ,xrange=[3,7], yrange=yrange, linestyle=0
    
    for i=1, numcases-1 do begin
        filename = strcompress('masstrans2_'+ks(i)+'.dat',/remove_all)    
        lines = file_lines(filepath(filename,root_dir='.',subdir=[strcompress(cases(i),/remove_all)]))
        temp = dblarr(3,lines)
        openr,1, filepath(filename,root_dir='.',subdir=[strcompress(cases(i),/remove_all)])
        readf,1,temp, format='(3(e22.15,2x))'
        close,1

        oplot,temp(0,*), temp(1,*), thick=4,linestyle=i
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
endif
end
