pro masstrans_plot, cases=cases, legend=legend, label=label, xrange=xrange,xtickinterval=xtickinterval,yrange=yrange, start=start $
                    ,ytickinterval=ytickinterval

cases=string(cases)
numcases=n_elements(cases)

temp=READ_ASCII(filepath('masstrans_ang_time.dat',root_dir='.',subdir=[strcompress(cases(0),/remove_all)]))
temp=double(temp.(0))


nel = n_elements(temp(0,*))
tmp = min(abs(temp(0,*) - 100.0), beg)
for i=beg,  nel-1 do temp(2,i) = mean(temp(2,beg:i))

; beg=0

set_plot, 'ps'
device, filename='masstrans_plot.ps' $
,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches

plot,temp(0,beg:nel-1),temp(2,beg:nel-1),xmargin=[8,2],ymargin=[3,2] $
  ,ytitle=textoidl('Mdot(r=5.5)\times10^{7}/M_d'),xtitle='t/orbits' $
  ,charsize=1.5,xminor=5, thick=4, xrange=xrange,xtickinterval=xtickinterval $
  ,yrange=yrange,ytickinterval=ytickinterval

for i=1, numcases-1 do begin

    temp=READ_ASCII(filepath('masstrans_ang_time.dat',root_dir='.',subdir=[strcompress(cases(i),/remove_all)]))
    temp=double(temp.(0))

    nel = n_elements(temp(0,*))
   tmp = min(abs(temp(0,*) - 100.0), beg)
   for j=beg,  nel-1 do temp(2,j) = mean(temp(2,beg:j))
    oplot,temp(0,beg:nel-1),temp(2,beg:nel-1),thick=4,linestyle=i
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

if keyword_set(start) then begin
    info = dblarr(11,start+1)
    openr,3,filepath('planet0.dat',root_dir='.',subdir=[cases(0)])
    readf,3,info
    close,3
    a0=info(1,0)
    p0=2.*!dpi*(a0)^(3./2.)
    time=string(info(7,start)/p0,format='(F7.2)')


    ks=string(start,format='(I03)')
    filename = strcompress('masstrans_'+ks+'.dat',/remove_all)
    temp=READ_ASCII(filepath(filename,root_dir='.',subdir=[strcompress(cases(0),/remove_all)]))
    temp=double(temp.(0))
    
    set_plot, 'ps'
    device, filename='masstrans_plot1d.ps' $
      ,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches
    
    plot,temp(0,*), temp(1,*), xmargin=[8,2],ymargin=[3,2] $
      ,ytitle=textoidl('M_{dot}\times10^7/M_d'),xtitle='r' $
      ,charsize=1.5, thick=4,xtickinterval=xtickinterval,ytickinterval=ytickinterval $
      ,xrange=xrange, yrange=yrange, title=time+' orbits', linestyle=0
    
    for i=1, numcases-1 do begin
        filename = strcompress('masstrans_'+ks+'.dat',/remove_all)
        temp=READ_ASCII(filepath(filename,root_dir='.',subdir=[strcompress(cases(i),/remove_all)]))
        temp=double(temp.(0))    
        oplot,temp(0,*), temp(1,*), thick=4,linestyle=i
    endfor
    
    ;oplot, [1,10],[0,0]

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

    
    filename = strcompress('masstrans_'+ks+'.dat',/remove_all)
    temp=READ_ASCII(filepath(filename,root_dir='.',subdir=[strcompress(cases(0),/remove_all)]))
    temp=double(temp.(0))
    
    set_plot, 'ps'
    device, filename='masstrans_plot1d_ang.ps' $
      ,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches
    
    plot,temp(0,*), temp(2,*)*1d9, xmargin=[8,2],ymargin=[3,2] $
      ,ytitle=textoidl('<\Sigma\deltau_r\deltau_\phi>\times10^9'),xtitle='r' $
      ,charsize=1.5, thick=4,xtickinterval=xtickinterval,ytickinterval=ytickinterval $
      ,xrange=xrange, yrange=yrange, title=time+' orbits', linestyle=0
    
    for i=1, numcases-1 do begin
        filename = strcompress('masstrans_'+ks+'.dat',/remove_all)
        temp=READ_ASCII(filepath(filename,root_dir='.',subdir=[strcompress(cases(i),/remove_all)]))
        temp=double(temp.(0))    
        oplot,temp(0,*), temp(2,*)*1d9, thick=4,linestyle=i
    endfor
    
    ;oplot, [1,10],[0,0]

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
