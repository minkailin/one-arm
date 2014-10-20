pro pattern_plot, cases=cases, legend=legend, label=label, xrange=xrange,xtickinterval=xtickinterval, modes=modes $
                  ,yrange=yrange, qmin=qmin, xminor=xminor, ytickinterval=ytickinterval, rvq=rvq, avgrange=avgrange $
                  , time = time 

cases=string(cases)
numcases=n_elements(cases)

temp=READ_ASCII(filepath('mode_amplitudes.dat',root_dir='.',subdir=[strcompress(cases(0),/remove_all)]))
temp=double(temp.(0))

;if(numcases eq 1) then begin
    
    set_plot, 'ps'
    device, filename='pattern_plot.ps' $
      ,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches
    

    slice = min(abs(time(0) - temp(0,*)), grid)

    plot,dindgen(15)+1, alog10(temp(1:2, grid)),xmargin=[8,2],ymargin=[3,2] $
      ,ytitle='log Amplitude',xtitle='m' $
      ,charsize=1.5, thick=4, xrange=xrange,xtickinterval=xtickinterval $
      ,yrange=yrange,xminor=xminor,ytickinterval=ytickinterval;,title=textoidl('Q_m=') + string(qmin, format='(f3.1)')
    
;    n_modes = n_elements(modes)

    for i=1, numcases-1 do begin
        temp=READ_ASCII(filepath('mode_amplitudes.dat',root_dir='.',subdir=[strcompress(cases(i),/remove_all)]))
        temp=double(temp.(0))
        slice = min(abs(time(i) - temp(0,*)), grid)

        oplot, dindgen(15)+1, alog10(temp(1:2, grid)),thick=4,linestyle=i
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
;endif

    temp=READ_ASCII(filepath('mode_amplitudes.dat',root_dir='.',subdir=[strcompress(cases(0),/remove_all)]))
    temp=double(temp.(0))

    set_plot, 'ps'
    device, filename='pattern_plot2.ps' $
      ,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches
    
    num = n_elements(temp(0,*))
    time_avg = dblarr(num)
    for j = 0, num-1 do time_avg(j) = mean(temp(2,0:j))
    
    plot,temp(0,*) , alog10(time_avg),xmargin=[8,2],ymargin=[3,2] $
      ,ytitle=textoidl('log |A|'),xtitle='t/orbits' $
      ,charsize=1.5, thick=6, xrange=xrange,xtickinterval=xtickinterval $
      ,yrange=yrange,xminor=xminor,ytickinterval=ytickinterval
    
    for j = 0, num-1 do time_avg(j) = mean(temp(3,0:j))
    oplot, temp(0,*), alog10(time_avg), thick=3

    for i=1, numcases-1 do begin
        temp=READ_ASCII(filepath('mode_amplitudes.dat',root_dir='.',subdir=[strcompress(cases(i),/remove_all)]))
        temp=double(temp.(0))

        num = n_elements(temp(0,*))
        time_avg = dblarr(num)
        for j = 0, num-1 do time_avg(j) = mean(temp(2,0:j))
        oplot, temp(0,*), alog10(time_avg),thick=6,linestyle=i


        if(i le 2) then begin
            for j = 0, num-1 do time_avg(j) = mean(temp(3,0:j))
            oplot, temp(0,*), alog10(time_avg),thick=3,linestyle=i
        endif

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











; if(numcases gt 1) then begin
    
;     if not keyword_set(rvq) then begin
;         set_plot, 'ps'
;         device, filename='pattern_plot.ps' $
;           ,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches
        
;         plot,temp(0,*), temp(modes + 1, *),xmargin=[8,2],ymargin=[3,2] $
;           ,ytitle='growth rate',xtitle=' time/orbits' $
;           ,charsize=1.5, thick=4, xrange=xrange,xtickinterval=xtickinterval $
;           ,yrange=yrange,xminor=xminor,ytickinterval=ytickinterval,title=strcompress("m="+string(modes),/remove_all)
        
;         for i=1, numcases-1 do begin
;             temp=READ_ASCII(filepath('growth_rates.dat',root_dir='.',subdir=[strcompress(cases(i),/remove_all)]))
;             temp=double(temp.(0))
;             oplot,temp(0,*), temp(modes+1, *),thick=4,linestyle=i
;         endfor
        
;         if keyword_set(legend) then begin
;             x0=legend(0)
;             x1=legend(1)
;             y0=legend(2)
;             dy=legend(3)
;             for j=0, numcases-1 do begin
;                 oplot, [x0,x1], [y0,y0]-dy*j, thick=4, linestyle=j
;                 xyouts, x1, y0-dy*j,label(j),charsize=1.5
;             endfor
;         endif
;         device,/close
;     endif
    
;     if keyword_set(rvq) then begin
;     n_modes = n_elements(modes)
    
;     set_plot, 'ps'
;     device, filename='pattern_plotQ.ps' $
;       ,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches
    
;     avgq = dblarr(numcases, n_modes + 1)
;     avgq(0:numcases-1 , 0) = rvq               ;qmin values
    
    
;     for i = 0, numcases-1 do begin
;         temp=READ_ASCII(filepath('growth_rates.dat',root_dir='.',subdir=[strcompress(cases(i),/remove_all)]))
;         temp=double(temp.(0))
;         aa = min(abs(temp(0,*)- avgrange(0)), range1)
;         aa = min(abs(temp(0,*)- avgrange(1)), range2)
;         for j = 1, n_modes do begin
;             avgq(i, j) = mean(temp(modes(j-1) + 1, range1:range2))*1d2
;         endfor
;     endfor

    

;     plot,avgq(*,0), avgq(*, 1),xmargin=[8,2],ymargin=[3,2] $
;       ,ytitle='average growth rate',xtitle=textoidl('Q_m') $
;       ,charsize=1.5, thick=4, xrange=xrange,xtickinterval=xtickinterval $
;       ,yrange=yrange,xminor=xminor,ytickinterval=ytickinterval
    
;     for i=2, n_modes do begin
;         oplot,avgq(*,0), avgq(*,i),thick=4,linestyle=i-1
;     endfor
       
;     if keyword_set(legend) then begin
;         x0=legend(0)
;         x1=legend(1)
;         y0=legend(2)
;         dy=legend(3)
;         for j=0, numcases-1 do begin
;             oplot, [x0,x1], [y0,y0]-dy*j, thick=4, linestyle=j
;             xyouts, x1, y0-dy*j,label(j),charsize=1.5
;         endfor
;     endif
;     device,/close
; endif


; endif

end
