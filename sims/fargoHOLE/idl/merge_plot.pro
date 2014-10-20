pro merge_plot, cases=cases, legend=legend, label=label, xrange=xrange,xtickinterval=xtickinterval, modes=modes $
                ,yrange=yrange, scale=scale, xminor=xminor, ytickinterval=ytickinterval, rvq=rvq $
                ,avgrange=avgrange, start=start

cases=string(cases)
numcases=n_elements(cases)



 ;  temp=READ_ASCII(filepath('vortex_time.dat',root_dir='.',subdir=[strcompress(cases(0),/remove_all)]))
;   temp=double(temp.(0))
;   n_modes = n_elements(modes)

;   n = n_elements(temp(0,*))
;   aa = min(where(temp(0,*) ge avgrange))

;   timeavg = dblarr(n_modes + 1, n-aa)

;   for i = 0, n-aa-1 do begin
;       timeavg(0,i) = temp(0,i+aa)
;       for j = 1, n_modes do timeavg(j, i) = mean(temp(modes(j-1) + 1,aa:aa+i))
;   endfor

;   set_plot, 'ps'
;   device, filename='merge_plot.ps' $
;   ,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches

; ;  plot,temp(0,*),alog10(temp(modes(0)+1,*)),xmargin=[8,2],ymargin=[3,1] $
;   plot,timeavg(0,*),alog10(timeavg(1,*)),xmargin=[8,2],ymargin=[3,1] $
;     ,ytitle='log(amplitude)',xtitle='t/orbits' $
;     ,charsize=1.5,xminor=4, thick=4, xrange=xrange,xtickinterval=xtickinterval $
;     ,yrange=yrange, ytickinterval=1

;   for i = 2, n_modes do begin
; ;      oplot,temp(0,*),alog10(temp(modes(i)+1,*)),thick=4,linestyle=i
;      oplot, timeavg(0,*),alog10(timeavg(i,*)), thick=4,linestyle=i-1

;   endfor

; ;  for i=1, numcases-1 do begin
; ;      temp=READ_ASCII(filepath('vortex_time.dat',root_dir='.',subdir=[strcompress(cases(i),/remove_all)]))
; ;      temp=double(temp.(0))
; ;      oplot,temp(0,*),alog10(temp(mode+1,*)),thick=4,linestyle=i
; ;  endfor




;   if keyword_set(legend) then begin
;   x0=legend(0)
;   x1=legend(1)
;   y0=legend(2)
;   dy=legend(3)
;   for j=0, numcases-1 do begin
;       oplot, [x0,x1], [y0,y0]-dy*j, thick=4, linestyle=j
;       xyouts, x1, y0-dy*j,label(j),charsize=1.5
;   endfor
;   endif

;   device,/close




; temp=READ_ASCII(filepath('merging_mode.dat',root_dir='.',subdir=[strcompress(cases(0),/remove_all)]))

 temp=READ_ASCII(filepath('vortex_time.dat',root_dir='.',subdir=[strcompress(cases(0),/remove_all)]))
 temp=double(temp.(0))

 num = n_elements(temp(0,*))
 time_avg = dblarr(num)
 for j = 0, num-1 do time_avg(j) = mean(temp(modes, 0:j))

; array = dindgen(modes) + 1
; aa = min(abs(temp(0,*) - start), range1)
; time = temp(0, range1)

 set_plot, 'ps'
 device, filename='merge_mode.ps' $
   ,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches

 plot,temp(0,*), alog10(time_avg) ,xmargin=[8,2],ymargin=[3,2] $
   ,ytitle=textoidl('log |A_2|'),xtitle='t/orbits' $
   ,charsize=1.5, thick=4, xrange=xrange,xtickinterval=xtickinterval $
   ,yrange=yrange,xminor=xminor,ytickinterval=ytickinterval ;$
   ;,title=string(time,format='(F7.2)')+' orbits'

 for i=1, numcases-1 do begin
     temp=READ_ASCII(filepath('vortex_time.dat',root_dir='.',subdir=[strcompress(cases(i),/remove_all)]))
     ;temp=READ_ASCII(filepath('merging_mode.dat',root_dir='.',subdir=[strcompress(cases(i),/remove_all)]))
      temp=double(temp.(0))
      num = n_elements(temp(0,*))
      time_avg = dblarr(num)
      for j = 0, num-1 do time_avg(j) = mean(temp(modes, 0:j))

;      aa = min(abs(start - temp(0,*)), range1)
      oplot,temp(0,*) ,alog10(time_avg),thick=4,linestyle=i
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


 ; if keyword_set(rvq) then begin
;      n_modes = n_elements(modes)
    
;      set_plot, 'ps'
;      device, filename='merge_plotQ.ps' $
;        ,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches
    
;      avgq = dblarr(numcases, n_modes + 1)
;      avgq(0:numcases-1 , 0) = rvq               ;qmin values
    
    
;      for i = 0, numcases-1 do begin
;          temp=READ_ASCII(filepath('vortex_time.dat',root_dir='.',subdir=[strcompress(cases(i),/remove_all)]))
;          temp=double(temp.(0))
;          aa = min(abs(temp(0,*)- avgrange(0)), range1)
;          aa = min(abs(temp(0,*)- avgrange(1)), range2)
         
;          for j = 1, n_modes do begin
;              avgq(i, j) = mean(temp(modes(j-1), range1:range2))
;          endfor
;      endfor

;      plot,avgq(*,0), alog10(avgq(*, 1)),xmargin=[8,2],ymargin=[3,2] $
;        ,ytitle='log(amplitude)',xtitle=textoidl('Q_m') $
;        ,charsize=1.5, thick=4, xrange=xrange,xtickinterval=xtickinterval $
;        ,yrange=yrange,xminor=xminor,ytickinterval=ytickinterval
    
;      for i=2, n_modes do begin
;          oplot,avgq(*,0), alog10(avgq(*,i)),thick=4,linestyle=i-1
;      endfor
       
;      if keyword_set(legend) then begin
;          x0=legend(0)
;          x1=legend(1)
;          y0=legend(2)
;          dy=legend(3)
;          for j=0, numcases-1 do begin
;              oplot, [x0,x1], [y0,y0]-dy*j, thick=4, linestyle=j
;              xyouts, x1, y0-dy*j,label(j),charsize=1.5
;          endfor
;      endif
;      device,/close
;  endif


end
