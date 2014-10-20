pro pdisk_kerz_cases, cases=cases, start=start, xrange=xrange, yrange=yrange, xtickinterval=xtickinterval, ytickinterval=ytickinterval, $
                         label=label, legend=legend, title=title, mode=mode, avg=avg
  if not keyword_set(title) then title=''   
 
  xtitle = textoidl('t/P_0')
  ytitle = textoidl('|W_'+strcompress(string(mode),/remove_all)+'|')
  
  filename = strcompress('pdisk_kerz.dat',/remove_all)
  fileloc  = filepath(filename,root_dir='.',subdir=[cases(0)])
  
  nlines = file_lines(fileloc)
  array = dblarr(3,nlines)
  
  openr,1, fileloc
  readf,1, array
  close,1

  if keyword_set(avg) then begin
         xaxis = interpol(array(0,*), avg*nlines,/spline)
         data  = interpol(array(1,*), avg*nlines,/spline)
         data2 = interpol(array(2,*), avg*nlines,/spline) 
;          xaxis = array(0,*)
;          data  = array(1,*)
          
;          for i=0, nlines-1 do begin
;          data(i) = mean(data(0:i)) 
;          endfor 
        

  endif else begin
         xaxis = array(0,*)
         data  = array(1,*)
         data2  = array(2,*)
  endelse 

  
  set_plot, 'ps'
  device, filename=strcompress('pdisk_kerz_cases.ps',/remove_all) $
          ,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches
  plot, xaxis, data,xmargin=[8.5,1.5],ymargin=[3.4,1.6], ystyle=0, xstyle=1  $
        ,charsize=1.5, thick=8, xrange=xrange, xtitle=xtitle, yrange=yrange,$
        linestyle = 0, ytitle =ytitle, xtickinterval=xtickinterval, ytickinterval=ytickinterval 
;  oplot, xaxis, data2, thick=4
  for n=1, n_elements(cases) - 1 do begin
     fileloc  = filepath(filename,root_dir='.',subdir=[cases(n)])
     nlines = file_lines(fileloc)
     array = dblarr(3,nlines)
     openr,1, fileloc
     readf,1, array
     close,1

     if keyword_set(avg) then begin
         xaxis = interpol(array(0,*), avg*nlines,/spline)
         data  = interpol(array(1,*), avg*nlines,/spline)
         data2 = interpol(array(2,*), avg*nlines,/spline)
;           xaxis = array(0,*)
;           data  = array(1,*)
;
;          for i=0, nlines-1 do begin
;          data(i) = mean(data(0:i))
;          endfor

    endif else begin
         xaxis = array(0,*)
         data  = array(1,*)
         data2 = array(2,*)
     endelse

     oplot, xaxis, data, linestyle=n, thick=8
;     oplot, xaxis, data2, linestyle=n, thick=4
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
     xyouts, x1, y0+dy,textoidl(title),charsize=1.5
  endif
  
  device,/close 
end
