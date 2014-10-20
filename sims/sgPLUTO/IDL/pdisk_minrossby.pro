pro pdisk_minrossby, cases=cases, xrange=xrange, yrange=yrange, xtickinterval=xtickinterval, ytickinterval=ytickinterval, $
                          label=label, legend=legend, title=title
  if not keyword_set(title) then title=''   
 
  xtitle = textoidl('t/P_0')
  ytitle=textoidl('min(Ro)')
  
  filename = strcompress('minrossby.dat',/remove_all)
  
  fileloc  = filepath(filename,root_dir='.',subdir=[cases(0)])
  
  nlines = file_lines(fileloc)
  array = dblarr(2,nlines)
  
  openr,1, fileloc
  readf,1, array
  close,1
  
  set_plot, 'ps'
  device, filename=strcompress('pdisk_minrossby.ps',/remove_all) $
          ,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches
  plot, array(0,*), array(1,*),xmargin=[8.5,1.5],ymargin=[3.4,1.6], ystyle=0, xstyle=1  $
        ,charsize=1.5, thick=4, xrange=xrange, title=textoidl(title), xtitle=xtitle, yrange=yrange,$
        linestyle = 0, ytitle =ytitle, xtickinterval=xtickinterval, ytickinterval=ytickinterval
  for n=1, n_elements(cases) - 1 do begin
     filename = strcompress('minrossby.dat',/remove_all)
     fileloc  = filepath(filename,root_dir='.',subdir=[cases(n)])
     nlines = file_lines(fileloc)
     array = dblarr(2,nlines)
     openr,1, fileloc
     readf,1, array
     close,1
     oplot, array(0,*), array(1,*), linestyle=n, thick=4
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
