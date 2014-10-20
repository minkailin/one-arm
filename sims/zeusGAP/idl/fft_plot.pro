pro fft_plot, cases=cases, mtot=mtot, mode=mode, legend=legend, label=label $
              ,xtickinterval=xtickinterval,xrange=xrange,yrange=yrange $
              ,ytickinterval=ytickinterval

ncase = n_elements(cases)

location =strcompress(cases(0),/remove_all)
nlines = file_lines(filepath('fft_time.dat',root_dir='.',subdir=[location]))
data = dblarr(mtot + 1,nlines)
openr, 1, filepath('fft_time.dat',root_dir='.',subdir=[location])
readf,1, data
close,1

ytitle = strcompress('log_{10}|C_'+string(mode)+'/C_0|',/remove_all)
;ytitle = strcompress('|C_'+string(mode)+'/C_0|',/remove_all)
ytitle = textoidl(ytitle)

data(1:mtot,*) = alog(data(1:mtot,*))

set_plot, 'ps'
device, filename='fft_mode.ps' $
  ,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches
plot, data(0,*), data(mode,*),xmargin=[8,2],ymargin=[3,1], ystyle=0  $
  ,charsize=1.5, thick=4, xrange=xrange, yrange=yrange, xtitle=textoidl('t/P_0'), ytickinterval=ytickinterval  $ 
  ,xtickinterval=xtickinterval, ytitle=ytitle

for k=1, ncase-1 do begin
    location =strcompress(cases(k),/remove_all)
    nlines = file_lines(filepath('fft_time.dat',root_dir='.',subdir=[location]))
    data = dblarr(mtot + 1,nlines)
    openr, 1, filepath('fft_time.dat',root_dir='.',subdir=[location])
    readf,1, data
    close,1    
    data(1:mtot,*) = alog(data(1:mtot,*))
    oplot, data(0,*), data(mode,*),thick=4,linestyle=k
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

