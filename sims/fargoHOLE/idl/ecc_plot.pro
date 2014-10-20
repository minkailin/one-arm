pro ecc_plot, cases=cases, legend=legend, label=label, xrange=xrange,xtickinterval=xtickinterval,yrange=yrange

cases=string(cases)
numcases=n_elements(cases)

temp=READ_ASCII(filepath('ecc_time.dat',root_dir='.',subdir=[strcompress(cases(0),/remove_all)]))
temp=double(temp.(0))


nel = n_elements(temp(0,*))
;tmp = min(abs(temp(0,*) - 50.0), beg)
;for i=beg,  nel-1 do temp(1,i) = mean(temp(1,beg:i))

beg=0

set_plot, 'ps'
device, filename='ecc_plot.ps' $
,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches

plot,temp(0,beg:nel-1),temp(1,beg:nel-1),xmargin=[6,2],ymargin=[3,2] $
  ,ytitle=textoidl('e_{edge}'),xtitle='t/orbits' $
  ,charsize=1.5,xminor=5, thick=4, xrange=xrange,xtickinterval=xtickinterval $
  ,yrange=yrange
;oplot,temp(0,beg:nel-1),temp(2,beg:nel-1)/!dpi, linestyle=1


for i=1, numcases-1 do begin

    temp=READ_ASCII(filepath('ecc_time.dat',root_dir='.',subdir=[strcompress(cases(i),/remove_all)]))
    temp=double(temp.(0))

    nel = n_elements(temp(0,*))
;    tmp = min(abs(temp(0,*) - 50.0), beg)
;    for j=beg,  nel-1 do temp(1,j) = mean(temp(1,beg:j))
    oplot,temp(0,beg:nel-1),temp(1,beg:nel-1),thick=4,linestyle=i
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

end
