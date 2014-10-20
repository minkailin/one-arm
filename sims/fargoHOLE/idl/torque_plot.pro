pro torque_plot, cases=cases, legend=legend, label=label, xrange=xrange,xtickinterval=xtickinterval, mode=mode $
,yrange=yrange, scale=scale, start=start, ytickinterval=ytickinterval, mp=mp, out=out

if not keyword_set(scale) then scale = 1d0

cases=string(cases)
numcases=n_elements(cases)

ks = string(start(0),format='(I03)') 
temp=READ_ASCII(filepath('torque_1d'+ks+'.dat',root_dir='.',subdir=[strcompress(cases(0),/remove_all)]))
temp=double(temp.(0))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;GET DIMENSIONS AND TIME UNIT;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
dims=(read_ascii(filepath('dims.dat',root_dir='.',subdir=[strcompress(cases(0),/remove_all)]))).(0)
nout=fix(dims(5))
nrad=fix(dims(6))
nsec=fix(dims(7))
if not keyword_set (out) then begin
    info=dblarr(11,nout+1)
endif else info=dblarr(11,out+1)
openr,3,filepath('planet0.dat',root_dir='.',subdir=[strcompress(cases(0),/remove_all)])
readf,3,info
close,3

plx=info(1,start(0))
ply=info(2,start(0))
rp = sqrt(plx*plx + ply*ply)
rhill = rp*(mp/3.0)^(1.0/3.0)


set_plot, 'ps'
device, filename='torque_plot.ps' $
  ,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches

rad = (temp(0,*) - rp)/rhill
plot,rad,temp(1,*)*scale,xmargin=[6,2],ymargin=[4,1] $
  ,ytitle=textoidl('Torque'),xtitle=textoidl('(r-r_p)/r_h') $
  ,charsize=1.5,xminor=5, thick=4, xrange=xrange,xtickinterval=xtickinterval $
  ,yrange=yrange, ytickinterval=ytickinterval

for i=1, numcases-1 do begin

    ks = string(start(i),format='(I03)') 
    temp=READ_ASCII(filepath('torque_1d'+ks+'.dat',root_dir='.',subdir=[strcompress(cases(i),/remove_all)]))
    temp=double(temp.(0))

    
    plx=info(1,start(i))
    ply=info(2,start(i))
    rp = sqrt(plx*plx + ply*ply)
    rhill = rp*(mp/3.0)^(1.0/3.0)

    rad = (temp(0,*) - rp)/rhill
    oplot,rad,temp(1,*)*scale,thick=4,linestyle=i
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

 temp=READ_ASCII(filepath('torque_time.dat',root_dir='.',subdir=[strcompress(cases(0),/remove_all)]))
 temp=double(temp.(0))

 set_plot, 'ps'
 device, filename='torque_plot2.ps' $
   ,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches

 plot,temp(0,*),temp(3,*)*1d4,xmargin=[8,2],ymargin=[3,2] $
   ,ytitle=textoidl('Torque'),xtitle='t/orbits' $
   ,charsize=1.5,xminor=4, thick=4, xrange=xrange,xtickinterval=xtickinterval $
   ,yrange=yrange

 for i=1, numcases-1 do begin
;     if keyword_set(scale) then temp(2,*) *= scale(i)

     ks = string(start(i),format='(I03)') 
     temp=READ_ASCII(filepath('torque_time.dat',root_dir='.',subdir=[strcompress(cases(i),/remove_all)]))
     temp=double(temp.(0))
    
     oplot,temp(0,*),temp(3,*)*1d4,thick=4,linestyle=i
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
