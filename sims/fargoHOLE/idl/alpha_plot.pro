pro alpha_plot, cases=cases, legend=legend, label=label, xtickinterval=xtickinterval, start=start, xrange=xrange $
,yrange=yrange,ytickinterval=ytickinterval

cases=string(cases)
numcases=n_elements(cases)

filename = strcompress('alpha_time.dat',/remove_all)
temp=READ_ASCII(filepath(filename,root_dir='.',subdir=[strcompress(cases(0),/remove_all)]))
temp=double(temp.(0))

nel = n_elements(temp(0,*))
a = min(abs(temp(0,*)-25.0), beg)
avg_alpha = dblarr(nel-beg)
for j = 0, nel-beg-1 do avg_alpha(j) = mean(temp(2,beg: beg+j))

set_plot, 'ps'
device, filename='alpha_plot.ps' $
  ,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches

plot,temp(0,*), temp(1,*)*1d3, xmargin=[6,2],ymargin=[3,2] $
  ,ytitle=textoidl('<\alpha>\times10^{3}'),xtitle='t/orbits' $
  ,charsize=1.5, thick=4,xtickinterval=xtickinterval,ytickinterval=ytickinterval $
  ,xrange=xrange,yrange=yrange; , /nodata,ystyle=4;, title=textoidl('Q_m=3')

;axis, /save, yaxis=0, ytitle=textoidl('<\alpha>\times10^{3}, solid'),charsize=1.5, yrange=[1,10]
;oplot,temp(0,*), temp(2,*)*1d3 , linestyle=0, thick=4

;axis, /save, yaxis=1, ytitle=textoidl('<\Delta\Sigma/\Sigma>_+, dotted'),charsize=1.5, yrange=[0.1,0.8] $
;  , ytickinterval=0.1
;oplot,temp(0,*), temp(1,*), linestyle=1,thick=4

for i=1, numcases-1 do begin
    temp=READ_ASCII(filepath(filename,root_dir='.',subdir=[strcompress(cases(i),/remove_all)]))
    temp=double(temp.(0))
    nel = n_elements(temp(0,*))
    a = min(abs(temp(0,*)-25.0), beg)
    avg_alpha = dblarr(nel-beg)
    for j=0, nel-beg-1 do avg_alpha(j) = mean(temp(2,beg: beg+j))
    
    oplot,temp(0,*), temp(1,*)*1d3, thick=4,linestyle=i
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


 numcases = n_elements(start)
 ks=string(start,format='(I03)')
 filename = strcompress('alpha_'+ks(0)+'.dat',/remove_all)
 temp=READ_ASCII(filepath(filename,root_dir='.',subdir=[strcompress(cases(0),/remove_all)]))
 temp=double(temp.(0))

 set_plot, 'ps'
 device, filename='alpha_plot2.ps' $
 ,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches

 plot,temp(0,*), temp(1,*)*1d3, xmargin=[6,2],ymargin=[3,2] $
 ,ytitle=textoidl('<\alpha>_\phi\times10^{3}'),xtitle='r' $
 ,charsize=1.5, thick=4,xtickinterval=xtickinterval,ytickinterval=ytickinterval $
   ,xrange=xrange, yrange=yrange;, title=time+' orbits'

 for i=1, numcases-1 do begin
     filename = strcompress('alpha_'+ks(i)+'.dat',/remove_all)
     temp=READ_ASCII(filepath(filename,root_dir='.',subdir=[strcompress(cases(i),/remove_all)]))
     temp=double(temp.(0))    
     oplot,temp(0,*), temp(1,*)*1d3, thick=4,linestyle=i
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


;filename2 = strcompress('shear_'+ks+'.dat',/remove_all)
;temp2=READ_ASCII(filepath(filename2,root_dir='.',subdir=[strcompress(cases(0),/remove_all)]))
;temp2=double(temp2.(0))

;sg_shear = abs(temp2(1,*))

;temp2(1,*)/=-1.5*temp2(0,*)^(-2.5)

;temp2(1,*) = deriv(temp2(0,*),temp2(0,*)^3d0*temp2(1,*))

; set_plot, 'ps'
; device, filename='shear_plot.ps' $
; ,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches

; ;ytitle = textoidl("<\Sigma\Omega'>_\phi/(\Sigma_{t=0}\Omega_k')")
; ytitle = textoidl("(r^3\Sigma\Omega^{'})^{'}\times10^3")

; plot,temp2(0,*), temp2(1,*)*1d3,xmargin=[6,2],ymargin=[3,2] $
; ,ytitle=ytitle,xtitle='r' $
; ,charsize=1.5, thick=4,xtickinterval=xtickinterval,ytickinterval=ytickinterval,title=time+' orbits' $
; ,xrange=xrange,yrange=yrange


; ; temp = min(abs(temp2(0,*)-xrange(0)),x0)
; ; temp = min(abs(temp2(0,*)-xrange(1)),x1)

; ; print, 'average viscous torque is', mean(temp2(1,x0:x1))

; for i=1, numcases-1 do begin
;     temp2=READ_ASCII(filepath(filename2,root_dir='.',subdir=[strcompress(cases(i),/remove_all)]))
;     temp2=double(temp2.(0))
;     ;temp2(1,*)/=-1.5*temp2(0,*)^(-2.5)

;     temp2(1,*) = deriv(temp2(0,*),temp2(0,*)^3d0*temp2(1,*))

;     oplot,temp2(0,*),temp2(1,*)*1d3,thick=4,linestyle=i

;     ; temp = min(abs(temp2(0,*)-xrange(0)),x0)
; ;     temp = min(abs(temp2(0,*)-xrange(1)),x1)
    
; ;     print, 'average viscous torque is', mean(temp2(1,x0:x1))

; endfor

; if keyword_set(legend) then begin
; x0=legend(0)
; x1=legend(1)
; y0=legend(2)
; dy=legend(3)
; for j=0, numcases-1 do begin
;     oplot, [x0,x1], [y0,y0]-dy*j, thick=4, linestyle=j
;     xyouts, x1, y0-dy*j,label(j),charsize=1.5
; endfor
; endif

; device,/close

end



;dims=(read_ascii(filepath('dims.dat',root_dir='.',subdir=[cases(0)]))).(0)
;info = dblarr(11,start+1)
;openr,3,filepath('planet0.dat',root_dir='.',subdir=[cases(0)])
;readf,3,info
;close,3
;a0=info(1,0)
;p0=2.*!dpi*(a0)^(3./2.)
;time=string(info(7,start)/p0,format='(F7.2)')

;ks=string(start,format='(I03)')
