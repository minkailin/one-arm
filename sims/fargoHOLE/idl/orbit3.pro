PRO orbit3, cases=cases, ntime=ntime, legend=legend, label=label, xrange=xrange, yrange=yrange $
            ,ytickinterval=ytickinterval, xtickinterval=xtickinterval, mp=mp

bigM = mp + 1d0
cases=string(cases)
numcases=n_elements(cases)
orbs=dblarr(8,ntime,numcases)
for i=0, numcases-1 do begin
    temp=READ_ASCII(filepath('bigplanet0.dat',root_dir='.',subdir=[strcompress(cases(i),/remove_all)]))
    temp=double(temp.(0))
    a=n_elements(temp(0,*))/ntime
    a=fix(a)
     for j=0L, ntime-1 do begin
     orbs(0,j,i) = temp(7,j*a);time
     orbs(1,j,i) = temp(1,j*a);x 
     orbs(2,j,i) = temp(2,j*a);y     
     orbs(3,j,i) = sqrt(orbs(1,j,i)^2.0+orbs(2,j,i)^2.0);rp(t)
     orbs(4,j,i) = temp(6,j*a)
     orbs(0,j,i)/= 2.0*!dpi*orbs(3,0,i)^1.5

     orbs(5,j,i) = temp(3,j*a);vx 
     orbs(6,j,i) = temp(4,j*a);vy  
     vsq = orbs(5,j,i)^2 + orbs(6,j,i)^2
     orbs(7,j,i) = bigM*orbs(3,j,i)/(2d0*bigM - orbs(3,j,i)*vsq) 

;h3 = x*vy - y*vx
;angmom = sqrt(h1^2 + h2^2 + h3^2)


;incl      = atan(sqrt(h1^2 +h2^2)/h3)
;semimajor = bigM*r/(2d0*bigM - r*vel^2)
;eccen     = sqrt(1d0 - angmom^2/(bigM*semimajor))

; semimajor = bigM*r/(2d0*bigM - r*vel^2)

 endfor
endfor



;;;;;;;;;;
;PLOTTING;
;;;;;;;;;;
set_plot, 'ps'
device, filename='orbit3_rp.ps', xsize=8, ysize=4.5, xoffset=0, yoffset=0, /inches,/color,bits_per_pixel=8
plot, orbs(0,*,0), orbs(3,*,0), xmargin=[6,2], ymargin=[3,2], xtitle='t/orbits' $
, ytitle=textoidl('r_p(t)'), yrange=yrange, ytickinterval=ytickinterval $
, charsize=1.5,yminor=10,thick=4,linestyle=0, xrange=xrange, xtickinterval=xtickinterval, xminor=5
for k=1, numcases-1 do begin
oplot, orbs(0,*,k), orbs(3,*,k),thick=4,linestyle=k
endfor

if keyword_set(legend) then begin
x0=legend(0)
x1=legend(1)
y1=legend(2)
dy=legend(3)
for k=0, numcases-1 do begin
oplot,[x0,x1],[y1,y1]-k*dy,thick=4,linestyle=k
xyouts,x1,y1-k*dy,textoidl(label(k)),charsize=1.5
endfor
endif

set_plot, 'ps'
device, filename='orbit3_a.ps', xsize=8, ysize=4.5, xoffset=0, yoffset=0, /inches,/color,bits_per_pixel=8
plot, orbs(0,*,0), orbs(7,*,0), xmargin=[6,2], ymargin=[3,2], xtitle='t/orbits' $
, ytitle=textoidl('a(t)'), yrange=yrange, ytickinterval=ytickinterval $
, charsize=1.5,yminor=10,thick=4,linestyle=0, xrange=xrange, xtickinterval=xtickinterval, xminor=5
for k=1, numcases-1 do begin
oplot, orbs(0,*,k), orbs(3,*,k),thick=4,linestyle=k
endfor

if keyword_set(legend) then begin
x0=legend(0)
x1=legend(1)
y1=legend(2)
dy=legend(3)
for k=0, numcases-1 do begin
oplot,[x0,x1],[y1,y1]-k*dy,thick=4,linestyle=k
xyouts,x1,y1-k*dy,textoidl(label(k)),charsize=1.5
endfor
endif




; set_plot, 'ps'
; device, filename='orbit3_massloss.ps', xsize=8, ysize=4.5, xoffset=0, yoffset=0, /inches,/color,bits_per_pixel=8
; plot, orbs(0,*,0), orbs(4,*,0)*1d6, xmargin=[6,2], ymargin=[3,2], xtitle='t/orbits' $
; , ytitle=textoidl('M_{loss}\times10^6'), yrange=yrange, ytickinterval=ytickinterval $
; , charsize=1.5,yminor=5,thick=4,linestyle=0, xrange=xrange, xtickinterval=xtickinterval, xminor=5
; for k=1, numcases-1 do begin
; oplot, orbs(0,*,k), orbs(4,*,k)*1d6,thick=4,linestyle=k
; endfor

; if keyword_set(legend) then begin
; x0=legend(0)
; x1=legend(1)
; y1=legend(2)
; dy=legend(3)
; for k=0, numcases-1 do begin
; oplot,[x0,x1],[y1,y1]-k*dy,thick=4,linestyle=k
; xyouts,x1,y1-k*dy,label(k),charsize=1.5
; endfor
; endif

 device,/close

END
