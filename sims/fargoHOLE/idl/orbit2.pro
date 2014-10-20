PRO orbit2, cases=cases,  ntime=ntime, smooth=smooth, segment=segment $
, adot=adot, a0=a0,ecc=ecc,xrange=xrange,yrange=yrange, val=val, order=order, legend=legend $
, label=label, xtickinterval=xtickinterval
;a0 can manually set initial semimajor axis to get orignal orbital
;time (standard unit of time). required if simulation doesn't run from t=0.
pi=3.141592654
;;;;;;;;;;;;;;;;;
;SURFACE DENSITY;
;;;;;;;;;;;;;;;;;
;sigma=dblarr(finish-start+1)
;openr,lun,'sigmas.dat', /get_lun
;readf, lun, sigma, format='(d5.1)'
;close,lun
;;;;;;;;;;;;;;;;;;
;SET UP TIME AXIS;
;;;;;;;;;;;;;;;;;;
;data=filepath('dims.dat',root_dir='.',subdir=['out1'])
;dims=(read_ascii(data)).(0)
;time=findgen(dims(0,6))
;time=(time+1)*(2.*3.141592654)
;ntime=n_elements(time)
;;;;;;;;;;;;;;;;;
;GET PLANET DATA;
;;;;;;;;;;;;;;;;;
cases=string(cases)
numcases=n_elements(cases)
orbs=dblarr(3,ntime,numcases)
for i=0, numcases-1 do begin
    temp=READ_ASCII(filepath('orbit0.dat',root_dir='.',subdir=[strcompress(cases(i),/remove_all)]))
    temp=double(temp.(0))
    a=n_elements(temp(0,*))/ntime
    a=fix(a)
     for j=0L, ntime-1 do begin
     orbs(0,j,i)=temp(0,j*a);time
     orbs(1,j,i)=temp(2,j*a);semimajor axis
     orbs(2,j,i)=temp(1,j*a);eccentricity
 endfor
endfor
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;REAL TIME TO ORBITS CONVERSION;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
if not keyword_set(a0) then a0=orbs(1,0,0)
p0=2.*pi*a0^(3./2.)
orbs(0,*,*)=orbs(0,*,*)/p0
total_orb=orbs(0,ntime-1,0)
;Polynomial curve fitting. Fit a polynomial to a segment of a
;curve. 'segment' is a array [a,b] indicating start and end of curve
;in orbits.
if keyword_set(segment) then begin
for k=0, numcases-1 do begin
temp1=min(abs(orbs(0,*,k)-segment(0)),seg1)
temp2=min(abs(orbs(0,*,k)-segment(1)),seg2)
xx=orbs(0,seg1:seg2,k)
yy=orbs(1,seg1:seg2,k)
result=poly_fit(xx,yy,order,/double,sigma=sigma)
xplot=segment(0)+(segment(1)-segment(0))*dindgen(100)/99.
yplot=dblarr(100)
for i=0, 99 do begin 
yplot(i)=0.0
for m=0, order do yplot(i)+=result(m)*xplot(i)^double(m)
endfor
print, strcompress(string(order,format='(I01)'))+'th order polynomial for curve '+string(k)
print, result
;print, 'errors are ', sigma
if keyword_set(adot) then begin
rate=0.0
for i=1, order do rate+=double(i)*result(i)*double(adot)^double(i-1)
print, 'migration rate at '+strcompress(string(adot,format='(F5.2)'))+' orbits is ='+string(rate,format='(E9.2)')
endif
endfor
endif
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;SMOOTHING OPTION, smooth is number of orbits to average over;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
if keyword_set(smooth) then begin
ncell=smooth/(total_orb/(ntime-1))
ntime1=total_orb/smooth
ncell=round(ncell)
ntime1=floor(ntime1)
orb_smooth=dblarr(2,ntime1,numcases)
for m=0, numcases-1 do begin
for n=0, ntime1-2 do begin
orb_smooth(0,n,m)=mean(orbs(0,n*ncell:ncell*(n+1),m))
orb_smooth(1,n,m)=mean(orbs(1,n*ncell:ncell*(n+1),m))
endfor
orb_smooth(0,ntime1-1,m)=mean(orbs(0,ncell*(ntime1-1):ntime-1,m))
orb_smooth(1,ntime1-1,m)=mean(orbs(1,ncell*(ntime1-1):ntime-1,m))
endfor
endif
;;;;;;;;;;
;PLOTTING;
;;;;;;;;;;
if not keyword_set(xrange) then xrange=[0.,175.]
if not keyword_set(yrange) then yrange=[min(orbs(1,*,*)),max(orbs(1,*,*))]
set_plot, 'ps'
device, filename='orbit2_a.ps', xsize=8, ysize=4.5, xoffset=0, yoffset=0, /inches,/color,bits_per_pixel=8
if keyword_set(smooth) then begin
loadct,1
plot, orb_smooth(0,*,0), orb_smooth(1,*,0), xmargin=[14,6], ymargin=[6,6], xtitle='Orbits', ytitle='a', yrange=[1.5, 2.0],/nodata $
,xtickinterval=50, xminor=10
for k=0, numcases-1 do begin
oplot, orb_smooth(0,*,k), orb_smooth(1,*,k),thick=5,color=128
endfor
endif else begin
;val=strcompress(string(val,format='(f4.2)'),/remove_all)
plot, orbs(0,*,0), orbs(1,*,0), xmargin=[8,2], ymargin=[3,1], xtitle='t/orbits', ytitle='a', yrange=yrange $
,charsize=1.5,yminor=5,thick=4,linestyle=0,xrange=xrange,xtickinterval=xtickinterval
;oplot,[100,125],[1.9,1.9],thick=4,linestyle=0
;xyouts,125,1.9,textoidl(' \nu_0=')+val(0),charsize=1.5
linestyle = 0
for k=1, numcases-1 do begin        
    if ((k mod 2) ne 0) then begin
        thick=2
        linestyle = linestyle
    endif else begin
        thick = 4
        linestyle =  linestyle+1        
    endelse

    print, k, linestyle
    oplot, orbs(0,*,k), orbs(1,*,k),thick=thick,linestyle=linestyle

;oplot,[100,125],[1.9-k*0.1,1.9-k*0.1],thick=4,linestyle=k
;xyouts,125,1.9-0.1*k,textoidl(' \nu_0=')+val(k),charsize=1.5
endfor

if keyword_set(legend) then begin
    x0=legend(0)
    x1=legend(1)
    y0=legend(2)
    dy=legend(3)
    for j=0, n_elements(label)-1 do begin
        ;if ((j mod 2) ne 0) then begin
        ;    thick=2
        ;    linestyle = j - 1 
        ;endif else begin
        ;    thick = 4
        ;    linestyle =  j
        ;endelse
        oplot, [x0,x1], [y0,y0]-dy*j, thick=4, linestyle=j
        xyouts, x1, y0-dy*j,label(j),charsize=1.5
    endfor
endif

endelse
device,/close
if keyword_set(ecc) then begin
set_plot, 'ps'
device, filename='orbit2_ecc.ps', xsize=8, ysize=4.5, xoffset=0, yoffset=0, /inches,/color,bits_per_pixel=8
plot, orbs(0,*,0), orbs(2,*,0), xmargin=[8,2], ymargin=[3,2], xtitle='Orbits', ytitle='e', yrange=[0.0, max(orbs(2,*,*))]$
,charsize=1.5,thick=4,linestyle=0,ystyle=2 
; oplot,[120,150],[2,2],thick=4,linestyle=0
; xyouts,150,2,textoidl('\nu_0=')+strcompress(string(val(0)),/remove_all)
; for k=1, numcases-1 do begin
; oplot, orbs(0,*,k), orbs(2,*,k),thick=4,linestyle=k
; oplot,[120,150],[2-k*0.1,2-k*0.1],thick=4,linestyle=k
; xyouts,150,2,textoidl('\nu_0=')+strcompress(string(val(k)),/remove_all)
; endfor
device,/close
endif

if(numcases eq 1) then begin
set_plot, 'ps'
device, filename='orbit2_.ps', xsize=8, ysize=4.5, xoffset=0, yoffset=0, /inches,/color,bits_per_pixel=8
plot, orbs(0,*,0), orbs(2,*,0), xmargin=[8,8], ymargin=[3,2], xtitle='Orbits', ytitle='e',/nodata $
,charsize=1.5,thick=4,linestyle=0,ystyle=4,xrange=xrange
axis,yaxis=0,charsize=1.5,ytitle='a, solid',/save,yrange=[min(orbs(1,*,0)),max(orbs(1,*,0))]
oplot, orbs(0,*,0), orbs(1,*,0), linestyle=0, thick=4
axis,yaxis=1,charsize=1.5,ytitle='e, dashed',/save,yrange=[min(orbs(2,*,0)),max(orbs(2,*,0))]
oplot, orbs(0,*,0), orbs(2,*,0), linestyle=2, thick=4
device,/close
endif

END
