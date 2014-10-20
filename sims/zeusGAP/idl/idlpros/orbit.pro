PRO orbit, start=start, finish=finish, ntime=ntime, smooth=smooth, segment=segment, torb=torb
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
orbs=dblarr(2,ntime,finish-start+1)
for i=start, finish do begin
    num=string(i)
    temp=READ_ASCII(filepath('orbit0.dat',root_dir='.',subdir=[strcompress('out'+num,/remove_all)]))
    temp=temp.(0)
    a=n_elements(temp(0,*))/ntime
    a=fix(a)
     for j=0, ntime-1 do begin
     orbs(0,j,i-start)=temp(0,j*a)
      orbs(1,j,i-start)=temp(2,j*a)
 endfor
endfor
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;REAL TIME TO ORBITS CONVERSION;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
a0=orbs(1,0,0)
p0=2.*pi*a0^(3./2.)
orbs(0,*,*)=orbs(0,*,*)/p0
total_orb=orbs(0,ntime-1,0)

;Polynomial curve fitting. Fit a quadratic curve to a segment of a
;curve. 'segment' is a array [a,b] indicating start and end of curve
;in orbits.
if keyword_set(segment) then begin
for k=start, finish do begin
temp1=abs(orbs(0,*,k-start)-segment(0))
temp2=abs(orbs(0,*,k-start)-segment(1))
seg1=where(temp1 eq min(temp1))
seg2=where(temp2 eq min(temp2))
xx=orbs(0,seg1:seg2,k-start)*p0
yy=orbs(1,seg1:seg2,k-start)
result=poly_fit(xx,yy,2,/double,sigma=sigma)
print, 'quadratic fit to curve '+string(k)+' is'
print, 'a(t)='+string(result(0))+string(result(1))+'t'+string(result(2))+'t^2'
print, 'errors are ', sigma
if keyword_set(torb) then begin
treal=double(p0*torb)
rate=result(1)+2.*result(2)*treal
print, 'migration rate at '+strcompress(string(torb,format='(F5.2)'))+' orbits is ='+string(rate,format='(E9.2)')
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
orb_smooth=dblarr(2,ntime1,finish-start+1)
for m=0, finish-start do begin
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
set_plot, 'ps'
device, filename='orbit.ps', xsize=8, ysize=6, xoffset=0, yoffset=0, /inches,/color,bits_per_pixel=8
if keyword_set(smooth) then begin
loadct,1
plot, orb_smooth(0,*,0), orb_smooth(1,*,0), xmargin=[14,6], ymargin=[6,6], xtitle='Orbits', ytitle='a', yrange=[1.5, 2.0],/nodata $
,xrange=[0,175],xtickinterval=25, xminor=5
for k=0, finish-start do begin
oplot, orb_smooth(0,*,k), orb_smooth(1,*,k),thick=5,color=128
endfor
endif else begin
;pos=320
plot, orbs(0,*,0), orbs(1,*,0), xmargin=[6,4], ymargin=[3,3], xtitle='Orbits', ytitle='a', yrange=[0.4, 2.4]$
,charsize=1.5,ytickinterval=0.4,yminor=2,xstyle=2,thick=2
;xyouts, orbs(0,pos,1), orbs(1,pos,0),textoidl('\nu_0=0'),charsize=1.5
for k=0, finish-start do begin
oplot, orbs(0,*,k), orbs(1,*,k),thick=2
endfor
; xyouts,150.,1., textoidl('\nu_0=0.1'),charsize=1.5
; arrow,orbs(0,400,1)+1,orbs(1,400,1),150.,1.,/data, hsize=0
; xyouts,orbs(0,ntime-1,2),orbs(1,ntime-1,2), textoidl('\nu_0=1.0'),charsize=1.5
; xyouts,orbs(0,ntime-1,3),orbs(1,ntime-1,3), textoidl('\nu_0=10'),charsize=1.5
; xyouts,orbs(0,ntime-1,4),orbs(1,ntime-1,4), textoidl('\nu_0=100'),charsize=1.5
; xyouts,orbs(0,ntime-1,5),orbs(1,ntime-1,5), textoidl('\nu_0=5.0'),charsize=1.5
; xyouts,orbs(0,ntime-1,6)+5,orbs(1,ntime-1,6)-0.1, textoidl('\nu_0=50'),charsize=1.5
; arrow,orbs(0,ntime-1,6)+1,orbs(1,ntime-1,6),orbs(0,ntime-1,6)+5,orbs(1,ntime-1,6)-0.1,/data, hsize=0
endelse
device,/close
END
