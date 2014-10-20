;THIS PROCEDURE CALCULATES AND PLOTS THE PLANETS' MIGRATION RATE,
;da(t)/dt AS A FUNCTION OF TIME. RAW DATA IS A HUGE ARRAY, NTIME IS
;THE DESIRED NUMBER OF SAMPLES TO TAKE FROM THIS ARRAY.

PRO adot, start=start, finish=finish, ntime=ntime, smooth=smooth
pi=3.141592654
Dt=(pi/10.)*2.0
;;;;;;;;;;;;;;;;;
;GET PLANET DATA;
;;;;;;;;;;;;;;;;;
orbs=dblarr(2,ntime,finish-start+1)
num=file_lines(filepath('orbit0.dat',root_dir='.',subdir=[strcompress('out14',/remove_all)]))/ntime
num=fix(num)
adot=dblarr(2,num)
for i=start, finish do begin
    name=string(i)
    data=READ_ASCII(filepath('orbit0.dat',root_dir='.',subdir=[strcompress('out'+name,/remove_all)]))
    data=data.(0)
     for j=0, ntime-1 do begin
         for k=0,num-1 do begin
             adot(0,k)=data(0,j*num+k+1)
             adot(1,k)=(data(2,j*num+k+2)-data(2,j*num+k))/Dt
         endfor
         orbs(0,j,i-start)=mean(adot(0,*))
         orbs(1,j,i-start)=mean(adot(1,*))
     endfor
 endfor
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;REAL TIME TO ORBITS CONVERSION (should really take data from planet0.dat;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
a0=data(2,0)
p0=2.*pi*a0^(3./2.)
orbs(0,*,*)=orbs(0,*,*)/p0
total_orb=orbs(0,ntime-1,0)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;SMOOTHING OPTION, smooth is number of orbits to average over;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
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
device, filename='adot.ps', xsize=8, ysize=6, xoffset=0, yoffset=0, /inches,/color,bits_per_pixel=8
if keyword_set(smooth) then begin
loadct,1
orb_smooth(1,*,*)=orb_smooth(1,*,*)*10000.
plot, orb_smooth(0,*,0), orb_smooth(1,*,0), xmargin=[8,2], ymargin=[3,3], xtitle='Orbits' $
, ytitle='(da/dt)/'+textoidl('10^{-4}'), yrange=[min(orb_smooth(1,*,*)), max(orb_smooth(1,*,*))],/nodata $
,xrange=[0,175],xtickinterval=25, xminor=5,charsize=1.5
for k=0, finish-start do begin
oplot, orb_smooth(0,*,k), orb_smooth(1,*,k),thick=2;,color=128
endfor
endif else begin
;pos=320
orbs(1,*,*)=orbs(1,*,*)*1000.
plot, orbs(0,*,0), orbs(1,*,0), xmargin=[8,2], ymargin=[3,3], xtitle='Orbits', ytitle='(da/dt)/'+textoidl('10^{-4}') $
, yrange=[min(orbs(1,*,*)),max(orbs(1,*,*))], charsize=1.5,xstyle=2,thick=2, xtickinterval=25, xminor=5, xrange=[0.,175.]
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
