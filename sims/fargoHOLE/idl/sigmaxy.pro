; THIS PROCEDURE DOES A CARTESIAN CONTOUR PLOT OF DENSITY (ANNULUS
; ABOUT PLANET)
; 20/12/2008: Added option to plot tracer fluid. Set keyword tracer to enable. 
pro sigmaxy, loc=loc, start=start, finish=finish, xs=xs, tracer=tracer $
,ct=ct,plotrange=plotrange, xrange=xrange, yrange=yrange,jump=jump,mass=mass $
, order=order
common consts, pi
common curve, coeff
pi=3.141592654
factor0=double((mass/3.)^(1./3.))
if n_elements(start) ne 1 then begin
cases=start
endif
if not keyword_set(finish) then begin
cases=[start]
endif else cases=start+indgen(finish-start+1)
numcases=n_elements(cases)
if not keyword_set(ct) then ct=5
if not keyword_set(xs) then xs=2.0
if not keyword_set(xrange) then xrange=[-xs,xs]
if not keyword_set(yrange) then yrange=[-1.0,1.0]
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;DATA LOCATION,GET RUN INFO, READ PLANET ORBIT INFO.;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
location=strcompress(loc,/remove_all)
dims=(read_ascii(filepath('dims.dat',root_dir='.',subdir=[location]))).(0)
nout=fix(dims(5))
nrad=fix(dims(6))
nsec=fix(dims(7))
info=dblarr(11,nout+1)
openr,3,filepath('planet0.dat',root_dir='.',subdir=[location])
readf,3,info
close,3
a0=info(1,0)
dt=info(7,1)
p0=2.*pi*(a0)^(3./2.)
;Create 1d arrays for radial and azimuthal values for contour
;plot,
azifix=dblarr(nsec)
radtmp=dblarr(nrad+1)
openr,1,filepath('used_rad.dat',root_dir='.',subdir=[location])
readf,1,radtmp
close,1
azifix=dindgen(nsec)*2.*pi/nsec
;Arrays to hold variables.Also create 2d array to hold radii values
;over grid.
sigma=dblarr(nsec,nrad)
label=dblarr(nsec,nrad)
rmed=radtmp(0:nrad-1)+radtmp(1:nrad)
rmed=rmed*0.5
;Cell sizes.
dtheta=2.*pi/nsec
dr=(radtmp(nrad)-radtmp(0))/nrad
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;SETUP DONE, START CALCULATION AND PLOTTING;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
xxtitle=textoidl('(r-r_p)/r_h')
yytitle=textoidl('(\phi-\phi_p)/\pi')
loadct,ct, bottom=0
for k=0, numcases-1 do begin
ks=string(cases(k),format='(I03)')
;Read raw data.
openr,2,filepath(strcompress('gasdens'+string(cases(k))+'.dat',/remove_all),root_dir='.',subdir=[location]) 
readu,2,sigma
close,2
if keyword_set(tracer) then begin
openr,3,filepath(strcompress('gaslabel'+string(cases(k))+'.dat',/remove_all),root_dir='.',subdir=[location])
readu,3,label
close,3
endif
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;NOW DO THE POLAR CONTOUR PLOTS;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
if keyword_set(tracer) then begin
data=label(0:nsec-1,1:nrad-2)
endif else data=sigma(0:nsec-1,1:nrad-2)
;Physical size of half width. Find planet radius and azimuth, and find
;range of grid points that correspond to sector
plx=info(1,cases(k))
ply=info(2,cases(k))
plrad=sqrt(plx*plx+ply*ply)
rhill=plrad*factor0
phi=pltphi(plx,ply)
rplot=(rmed-plrad)/rhill
arr1=min(abs(rplot-xrange(0)),rmns)
arr2=min(abs(rplot-xrange(1)),rpls)
azi=azifix-phi
data2=transpose(data(0:nsec-1,rmns:rpls))
; if not keyword_set(plotrange) then plotrange=[alog(min(data2(where(data2 gt 0.)))),alog(max(data2))]
; if not keyword_set(tracer) then data2=alog(data2)
plotrange=[min(data2),max(data2)]
levels=plotrange(0)+(plotrange(1)-plotrange(0))*(dindgen(128)/128.)
time=string(info(7,cases(k))/p0,format='(F6.2)')
data3=dblarr(rpls-rmns+1,nsec)
azinew=dblarr(nsec)
if phi gt pi then begin
excess=where(azi lt -pi)
azi(excess)=2.*pi+azi(excess)
seg1=where(azi eq min(azi))
if seg1 ne nsec-1 then begin
data3(*,0:nsec-seg1-1)=data2(*,seg1:nsec-1)
azinew(0:nsec-seg1-1)=azi(seg1:nsec-1)/pi
endif
if seg1 ne 1 then begin
data3(*,nsec-seg1:nsec-1)=data2(*,0:seg1-1)
azinew(nsec-seg1:nsec-1)=azi(0:seg1-1)/pi
endif
endif else if (0. le phi and  phi lt pi and (pi-phi) gt dtheta ) then begin
excess=where(azi ge pi)
azi(excess)=azi(excess)-2.*pi
seg2=where(azi eq max(azi))
if seg2 ne nsec-2 then begin
data3(*,0:nsec-seg2-2)=data2(*,seg2+1:nsec-1)
azinew(0:nsec-seg2-2)=azi(seg2+1:nsec-1)/pi
endif
if seg2 ne 0 then begin
data3(*,nsec-seg2-1:nsec-1)=data2(*,0:seg2)
azinew(nsec-seg2-1:nsec-1)=azi(0:seg2)/pi
endif
endif else begin
data3=data2
azinew(0:nsec-1)=azi(0:nsec-1)/pi
endelse
;The actual plotting part.
set_plot, 'ps'
device, filename=filepath(strcompress('sigmaxy_'+ks+'.ps',/remove_all),root_dir='.',subdir=[location])$
,/color, bits_per_pixel=8,xsize=12, ysize=12
temp=strcompress(string(plrad,format='(f5.2)'),/remove_all)
contour,data3,rplot(rmns:rpls),azinew(0:nsec-1) ,/fill $
,levels=levels,title=time+' orbits,'+textoidl(' r_p=')+temp,ymargin=[4,2],xmargin=[8,10],color=0$
,xtitle=xxtitle,ytitle=yytitle,charsize=1.5,yrange=yrange,xrange=xrange,xtickinterval=2.0
colorbar, position=[0.95, 0.18, 1, 0.88],/vertical,/left,range=plotrange,format='(e10.1)',color=0,charsize=1.5
oplot,[0.,0.],[0.,0.],psym=6,symsize=1,color=-1
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;Option to indicate shock location and pre-post sampling points;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
if keyword_set(jump) then begin
if n_elements(jump) eq 1 then jump=[jump,jump]
shock_fit,loc=loc,start=cases(k),xrange=[0.0,xrange(1)] ,yrange=yrange $
, eta=0.9, order=order,/noplot,mass=mass
ord=n_elements(coeff)
;Axes and arrays
pts=floor(xrange(1)*rhill/dr)
xshock=rhill*xrange(1)*dindgen(pts)/(pts-1.0)
y_shock=dblarr(pts)
pre_shock=dblarr(2,pts)
post_shock=dblarr(2,pts)
;Shock curve y=y(x) and it's derivative
for i=0, pts-1 do begin
dys=0.0
ys=0.0
for j=0, ord-1 do ys+=coeff(j)*xshock(i)^double(j)
for j=1, ord-1 do dys+=double(j)*coeff(j)*xshock(i)^double(j-1)
x0=xshock(i)
theta=atan(-1.0/dys)
dx=jump*rhill*cos(theta)
dy=jump*rhill*sin(theta)
x1=x0+dx(0);pre-shock sample x-coordinate
y1=ys+dy(0)
x2=x0-dx(1);post-shock sample x-coordinate
y2=ys-dy(1)
y_shock(i)=ys
pre_shock(0,i)=x1
pre_shock(1,i)=y1
post_shock(0,i)=x2
post_shock(1,i)=y2
;print, ys, dy(1), y2, jump(1), rhill, sin(theta)
endfor
xshock/=rhill
y_shock/=plrad*pi
pre_shock(0,*)/=rhill
pre_shock(1,*)/=plrad*pi
post_shock(0,*)/=rhill
post_shock(1,*)/=plrad*pi
oplot,xshock,y_shock,color=-1,thick=2
oplot,pre_shock(0,*),pre_shock(1,*),color=-1,linestyle=2,thick=2
oplot,post_shock(0,*),post_shock(1,*),linestyle=2,thick=2
endif
device,/close
print, 'done'+ks
endfor
end




































;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;OLD CODE FOR CONTOUR PLOT (CHANGING PLOTTING ORDER)
; if seg1 ne nsec-1 then begin
; contour,data2(*,seg1:nsec-1),rplot(rmns:rpls),azi(seg1:nsec-1)/pi ,/fill $
; ,levels=levels,title=time+' orbits',ymargin=[4,2],xmargin=[8,10],color=0$
; ,xtitle=xxtitle,ytitle=yytitle,charsize=1.5,yrange=yrange,xrange=xrange
; endif
; if seg1 ne 1 then begin
; contour,data2(*,0:seg1-1),rplot(rmns:rpls),azi(0:seg1-1)/pi ,/fill $
; ,levels=levels,/overplot,yrange=yrange,title=time+' orbits',ymargin=[4,2],xmargin=[8,10],color=0$
; ,xtitle=xxtitle,ytitle=yytitle,charsize=1.5,xrange=xrange
; endif
; if seg2 ne nsec-2 then begin
; contour,data2(*,seg2+1:nsec-1),rplot(rmns:rpls),azi(seg2+1:nsec-1)/pi ,/fill $
; ,levels=levels,title=time+' orbits',ymargin=[4,2],xmargin=[8,10],color=0$
; ,xtitle=xxtitle,ytitle=yytitle,charsize=1.5,yrange=yrange,xrange=xrange
; endif
; if seg2 ne 0 then begin
; contour,data2(*,0:seg2),rplot(rmns:rpls),azi(0:seg2)/pi ,/fill $
; ,levels=levels,/overplot,yrange=yrange,title=time+' orbits',ymargin=[4,2],xmargin=[8,10],color=0$
; ,xtitle=xxtitle,ytitle=yytitle,charsize=1.5,xrange=xrange
; endif
;endelse
;do not plot points in order of index. may be gap so manually insert
;line segment
;arrow, azi(0),data1d(0),azi(nsec-1),data1d(nsec-1),hsize=0,/data
;data3=dblarr(rpls-rmns+1,2)
;data3(*,0)=data2(*,nsec-1)
;data3(*,1)=data2(*,0)
;if keyword_set(upper) then begin
;contour,data3(*,0:1),rplot(rmns:rpls),[azi(nsec-1),azi(0)]/pi,/fill,levels=levels,/overplot
;endif else begin
;contour,data3(*,0:1),rplot(rmns:rpls),[azi(0),azi(0)]/pi,/fill,levels=levels,/overplot
;endelse
