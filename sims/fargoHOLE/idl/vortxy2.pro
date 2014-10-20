; THIS PROCEDURE CALCULATES SIGMA/OMEGA ON THE GRID POINTS USING
; CENTRAL DIFFERENCING, AND DOES A CARTESIAN CONTOUR PLOT CENTERED
; ABOUT PLANET.
; 21/12/2008: Added option to overplot tracer.
; 25/2/2009: Option to overplot velocity vectors. WARNING: HYDRO
; VELOCITIES ARE TRANSFORMED TO THOSE CORRESPONDING TO PLOT
; CO-ORDINATES. THE PLOT VECTORS ONLY MAKE SENSE IF YOU TAKE A SMALL
; PATCH (SHEARING BOX). OTHERWISE YOU ARE COMPARING VECTORS AT
; DIFFERENT POINTS. 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;                     SOME FUNCTIONS NEEDED                      ;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;Function decide how to shift azimuthal indices/values
function whichcase, phi
common consts, pi
common grid, nsec, nrad, azi, rmed, dtheta, rmns, rpls, azinew, rplot
common hydro, vtheta, vrad
common planet, info, plx, ply, plrad
common output, seg
if phi gt pi then begin
excess=where(azi lt -pi)
azi(excess)=2.*pi+azi(excess)
seg=where(azi eq min(azi))
option=1
endif else if (0. le phi and  phi lt pi and (pi-phi) gt dtheta ) then begin
excess=where(azi ge pi)
azi(excess)=azi(excess)-2.*pi
seg=where(azi eq max(azi))
option=2
endif else begin
option=3
endelse
return, option
end
;Actual function that does the shift.
function shiftazi, input, opt
common consts, pi
common grid, nsec, nrad, azi, rmed, dtheta, rmns, rpls, azinew, rplot
common hydro, vtheta, vrad
common planet, info, plx, ply, plrad
common output, seg
data3=dblarr(rpls-rmns+1,nsec)
in=dblarr(rpls-rmns+1,nsec)
s=size(input)
if s(0) eq 1 then begin 
for i=0, rpls-rmns do in(i,*)=input(*)
endif else in=input
case opt of
1: begin
 if seg ne nsec-1 then begin
 data3(*,0:nsec-seg-1)=in(*,seg:nsec-1)
 endif
 if seg ne 1 then begin
 data3(*,nsec-seg:nsec-1)=in(*,0:seg-1)
 endif
end
2: begin
 if seg ne nsec-2 then begin
 data3(*,0:nsec-seg-2)=in(*,seg+1:nsec-1)
 endif
 if seg ne 0 then begin
 data3(*,nsec-seg-1:nsec-1)=in(*,0:seg)
endif
end
3: begin
data3=in
end
endcase
return, data3
end
;Conversion of planet vx,vy velocity to v_r, omega_p
function convert, r, x_p, y_p, vx_p, vy_p
sine=double(y_p/r)
cosine=double(x_p/r)
vradial=double(vx_p*cosine+vy_p*sine)
omega_p=double(-vx_p*sine+vy_p*cosine)/r
return, [vradial,omega_p]
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;           SUB-PROGRAM FOR PLOTTING VELOCITY FIELD            ;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro velocity, slice=slice, sample=sample,length=length, yrange=yrange, color=color
common consts, pi
common grid, nsec, nrad, azi, rmed, dtheta, rmns, rpls, azinew, rplot
common hydro, vtheta, vrad
common planet, info, plx, ply, plrad
common misc, opt,factor0
;azinew=azinew-dtheta/2.0
;Vrad and vtheta are not defined at the same point. Do an average to get appropriate value of vrad at each vtheta.
;Put these values in the old array. Outer radius boundary values are
;unchanged (manipulations involving these values will result in wrong
;values). These vrad are at cell edge (rather than cell centre) so
;they are defined -dtheta/2.0 in the azimuth relative to the cell-centrered 
vrad(0,0:nrad-2)=(vrad(0,0:nrad-2)+vrad(nsec-1,0:nrad-2)+vrad(nsec-1,1:nrad-1)+vrad(0,1:nrad-1))/4.0d
vrad(1:nsec-1,0:nrad-2)=(vrad(1:nsec-1,0:nrad-2)+vrad(0:nsec-2,0:nrad-2)+vrad(0:nsec-2,1:nrad-1)+vrad(1:nsec-1,1:nrad-1))/4.0d
;Planet velocity
polar_vel=convert(plrad,plx,ply,info(3,slice),info(4,slice))
;Shift and scale velocity field appropriate to plot units (which is
;relative to planet)
for i=0, nsec-1 do begin
vrad(i,*)=vrad(i,*)-polar_vel(0)*rmed(*)/plrad
vtheta(i,*)=vtheta(i,*)/rmed(*)
endfor
vrad=vrad/(factor0*plrad)
vtheta=(vtheta-polar_vel(1))/pi
;Extract relavent sub-set of data.
vx=transpose(vrad(0:nsec-1,rmns:rpls))
vy=transpose(vtheta(0:nsec-1,rmns:rpls))
;Azimuthal shift
vx_new=shiftazi(vx,opt)
vy_new=shiftazi(vy,opt)
;Sample data
vx_plot=congrid(vx_new,sample,sample)
vy_plot=congrid(vy_new,sample,sample)
ll=min(abs(azinew-yrange(0)),lll)
uu=min(abs(azinew-yrange(1)),uuu)
xaxis=congrid(rplot(rmns:rpls),sample,sample)
yaxis=congrid(azinew(lll:uuu),sample,sample)
;Plot
velovect2, vx_plot(1:sample-1,1:sample-1), vy_plot(1:sample-1,1:sample-1) $
, xaxis(1:sample-1), yaxis(1:sample-1),/overplot,length=length,color=color
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;                           MAIN PROGRAM                           ;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro vortxy2, loc=loc, start=start, finish=finish, xs=xs, tracer=tracer $
,plotrange=plotrange, xrange=xrange, ct=ct,yrange=yrange,stream=stream $
,vcol=vcol, vlen=vlen, sample=sample, neg=neg, jump=jump, mass=mass, h=h $
,order=order, eta=eta, out=out
common consts, pi
common grid, nsec, nrad, azi, rmed, dtheta, rmns, rpls, azinew, rplot
common hydro, vtheta, vrad
common planet, info,plx, ply, plrad
common misc, opt, factor0
common curve, coeff
;Constants
pi=3.141592654
factor0=double((mass/3.)^(1./3.))
;Options
if not keyword_set(ct) then ct=5
if not keyword_set(yrange) then yrange=[-1.0,1.0]
if not keyword_set(xs) then xs=2.0
if not keyword_set(xrange) then xrange=[-xs,xs]
if keyword_set(plotrange) then plotrange0=plotrange
if not keyword_set(vcol) then vcol=-1
if not keyword_set(vlen) then vlen=2
if not keyword_set(sample) then sample=30
if n_elements(start) ne 1 then begin
cases=start
endif
if not keyword_set(finish) then begin
cases=[start]
endif else cases=start+indgen(finish-start+1)
numcases=n_elements(cases)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;DATA LOCATION,GET RUN INFO, READ PLANET ORBIT INFO.;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
location=strcompress(loc,/remove_all)
dims=(read_ascii(filepath('dims.dat',root_dir='.',subdir=[location]))).(0)
nout=fix(dims(5))
nrad=fix(dims(6))
nsec=fix(dims(7))
if not keyword_set(out) then begin
    info=dblarr(11,nout+1)
endif else info=dblarr(11,out+1)
openr,3,filepath('planet0.dat',root_dir='.',subdir=[location])
readf,3,info
close,3
a0=info(1,0)
dt=info(7,1)
p0=2.*pi*(a0)^(3./2.)
;Create 1d arrays for radial and azimuthal values for contour
;plot,
azifix=dblarr(nsec)
azinew=dblarr(nsec)
radtmp=dblarr(nrad+1)
openr,1,filepath('used_rad.dat',root_dir='.',subdir=[location])
readf,1,radtmp
close,1
azifix=dindgen(nsec)*2.*pi/nsec
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;SETUP FOR VORTICITY CALCULATION;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;Arrays to hold variables.Also create 2d array to hold radii values
;over grid.
sigma=dblarr(nsec,nrad)
vtheta= dblarr(nsec,nrad)
vrad=dblarr(nsec,nrad)
label=dblarr(nsec,nrad)
rmed=0.5*(radtmp(0:nrad-1)+radtmp(1:nrad))
radius=dblarr(nsec,nrad)
for l=0, nsec-1 do radius(l,0:nrad-1)=rmed(0:nrad-1)
;Cell sizes.
dtheta=2.*pi/nsec
dr=(radtmp(nrad)-radtmp(0))/nrad
;Array to hold results.Don't calculate vorticity for first and last
;radius (boundaries).
vort=dblarr(nsec,nrad-2)
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
openr,3,filepath(strcompress('gasvtheta'+string(cases(k))+'.dat',/remove_all),root_dir='.',subdir=[location]) 
readu,3,vtheta
close,3
openr,4,filepath(strcompress('gasvrad'+string(cases(k))+'.dat',/remove_all),root_dir='.',subdir=[location]) 
readu,4,vrad
close,4
;Now calculate vorticity. We need r*vtheta. For convenience define
;this array. Then begin finite differencing, Convention is (i,j) for
;(theta,radius indicies). Start with boundary values (theta=0 and
;theta=2pi grid points). We on working on the data array of size
;(nsec,nrad), but calculating vorticity everywhere except first and
;last radius value. Then put results into the smaller array, vort.
rvtheta=radius*vtheta
;When i=0, j=1 to nrad-2
dr=radius(0,2:nrad-1)-radius(0,0:nrad-3)
vort(0,0:nrad-3)=(0.5*(rvtheta(0,2:nrad-1)+rvtheta(1,2:nrad-1))-0.5*(rvtheta(0,0:nrad-3)+rvtheta(1,0:nrad-3)))/(2.*dr) $
-(0.5*(vrad(1,1:nrad-2)+vrad(1,2:nrad-1))-0.5*(vrad(nsec-1,1:nrad-2)+vrad(nsec-1,2:nrad-1)))/(2.*dtheta)
;When i=nsec-1, j=1 to nrad-2
dr=radius(0,2:nrad-1)-radius(0,0:nrad-3)
vort(nsec-1,0:nrad-3)=(0.5*(rvtheta(nsec-1,2:nrad-1)+rvtheta(0,2:nrad-1))-0.5*(rvtheta(nsec-1,0:nrad-3)+rvtheta(0,0:nrad-3)))/(2.*dr) $
-(0.5*(vrad(0,1:nrad-2)+vrad(0,2:nrad-1))-0.5*(vrad(nsec-2,1:nrad-2)+vrad(nsec-2,2:nrad-1)))/(2.*dtheta)
;Rest of grid, when i=1 to nsec-2 and j=1 to nrad-2
dr=radius(1:nsec-2,2:nrad-1)-radius(1:nsec-2,0:nrad-3)
vort(1:nsec-2,0:nrad-3)=(0.5*(rvtheta(1:nsec-2,2:nrad-1)+rvtheta(2:nsec-1,2:nrad-1))-0.5*(rvtheta(1:nsec-2,0:nrad-3)+ $
                                                                                          rvtheta(2:nsec-1,0:nrad-3)))/(2.*dr) $
-(0.5*(vrad(2:nsec-1,1:nrad-2)+vrad(2:nsec-1,2:nrad-1))-0.5*(vrad(0:nsec-3,1:nrad-2)+vrad(0:nsec-3,2:nrad-1)))/(2.*dtheta)
;Still need to divide by radius at each grid to get vorticity,
vort(0:nsec-1,0:nrad-3)=vort(0:nsec-1,0:nrad-3)/radius(0:nsec-1,1:nrad-2)
;Get planet info, and get grid points corresponding to radii of
;interest
time=string(cases(k)*dt/p0,format='(F6.2)')
plx=info(1,cases(k))
ply=info(2,cases(k))
plrad=sqrt(plx*plx+ply*ply)
rhill=plrad*factor0
phi=pltphi(plx,ply)
rplot=(rmed-plrad)/rhill
azi=azifix-phi
arr1=min(abs(rplot-xrange(0)),rmns)
arr2=min(abs(rplot-xrange(1)),rpls)
arr3=min(abs(azinew-yrange(0)),phimns)
arr4=min(abs(azinew-yrange(1)),phipls)
;Now extract relavent data set
data=vort(0:nsec-1,0:nrad-30)/sigma(0:nsec-1,1:nrad-2);/vort(0:nsec-1,0:nrad-3);inverse vortensity
data2=transpose(data(0:nsec-1,rmns-1:rpls-1))
if not keyword_set(plotrange0) then begin
;plotrange=[alog(min(data2(where(data2 gt 0.)))),alog(max(data2))]
plotrange=[min(data2),max(data2)]*1d-2
endif else if plotrange0(0) eq 0 then begin 
plotrange=[alog(min(data2(where(data2 gt 0.)))),plotrange(1)]
endif
;Shift data so plot is centered about planet.
opt=whichcase(phi)
data3=shiftazi(data2,opt)
azitmp=shiftazi(azi,opt)
azinew(0:nsec-1)=azitmp(0,0:nsec-1)/pi
;Deal with negative sigma/omega by assigning them small values, save
;un-logged data:
invort=data3
zeros=where(data3 le 0.)
s=size(ind)
;if min(zeros) ne -1 then begin
;data3(zeros)=exp(-100.)
;ind=array_indices(data3,zeros)
;endif
;data3=alog(data3)
levels=plotrange(0)+(plotrange(1)-plotrange(0))*(dindgen(64)/63.)
;PLOT THIS FUCKER
set_plot, 'ps'
device, filename=filepath(strcompress('vortxy2_'+ks+'.ps',/remove_all),root_dir='.',subdir=[location])$
,/color, bits_per_pixel=8,xsize=12, ysize=15
temp=strcompress(string(plrad,format='(f5.2)'),/remove_all)
contour,data3*1d-2,rplot(rmns:rpls),azinew(0:nsec-1) ,/fill $
,levels=levels,title=time+' orbits,'+textoidl(' r_p=')+temp,ymargin=[4,2],xmargin=[8,10],color=0,xtickinterval=2.0$
,xtitle=xxtitle,ytitle=yytitle,charsize=1.5,yrange=yrange,xrange=xrange
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;Option to indicate negative values;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
if keyword_set(neg) then begin
for j=0, s(2)-1 do begin
oplot, [rplot(ind(0,j)+rmns),rplot(ind(0,j)+rmns)] $
,[azinew(ind(1,j)+phimns),azinew(ind(1,j)+phimns)], color=0, psym=1
endfor
endif
colorbar, position=[0.80, 0.18, 0.85, 0.88],/vertical,/right,range=plotrange,format='(F5.1)',color=0,charsize=1.5
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;Option to plot tracer fluid;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
if keyword_set(tracer) then begin
openr,5,filepath(strcompress('gaslabel'+string(cases(k))+'.dat',/remove_all),root_dir='.',subdir=[location])
readu,5,label
close,5
label2=label(0:nsec-1,1:nrad-2)
label3=transpose(label2(0:nsec-1,rmns:rpls))
beta=(max(label3)-min(label3))/(max(data2)-min(data2))
alpha=min(label3)-beta*min(data2)
label3=(label3-alpha)/beta
data3=rotate(label3,opt)
contour,data3,rplot(rmns:rpls),azinew(0:nsec-1),/overplot,color=-1, thick=0.5,nlevels=3;,/fill
endif
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;Option to overplot velocity field;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
if keyword_set(stream) then begin
velocity, slice=cases(k), sample=sample, length=vlen, yrange=yrange, color=vcol
endif
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;Option to get vortensity jump across shock (outer);
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
if keyword_set(jump) then begin
if n_elements(jump) eq 1 then jump=[jump,jump]
if not keyword_set(order) then order=5
shock_fit,loc=loc,start=cases(k),xrange=xrange ,yrange=yrange, eta=eta, order=order,/noplot, mass=mass, out=out
ord=n_elements(coeff)
x=rplot(rmns:rpls)*rhill
y=plrad*azinew*pi
;Axes and arrays
pts=floor((xrange(1)-xrange(0))*factor0*nrad/alog(10.0))
xshock=rhill*(xrange(0)+(xrange(1)-xrange(0))*dindgen(pts)/(pts-1.0))
y_s=dblarr(pts)
dys=dblarr(pts)
vortensity_jump=dblarr(pts)
v_normal=dblarr(pts)
pre_x=dblarr(pts)
xxx=dblarr(pts)
yyy=dblarr(pts)
;Get pre-shock flow field from data
polar_vel=convert(plrad,plx,ply,info(3,cases(k)),info(4,cases(k)))
v_x=shiftazi(transpose(vrad(0:nsec-1,rmns:rpls)),opt)-polar_vel(0)
v_y=shiftazi(transpose(vtheta(0:nsec-1,rmns:rpls)/radius(0:nsec-1,rmns:rpls)),opt)-polar_vel(1)
v_y*=plrad
vorticity=shiftazi(transpose(vort(0:nsec-1,rmns-1:rpls-1)),opt)
presigma=shiftazi(transpose(sigma(0:nsec-1,rmns:rpls)),opt)
;Shock curve y=y(x) and it's derivative
for i=0, pts-1 do begin
y_s(i)=0.0
dys(i)=0.0
for j=0, ord-1 do y_s(i)=y_s(i)+coeff(j)*xshock(i)^double(j)
for j=1, ord-1 do dys(i)=dys(i)+double(j)*coeff(j)*xshock(i)^double(j-1)
endfor

for i=0, pts-1 do begin
x0=xshock(i)
y0=y_s(i)
theta=atan(-1.0/dys(i))
dx=jump*rhill*cos(theta)
dy=jump*rhill*sin(theta)
x1=x0+dx(0);pre-shock sample x-coordinate
y1=y0+dy(0)
x2=x0-dx(1);post-shock sample x-coordinate
y2=y0-dy(1)
tmp1=min(abs(x-x1),xx1);pre-shock sampling grid point
tmp2=min(abs(y-y1),yy1)
tmp3=min(abs(x-x2),xx2);post-shock sampling grid point
tmp4=min(abs(y-y2),yy2)
oplot,[rplot(xx1+rmns),rplot(xx1+rmns)],[azinew(yy1),azinew(yy1)],color=-1
oplot,[rplot(xx2+rmns),rplot(xx2+rmns)],[azinew(yy2),azinew(yy2)],color=-1
;vortensity_jump(i)=1.0/invort(xx2,yy2)-1.0/invort(xx1,yy1)
vortensity_jump(i)=invort(xx2,yy2)-invort(xx1,yy1)
v_normal(i)=(-dys(i)*v_x(xx1,yy1)+v_y(xx1,yy1))/sqrt(1.0+dys(i)^2.0)
;Output information needed to take pre-shock flow from data
pre_x(i)=x1
xxx(i)=xx1
yyy(i)=yy1
endfor
vn_coeff=svdfit(pre_x,v_normal,11,/double);normal velocity as a function of x
;dvnormal=deriv(pre_x,v_normal)
openw,1,filepath(strcompress('vortensity_jump.dat',/remove_all),root_dir='.',subdir=[location])
for i=0, pts-1 do begin
 dvn=0.0
 for j=1, 10 do dvn=dvn+double(j)*vn_coeff(j)*pre_x(i)^double(j-1)
;dvn=dvnormal(i)
cs=h*(plrad+pre_x(i))^(-0.5)
dcs2=-h*h*(plrad+pre_x(i))^(-2.0)
machsq=(v_normal(i)/cs)^2.0
xdot=1.0/sqrt(1.0+dys(i)^2.0)
; domega=-xdot*dvn*(machsq-1.0)^2.0/machsq $
; +(machsq-1.0)*vorticity(xxx(i),yyy(i)) $
; -xdot*(machsq-1.0)*dcs2/v_normal(i)
; predicted=(vorticity(xxx(i),yyy(i))+domega)/(machsq*presigma(xxx(i),yyy(i)))-1.0/invort(xxx(i),yyy(i))
predicted=-(machsq-1.0)^2.0*dvn*xdot/(presigma(xxx(i),yyy(i))*machsq^2.0) $
-(machsq-1.0)*dcs2*xdot/(presigma(xxx(i),yyy(i))*machsq*v_normal(i))
printf,1,xshock(i),machsq,vortensity_jump(i),predicted,format='(4(e13.6E2,1x))'
endfor
close,1
endif
device,/close
print, 'done'+ks
endfor
end






































;Shift azimuthal values so that planet is at centre of plot
; data3=dblarr(rpls-rmns+1,nsec)
; azinew=dblarr(nsec)
; if phi gt pi then begin
; excess=where(azi lt -pi)
; azi(excess)=2.*pi+azi(excess)
; seg1=where(azi eq min(azi))
; if seg1 ne nsec-1 then begin
; data3(*,0:nsec-seg1-1)=data2(*,seg1:nsec-1)
; azinew(0:nsec-seg1-1)=azi(seg1:nsec-1)/pi
; endif
; if seg1 ne 1 then begin
; data3(*,nsec-seg1:nsec-1)=data2(*,0:seg1-1)
; azinew(nsec-seg1:nsec-1)=azi(0:seg1-1)/pi
; endif
; endif else if (0. le phi and  phi lt pi and (pi-phi) gt dtheta ) then begin
; excess=where(azi ge pi)
; azi(excess)=azi(excess)-2.*pi
; seg2=where(azi eq max(azi))
; if seg2 ne nsec-2 then begin
; data3(*,0:nsec-seg2-2)=data2(*,seg2+1:nsec-1)
; azinew(0:nsec-seg2-2)=azi(seg2+1:nsec-1)/pi
; endif
; if seg2 ne 0 then begin
; data3(*,nsec-seg2-1:nsec-1)=data2(*,0:seg2)
; azinew(nsec-seg2-1:nsec-1)=azi(0:seg2)/pi
; endif
; endif else begin
; data3=data2
; azinew(0:nsec-1)=azi(0:nsec-1)/pi
; endelse

;OLD PLOTTING ROUTINE
; if phi gt pi then begin
; excess=where(azi lt -pi)
; azi(excess)=2.*pi+azi(excess)
; seg1=where(azi eq min(azi))
; if seg1 ne nsec-1 then begin
; contour,data2(*,seg1:nsec-1),rplot(rmns:rpls),azi(seg1:nsec-1)/pi,/fill $
; ,levels=levels,title=time+' orbits',ymargin=[4,2],xmargin=[8,10],color=0$
; ,xtitle=xxtitle,ytitle=yytitle,charsize=1.5,yrange=[-1.,1.], xrange=xrange
; endif
; if seg1 ne 1 then begin
; contour,data2(*,0:seg1-1),rplot(rmns:rpls),azi(0:seg1-1)/pi,/fill $
; ,levels=levels,/overplot,yrange=[-1.,1.],title=time+' orbits',ymargin=[4,2],xmargin=[8,10],color=0$
; ,xtitle=xxtitle,ytitle=yytitle,charsize=1.5, xrange=xrange
; endif
; endif else if (0. le phi and  phi lt pi and (pi-phi) gt dtheta ) then begin
; excess=where(azi ge pi)
; azi(excess)=azi(excess)-2.*pi
; seg2=where(azi eq max(azi))
; if seg2 ne nsec-2 then begin
; contour,data2(*,seg2+1:nsec-1),rplot(rmns:rpls),azi(seg2+1:nsec-1)/pi,/fill $
; ,levels=levels,title=time+' orbits',ymargin=[4,2],xmargin=[8,10],color=0$
; ,xtitle=xxtitle,ytitle=yytitle,charsize=1.5,yrange=[-1.,1.], xrange=xrange
; endif
; if seg2 ne 0 then begin
; contour,data2(*,0:seg2),rplot(rmns:rpls),azi(0:seg2)/pi,/fill $
; ,levels=levels,/overplot,yrange=[-1.,1.],title=time+' orbits',ymargin=[4,2],xmargin=[8,10],color=0$
; ,xtitle=xxtitle,ytitle=yytitle,charsize=1.5, xrange=xrange
; endif
; endif else begin
; contour,data2,rplot(rmns:rpls),azi(0:nsec-1)/pi,/fill $
; ,levels=levels,title=time+' orbits',ymargin=[4,2],xmargin=[8,10],color=0$
; ,xtitle=xxtitle,ytitle=yytitle,charsize=1.5,yrange=[-1.,1.], xrange=xrange
; endelse
; ;did not plot points in order of index. may be gap so manually insert
; ;line segment
; ;arrow, azi(0),data1d(0),azi(nsec-1),data1d(nsec-1),hsize=0,/data
; data3=dblarr(rpls-rmns+1,2)
; data3(*,0)=data2(*,nsec-1)
; data3(*,1)=data2(*,0)
; contour,data3(*,0:1),rplot(rmns:rpls),[azi(nsec-1),azi(0)]/pi,/fill,levels=levels,/overplot
; if phi gt pi then begin
; if seg1 ne nsec-1 then begin
; contour,label3(*,seg1:nsec-1),rplot(rmns:rpls),azi(seg1:nsec-1)/pi,/overplot,color=-1,nlevels=3
; endif
; if seg1 ne 1 then begin
; contour,label3(*,0:seg1-1),rplot(rmns:rpls),azi(0:seg1-1)/pi,/overplot,color=-1,nlevels=3
; endif
; endif else if (0. le phi and  phi lt pi and (pi-phi) gt dtheta ) then begin
; if seg2 ne nsec-2 then begin
; contour,label3(*,seg2+1:nsec-1),rplot(rmns:rpls),azi(seg2+1:nsec-1)/pi,/overplot,color=-1,nlevels=3
; endif
; if seg2 ne 0 then begin
; contour,label3(*,0:seg2),rplot(rmns:rpls),azi(0:seg2)/pi,/overplot,color=-1,nlevels=3
; endif
; endif else begin
; contour,label3,rplot(rmns:rpls),azi(0:nsec-1)/pi,/overplot,color=-1,nlevels=3
; endelse
; endif
