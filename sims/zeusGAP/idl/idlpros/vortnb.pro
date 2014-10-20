; THIS PROCEDURE CALCULATES SIGMA/OMEGA ON THE GRID POINTS USING
; CENTRAL DIFFERENCING
pro vortnb, loc=loc, start=start, finish=finish, pframe=pframe
common consts, pi
pi=3.141592654
;;;;;;;;;;;;;;;
;WHERE IS DATA;
;;;;;;;;;;;;;;;
location=strcompress('out'+loc,/remove_all)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;GET DIMENSIONS AND TIME UNIT;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
dims=(read_ascii(filepath('dims.dat',root_dir='.',subdir=[location]))).(0)
nout=fix(dims(5))
nrad=fix(dims(6))
nsec=fix(dims(7))
info=dblarr(9,nout+1)
openr,3,filepath('planet0.dat',root_dir='.',subdir=[location])
readf,3,info
close,3
a0=info(1,0)
dt=info(7,1)
p0=2.*pi*(a0)^(3./2.)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;CREATE ARRAY FOR AZIMUTH AND RADIAL (for later polar plot);
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
in=4
azi=dblarr(nsec)
radtmp=dblarr(nrad+1)
rad=dblarr(in+nrad)
openr,1,filepath('used_rad.dat',root_dir='.',subdir=[location])
readf,1,radtmp
close,1
;for i=0, nsec-1 do azi(i)=2.*pi*i/nsec
azi=dindgen(nsec)*2.*pi/nsec
rad(0:in-1)=(radtmp(0)/in)*(dindgen(in)+0.5)
for j=in, nrad+in-1 do rad(j)=(radtmp(j-in)+radtmp(j-in+1))/2.
;;;;;;;;;;;;;;;;;;;;;;;;
;DO POLAR CONTOUR PLOTS;
;;;;;;;;;;;;;;;;;;;;;;;;
;arrays to hold variables. vtheta and vrad is bigger than array of raw data because
;we need some ghost cells for numerical differentiation
sigma=dblarr(nsec,nrad)
vtheta= dblarr(nsec,nrad)
vrad=dblarr(nsec,nrad)
;radius array. radius of the (k,l) grid only depends on l
rmed=radtmp(0:nrad-1)+radtmp(1:nrad)
rmed=rmed*0.5
;;;;;stuff for  numerical differentiation. 
dtheta=2.*pi/nsec
dr=(radtmp(nrad)-radtmp(0))/nrad
;;;;array for results
vort=dblarr(nsec-2,nrad-2)
loadct,3, bottom=0
for k=start, finish do begin
ks=string(k,format='(I03)')
;;;;;;;;getting data
openr,2,filepath(strcompress('gasdens'+string(k)+'.dat',/remove_all),root_dir='.',subdir=[location]) 
readu,2,sigma
close,2
openr,3,filepath(strcompress('gasvtheta'+string(k)+'.dat',/remove_all),root_dir='.',subdir=[location]) 
readu,3,vtheta
close,3
openr,2,filepath(strcompress('gasvrad'+string(k)+'.dat',/remove_all),root_dir='.',subdir=[location]) 
readu,2,vrad
close,2
radius=dblarr(nsec,nrad)
;radius is 2d array to hold radius of each grid point. these do not
;vary with the first index (azimuth)
for l=0, nsec-1 do radius(l,0:nrad-1)=rmed(0:nrad-1)
;print, 'created arrays with ghost cells for'
;can now calculate vorticity
vort(0:nsec-3,0:nrad-3)=(radius(1:nsec-2,2:nrad-1)*vtheta(1:nsec-2,2:nrad-1)-radius(1:nsec-2,0:nrad-3)*vtheta(1:nsec-2,0:nrad-3))/(2.*dr)$
-(vrad(2:nsec-1,1:nrad-2)-vrad(0:nsec-3,1:nrad-2))/(2.*dtheta)
vort(0:nsec-3,0:nrad-3)=vort(0:nsec-3,0:nrad-3)/radius(1:nsec-2,1:nrad-2)
;;;;;below is just a different method to calculate vort using one loop
; for j=1, nrad do begin
; vort(0:nsec-1,j-1)=(rad2(j+1)*vtheta(1:nsec,j+1)-rad2(j-1)*vtheta(1:nsec,j-1))/(2.*dr)$
; -(vrad(2:nsec+1,j)-vrad(0:nsec-1,j))/(2.*dtheta)
; vort(0:nsec-1,j-1)=vort(0:nsec-1,j-1)/rad2(j)
; endfor
;print, 'calculated vorcitity array'
;;;;;;;;;plotting vort;;;;;;;;;;;;;
data2=dblarr(nsec,nrad+in-1)
azi1=dblarr(nsec)
data=sigma(1:nsec-2,1:nrad-2)/vort(0:nsec-3,0:nrad-3)
data2(*,0:in)=1.
data2(0,*)=1.
data2(nsec-1,*)=1.
data2(1:nsec-2,in+1:nrad+in-2)=data
plotrange=[0.,0.01];fixed-value entries for fixed colorbar, set to data(min), data(max) for variable
levels=plotrange(0)+(plotrange(1)-plotrange(0))*(dindgen(64)/64.)
time=string(k*dt/p0,format='(F5.1)')
if keyword_set(pframe) then begin
azi1=dblarr(nsec)
plx=info(1,k)
ply=info(2,k)
phi=pltphi(plx,ply)
array=abs(azi-phi)
grid=where(array eq min(array))
for l=0, nsec-1 do azi1(l)=azi(l)-azi(grid)
endif else  azi1(0:nsec-1)=azi(0:nsec-1)
set_plot, 'ps'
device, filename=filepath(strcompress('sigvort'+ks+'.ps',/remove_all),root_dir='.',subdir=[location])$
,/color, bits_per_pixel=8,xsize=12, ysize=9
 polar_contour,data2,azi1,rad(0:nrad+in-2),/isotropic,/fill,levels=levels,title=time+' orbits',ymargin=[2.5,2.5]
 colorbar, position=[0.85, 0.1, 0.9, 0.9],/vertical,/right,range=plotrange,format='(F5.3)'
 device,/close
 print, 'done '+ks
endfor
end
