; THIS PROCEDURE CALCULATES SIGMA/OMEGA ON THE GRID POINTS USING
; CENTRAL DIFFERENCING, AND DOES A POLAR CONTOUR PLOT. OPTION TO
; ROTATE DISC SUCH THAT PLANET IS AT AZIMUTH=0. THIS VERSION
; DISCRETISES EXPRESSION FOR VORTICITY IN THE FORM: 
; omega=(1/r)*(utheta+r*d(utheta)/dr)-(1/r)*d(urad)/dtheta
pro vort2, loc=loc, start=start, finish=finish, pframe=pframe
common consts, pi
pi=3.141592654
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;DATA LOCATION,GET RUN INFO, READ PLANET ORBIT INFO.;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
location=strcompress('out'+loc,/remove_all)
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
;Create 1d arrays for radial and azimuthal values for contour
;plot. Extend radius down to zero to allow contour plot to have a
;'hole' in the middle.
in=4
azi=dblarr(nsec)
radtmp=dblarr(nrad+1)
rad=dblarr(in+nrad)
openr,1,filepath('used_rad.dat',root_dir='.',subdir=[location])
readf,1,radtmp
close,1
azi=dindgen(nsec)*2.*pi/nsec
rad(0:in-1)=(radtmp(0)/in)*(dindgen(in)+0.5)
for j=in, nrad+in-1 do rad(j)=(radtmp(j-in)+radtmp(j-in+1))/2.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;SETUP FOR VORTICITY CALCULATION;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;Arrays to hold variables.Also create 2d array to hold radii values
;over grid.
sigma=dblarr(nsec,nrad)
vtheta= dblarr(nsec,nrad)
vrad=dblarr(nsec,nrad)
rmed=radtmp(0:nrad-1)+radtmp(1:nrad)
rmed=rmed*0.5
radius=dblarr(nsec,nrad)
for l=0, nsec-1 do radius(l,0:nrad-1)=rmed(0:nrad-1)
;Cell sizes.
dtheta=2.*pi/nsec
dr=(radtmp(nrad)-radtmp(0))/nrad
;Array to hold results.Don't calculate vorticity for first and last
;radius (boundaries).
vort=dblarr(nsec,nrad-2)
data2=dblarr(nsec,nrad-2+in+1)
azi1=dblarr(nsec)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;SETUP DONE, START CALCULATION AND PLOTTING;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
loadct,3, bottom=0
for k=start, finish do begin
ks=string(k,format='(I03)')
;Read raw data.
openr,2,filepath(strcompress('gasdens'+string(k)+'.dat',/remove_all),root_dir='.',subdir=[location]) 
readu,2,sigma
close,2
openr,3,filepath(strcompress('gasvtheta'+string(k)+'.dat',/remove_all),root_dir='.',subdir=[location]) 
readu,3,vtheta
close,3
openr,4,filepath(strcompress('gasvrad'+string(k)+'.dat',/remove_all),root_dir='.',subdir=[location]) 
readu,4,vrad
close,4
;Now calculate vorticity.Convention is (i,j) for
;(theta,radius) indicies. Start with boundary values (theta=0 and
;theta=2pi grid points). We on working on the data array of size
;(nsec,nrad), but calculating vorticity everywhere except first and
;last radius value. Then put results into the smaller array, vort.
;When i=0, j=1 to nrad-2
vort(0,0:nrad-3)=(0.5*(vtheta(0,2:nrad-1)+vtheta(1,2:nrad-1))-0.5*(vtheta(0,0:nrad-3)+vtheta(1,0:nrad-3)))/(2.*dr) $
+0.5*(vtheta(0,1:nrad-2)+vtheta(1,1:nrad-2))/radius(0,1:nrad-2) $
-(0.5*(vrad(1,1:nrad-2)+vrad(1,2:nrad-1))-0.5*(vrad(nsec-1,1:nrad-2)+vrad(nsec-1,2:nrad-1)))/(2.*dtheta*radius(0,1:nrad-2))
;When i=nsec-1, j=1 to nrad-2
vort(nsec-1,0:nrad-3)=(0.5*(vtheta(nsec-1,2:nrad-1)+vtheta(0,2:nrad-1))-0.5*(vtheta(nsec-1,0:nrad-3)+vtheta(0,0:nrad-3)))/(2.*dr) $
+0.5*(vtheta(nsec-1,1:nrad-2)+vtheta(0,1:nrad-2))/radius(nsec-1,1:nrad-2) $
-(0.5*(vrad(0,1:nrad-2)+vrad(0,2:nrad-1))-0.5*(vrad(nsec-2,1:nrad-2)+vrad(nsec-2,2:nrad-1)))/(2.*dtheta*radius(nsec-1,1:nrad-2))
;Rest of grid, when i=1 to nsec-2 and j=1 to nrad-2
vort(1:nsec-2,0:nrad-3)=(0.5*(vtheta(1:nsec-2,2:nrad-1)+vtheta(2:nsec-1,2:nrad-1))-0.5*(vtheta(1:nsec-2,0:nrad-3)+ $
vtheta(2:nsec-1,0:nrad-3)))/(2.*dr) $
+0.5*(vtheta(1:nsec-2,1:nrad-2)+vtheta(2:nsec-1,1:nrad-2))/radius(1:nsec-2,1:nrad-2) $
-(0.5*(vrad(2:nsec-1,1:nrad-2)+vrad(2:nsec-1,2:nrad-1))-0.5*(vrad(0:nsec-3,1:nrad-2)+vrad(0:nsec-3,2:nrad-1)))/(2.*dtheta*radius(1:nsec-2,1:nrad-2))
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;NOW DO THE POLAR CONTOUR PLOTS;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
data=sigma(0:nsec-1,1:nrad-2)/vort(0:nsec-1,0:nrad-3)
plotrange=[alog(min(data(where(data gt 0.)))),alog(max(data))]
;Deal with negative sigma/omega by assigning them small values:
zeros=where(data le 0.)
if min(zeros) ne -1 then data(zeros)=exp(-10.)
;Fill in data2, which is the array used for plotting.
data=alog(data)
data2(*,0:in)=10.
data2(*,in+1:nrad+in-2)=data
levels=plotrange(0)+(plotrange(1)-plotrange(0))*(dindgen(64)/64.)
time=string(k*dt/p0,format='(F5.1)')
;Set pframe to rotate disc so planet stays at x=0. (Should set to
;default in future.
if keyword_set(pframe) then begin
azi1=dblarr(nsec)
plx=info(1,k)
ply=info(2,k)
phi=pltphi(plx,ply)
array=abs(azi-phi)
grid=where(array eq min(array))
for l=0, nsec-1 do azi1(l)=azi(l)-azi(grid)
endif else  azi1(0:nsec-1)=azi(0:nsec-1)
;The actual plotting part.
set_plot, 'ps'
device, filename=filepath(strcompress('vort2_'+ks+'.ps',/remove_all),root_dir='.',subdir=[location])$
,/color, bits_per_pixel=8,xsize=12, ysize=9
 polar_contour,data2,azi1,rad(0:nrad+in-2),/isotropic,/fill,levels=levels,title=time+' orbits',ymargin=[2.5,2.5],xmargin=[6,6]
 colorbar, position=[0.75, 0.1, 0.80, 0.9],/vertical,/right,range=plotrange,format='(F4.1)'
;Option to explicity show where planet is:
  plrad=sqrt(plx*plx+ply*ply)
  oplot,[plrad,plrad],[0,0],psym=6,symsize=1
 device,/close
 print, 'done '+ks
endfor
; rr=abs(rmed-plrad)
; print, mean(data(*,where(rr eq min(rr))))
end
