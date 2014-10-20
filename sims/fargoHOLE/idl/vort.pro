; THIS PROCEDURE CALCULATES SIGMA/OMEGA ON THE GRID POINTS USING
; CENTRAL DIFFERENCING, AND DOES A POLAR CONTOUR PLOT. OPTION TO
; ROTATE DISC SUCH THAT PLANET IS AT AZIMUTH=0.
; 20/12/2008: Added option to overplot tracer fluid density. But currently disabled (commented out).
pro vort, loc=loc, start=start, finish=finish,plotrange0=plotrange0,ct=ct,out=out
common consts, pi
pi=3.141592654
if not keyword_set(ct) then ct=5
if not keyword_set(finish) then finish=start
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
;for j=in, nrad+in-1 do rad(j)=(radtmp(j-in)+radtmp(j-in+1))/2.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;SETUP FOR VORTICITY CALCULATION;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;Arrays to hold variables.Also create 2d array to hold radii values
;over grid.
sigma=dblarr(nsec,nrad)
vtheta= dblarr(nsec,nrad)
vrad=dblarr(nsec,nrad)
dr2d=dblarr(nsec,nrad)
;label=dblarr(nsec,nrad)
;rmed=radtmp(0:nrad-1)+radtmp(1:nrad)
;rmed=rmed*0.5
rmed = 2.0/3.0*(radtmp(1:nrad)^3.-radtmp(0:nrad-1)^3.);
rmed(0:nrad-1) = rmed(0:nrad-1) / (radtmp(1:nrad)^2.-radtmp(0:nrad-1)^2.)
for j=in, nrad+in-1 do rad(j)=rmed(j-in)
radius=dblarr(nsec,nrad)
for l=0, nsec-1 do radius(l,0:nrad-1)=rmed(0:nrad-1)
;Cell sizes.
dtheta=2.*pi/nsec
dr=(radtmp(nrad)-radtmp(0))/nrad
for i=0, nsec-1 do dr2d(i,*)=dr(*)
;Array to hold results.Don't calculate vorticity for first and last
;radius (boundaries).
vort=dblarr(nsec,nrad-2)
data2=dblarr(nsec,nrad-2+in+1)
azi1=dblarr(nsec)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;SETUP DONE, START CALCULATION AND PLOTTING;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
loadct,ct, bottom=0
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
;if keyword_set(tracer) then begin
;openr,5,filepath(strcompress('gaslabel'+string(k)+'.dat',/remove_all),root_dir='.',subdir=[location]) 
;readu,5,label
;close,5
;endif
;Now calculate vorticity. We need r*vtheta. For convenience define
;this array. Then begin finite differencing, Convention is (i,j) for
;(theta,radius indicies). Start with boundary values (theta=0 and
;theta=2pi grid points). We on working on the data array of size
;(nsec,nrad), but calculating vorticity everywhere except first and
;last radius value. Then put results into the smaller array, vort.
rvtheta=radius*vtheta
;When i=0, j=1 to nrad-2
dr=radius(0,2:nrad-1)-radius(0,0:nrad-3)
vort(0,0:nrad-3)=(0.5*(rvtheta(0,2:nrad-1)+rvtheta(1,2:nrad-1))-0.5*(rvtheta(0,0:nrad-3)+rvtheta(1,0:nrad-3)))/(dr) $
-(0.5*(vrad(1,1:nrad-2)+vrad(1,2:nrad-1))-0.5*(vrad(nsec-1,1:nrad-2)+vrad(nsec-1,2:nrad-1)))/(2.*dtheta)
;When i=nsec-1, j=1 to nrad-2
dr=radius(0,2:nrad-1)-radius(0,0:nrad-3)
vort(nsec-1,0:nrad-3)=(0.5*(rvtheta(nsec-1,2:nrad-1)+rvtheta(0,2:nrad-1))-0.5*(rvtheta(nsec-1,0:nrad-3)+rvtheta(0,0:nrad-3)))/(dr) $
-(0.5*(vrad(0,1:nrad-2)+vrad(0,2:nrad-1))-0.5*(vrad(nsec-2,1:nrad-2)+vrad(nsec-2,2:nrad-1)))/(2.*dtheta)
;Rest of grid, when i=1 to nsec-2 and j=1 to nrad-2
dr=radius(1:nsec-2,2:nrad-1)-radius(1:nsec-2,0:nrad-3)
vort(1:nsec-2,0:nrad-3)=(0.5*(rvtheta(1:nsec-2,2:nrad-1)+rvtheta(2:nsec-1,2:nrad-1))-0.5*(rvtheta(1:nsec-2,0:nrad-3)+ $
rvtheta(2:nsec-1,0:nrad-3)))/dr $
-(0.5*(vrad(2:nsec-1,1:nrad-2)+vrad(2:nsec-1,2:nrad-1))-0.5*(vrad(0:nsec-3,1:nrad-2)+vrad(0:nsec-3,2:nrad-1)))/(2.*dtheta)
;Still need to divide by radius at each grid to get vorticity,
vort(0:nsec-1,0:nrad-3)=vort(0:nsec-1,0:nrad-3)/radius(0:nsec-1,1:nrad-2)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;NOW DO THE POLAR CONTOUR PLOTS;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;data=sigma(0:nsec-1,1:nrad-2)/vort(0:nsec-1,0:nrad-3)
data=vort(0:nsec-1,0:nrad-3)/sigma(0:nsec-1,1:nrad-2)

if not keyword_set(plotrange0) then begin
    plotrange=[alog10(min(data(where(data gt 0.)))),alog10(max(data))]
    ;plotrange=[min(data), max(data)]
endif else plotrange=plotrange0

;Deal with negative sigma/omega by assigning them small values:
zeros=where(data le 0.)
if min(zeros) ne -1 then data(zeros)=1d-10
;Fill in data2, which is the array used for plotting.
data=alog10(data)

data2=dblarr(nsec,nrad)
rad2 = dblarr(nrad)

for i=0, nsec-1 do data2(i,0:1)=10.0+10.0*dindgen(2)/2.0
data2(*,2:nrad-1)=data
rad2(0:1)=[0.0,rmed(1)/2.0]
rad2(2:nrad-1) = rmed(1:nrad-2)

levels=plotrange(0)+(plotrange(1)-plotrange(0))*(dindgen(16)/16.)
time=string(k*dt/p0,format='(F7.2)')
;Set pframe to rotate disc so planet stays at y=0.
plx=info(1,k)
ply=info(2,k)
phi=pltphi(plx,ply)
;array=abs(azi-phi)
;grid=where(array eq min(array))
for l=0, nsec-1 do azi1(l)=azi(l)-phi
;The actual plotting part.
set_plot, 'ps'
device, filename=filepath(strcompress('vort_'+ks+'.ps',/remove_all),root_dir='.',subdir=[location])$
,/color, bits_per_pixel=8,xsize=14, ysize=12
 polar_contour,data2(*,*),azi1(*),rad2,/isotropic,levels=levels,title=time+' orbits',ymargin=[2.5,2.5],xmargin=[6,6],color=0 , xtickinterval=5, ytickinterval=5, /dither;,/fill

;if keyword_set(tracer) then begin
;polar_contour,label,azi1,rmed,/isotropic,/fill,levels=levels,/overplot
;endif
 colorbar, position=[0.85, 0.1, 0.90, 0.9],/vertical,/right,range=plotrange,format='(F5.1)',color=0
;Show where planet is:
plrad=sqrt(plx*plx+ply*ply)
oplot,[plrad,plrad],[0,0],psym=6,symsize=1,color=-1
device,/close
print, 'done '+ks
endfor
end
