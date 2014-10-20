; THIS PROCEDURE CALCULATES OMEGA ON THE GRID POINTS USING
; CENTRAL DIFFERENCING. THIS VERSION IS FOR 1D PLOT OF
; VORTICITY V.S. RADIUS
pro omgrad, loc=loc, start=start, finish=finish, plazi=plazi, opp=opp
common consts, pi, nrad
common array, avgsig, data1d, rad
pi=3.141592654
xxtitle='r'
if keyword_set(plazi) then begin
yytitle=textoidl('\omega')
endif else yytitle=textoidl('<\omega>_{\phi}')
if keyword_set(opp) then begin
opp=double(opp*pi)
endif else opp=0.0
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;DATA LOCATION, GET INFO ABOUT RUN, READ PLANET ORBIT INFO, AND GET;
;RADIUS AND AZIMUTHAL ARRAY VALUES OF THE GRID;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
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
radtmp=dblarr(nrad+1)
openr,1,filepath('used_rad.dat',root_dir='.',subdir=[location])
readf,1,radtmp
close,1
rmed=radtmp(0:nrad-1)+radtmp(1:nrad)
rmed=rmed*0.5
azi=dindgen(nsec)*2.*pi/nsec
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;SETUP FOR VORTICITY CALCULATION;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;Arrays to hold variables. vtheta and vrad is bigger than array of raw data because
;we need some ghost cells for numerical differentiation.
sigma=dblarr(nsec,nrad)
vtheta= dblarr(nsec,nrad)
vrad=dblarr(nsec,nrad)
;Grid size for numerical differentiation. 
dtheta=2.*pi/nsec
dr=(radtmp(nrad)-radtmp(0))/nrad
;Array to hold results and that for plotting. ignore innermost and
;outermost annulus. also create 2d array to hold radii.
vort=dblarr(nsec,nrad-2)
data1d=dblarr(nrad-2)
radius=dblarr(nsec,nrad)
for l=0, nsec-1 do radius(l,0:nrad-1)=rmed(0:nrad-1)
;SETUP DONE, START CALCULATION AND PLOTTING.
for k=start, finish do begin
ks=string(k,format='(I03)')
;Read raw data.
openr,2,filepath(strcompress('gasdens'+string(k)+'.dat',/remove_all),root_dir='.',subdir=[location]) 
readu,2,sigma
close,2
openr,3,filepath(strcompress('gasvtheta'+string(k)+'.dat',/remove_all),root_dir='.',subdir=[location]) 
readu,3,vtheta
close,3
openr,2,filepath(strcompress('gasvrad'+string(k)+'.dat',/remove_all),root_dir='.',subdir=[location]) 
readu,2,vrad
close,2
;can now calculate vorticity. first do boundary values. convention is
;(i,j) for (theta, radius) indicies. we need r*vtheta. for convenience
;define this array.
rvtheta=radius*vtheta
;when i=0, j=1 to nrad-2
vort(0,0:nrad-3)=(0.5*(rvtheta(0,2:nrad-1)+rvtheta(1,2:nrad-1))-0.5*(rvtheta(0,0:nrad-3)+rvtheta(1,0:nrad-3)))/(2.*dr) $
-(0.5*(vrad(1,1:nrad-2)+vrad(1,2:nrad-1))-0.5*(vrad(nsec-1,1:nrad-2)+vrad(nsec-1,2:nrad-1)))/(2.*dtheta)
;when i=nsec-1, j=1 to nrad-2
vort(nsec-1,0:nrad-3)=(0.5*(rvtheta(nsec-1,2:nrad-1)+rvtheta(0,2:nrad-1))-0.5*(rvtheta(nsec-1,0:nrad-3)+rvtheta(0,0:nrad-3)))/(2.*dr) $
-(0.5*(vrad(0,1:nrad-2)+vrad(0,2:nrad-1))-0.5*(vrad(nsec-2,1:nrad-2)+vrad(nsec-2,2:nrad-1)))/(2.*dtheta)
;rest of grid, when i=1 to nsec-2 and j=1 to nrad-2
vort(1:nsec-2,0:nrad-3)=(0.5*(rvtheta(1:nsec-2,2:nrad-1)+rvtheta(2:nsec-1,2:nrad-1))-0.5*(rvtheta(1:nsec-2,0:nrad-3)+ $
rvtheta(2:nsec-1,0:nrad-3)))/(2.*dr) $
-(0.5*(vrad(2:nsec-1,1:nrad-2)+vrad(2:nsec-1,2:nrad-1))-0.5*(vrad(0:nsec-3,1:nrad-2)+vrad(0:nsec-3,2:nrad-1)))/(2.*dtheta)
;actually still need to divide by radius at each grid to get vorticity,
vort(0:nsec-1,0:nrad-3)=vort(0:nsec-1,0:nrad-3)/radius(0:nsec-1,1:nrad-2)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;PLOT AZIMUTHAL-AVERAGED VORTICITY VS RADIUS;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
data=vort(0:nsec-1,0:nrad-3)
plx=info(1,k)
ply=info(2,k)
if keyword_set(plazi) then begin
phi=(pltphi(plx,ply)+opp) mod (2.*pi)
array=abs(azi-phi)
grid=where(array eq min(array))
data1d(0:nrad-3)=data(grid,0:nrad-3)
endif else begin
for m=0, nrad-3 do data1d(m)=mean(data(0:nsec-1,m))
endelse
plotrange=[0.,max(data1d)]
;plotrange=[0.05,0.5]
time=string(k*dt/p0,format='(F5.1)')
set_plot, 'ps'
device, filename=filepath(strcompress('omgrad'+ks+'.ps',/remove_all),root_dir='.',subdir=[location])$
,bits_per_pixel=8,xsize=8, ysize=6,xoffset=0,yoffset=0,/inches,/color
plot,rmed(1:nrad-2),data1d(0:nrad-3),title=time+' orbits',xmargin=[8,6],ymargin=[3,3],xtitle=xxtitle,ytitle=yytitle $
,charsize=1.5,xrange=[min(rmed),max(rmed)],xtickinterval=0.4,xminor=4,yrange=plotrange
;;;;;;;;;;;;;;;;;;;;
;MARK PLANET RADIUS;
;;;;;;;;;;;;;;;;;;;;
pos=0.75*(plotrange(1)-plotrange(0))
plrad=sqrt(plx*plx+ply*ply)
array2=abs(rmed-plrad)
vline, plrad
device,/close
print, 'done '+ks
endfor
;print, data(grid,where(array2 eq min(array2)))
end

