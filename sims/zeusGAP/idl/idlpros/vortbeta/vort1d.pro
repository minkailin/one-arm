; THIS PROCEDURE CALCULATES SIGMA/OMEGA ON THE GRID POINTS USING
; CENTRAL DIFFERENCING. THIS VERSION IS FOR 1D PLOTS.
pro vort1d, loc=loc, start=start, finish=finish, planet=planet, plazi=plazi
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
;;;;;;;;;;;;;;;;;;;;;;;;;;;
;READ THE VALUES OF RADIUS;
;;;;;;;;;;;;;;;;;;;;;;;;;;;
radtmp=dblarr(nrad+1)
openr,1,filepath('used_rad.dat',root_dir='.',subdir=[location])
readf,1,radtmp
close,1
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;CREATE ARRAY OF AZIMUTHAL VALUES;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
azi=dindgen(nsec)*2.*pi/nsec
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;CALCULATE VORTICITY FOR EACH TIMESTEP;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;Arrays to hold variables. vtheta and vrad is bigger than array of raw data because
;we need some ghost cells for numerical differentiation.
sigma=dblarr(nsec,nrad)
vtheta= dblarr(nsec+2,nrad+2)
vrad=dblarr(nsec+2,nrad+2)
vthetaraw=dblarr(nsec,nrad)
vradraw=dblarr(nsec,nrad)
;Radius array. Radius of the (m,n) grid only depends on n.
rmed=radtmp(0:nrad-1)+radtmp(1:nrad)
rmed=rmed*0.5
;Grid size for numerical differentiation. 
dtheta=2.*pi/nsec
dr=(radtmp(nrad)-radtmp(0))/nrad
;Array to hold results.
vort=dblarr(nsec,nrad)
;;;;;;SETUP DONE;;;;;
xxtitle='r'
if keyword_set(plazi) then begin
yytitle=textoidl('\Sigma/\omega')
endif else yytitle=textoidl('<\Sigma/\omega>_{\phi}')
for k=start, finish do begin
ks=string(k,format='(I03)')
;Read raw data.
openr,2,filepath(strcompress('gasdens'+string(k)+'.dat',/remove_all),root_dir='.',subdir=[location]) 
readu,2,sigma
close,2
openr,3,filepath(strcompress('gasvtheta'+string(k)+'.dat',/remove_all),root_dir='.',subdir=[location]) 
readu,3,vthetaraw
close,3
vtheta(1:nsec,1:nrad)=vthetaraw(0:nsec-1,0:nrad-1)
openr,2,filepath(strcompress('gasvrad'+string(k)+'.dat',/remove_all),root_dir='.',subdir=[location]) 
readu,2,vradraw
close,2
vrad(1:nsec,1:nrad)=vradraw(0:nsec-1,0:nrad-1)
;Fill in ghost cells that is needed when doing azimuthal gradient at
;theta=0 or theta=2pi. 
vtheta(0,1:nrad)=vtheta(nsec,1:nrad)
vtheta(nsec+1,1:nrad)=vtheta(0,1:nrad)
vrad(0,1:nrad)=vrad(nsec,1:nrad)
vrad(nsec+1,1:nrad)=vrad(0,1:nrad)
;Fill in ghost cells beyond rmax and below rmin. This is needed when
;doing radialg radients at rmax and rmin. The values are determined
;such that the boundary value (data we have) is the average of it's 4
;neighbours. One of the neighbour is beyond the boundary.
vtheta(1:nsec,0)=4.*vtheta(1:nsec,1)-vtheta(0:nsec-1,1)-vtheta(1:nsec,2)-vtheta(2:nsec+1,1)
vtheta(1:nsec,nrad+1)=4.*vtheta(1:nsec,nrad)-vtheta(2:nsec+1,nrad)-vtheta(0:nsec-1,nrad)-vtheta(1:nsec,nrad-1)
vrad(1:nsec,0)=4.*vrad(1:nsec,1)-vrad(0:nsec-1,1)-vrad(1:nsec,2)-vrad(2:nsec+1,1)
vrad(1:nsec,nrad+1)=4.*vrad(1:nsec,nrad)-vrad(2:nsec+1,nrad)-vrad(0:nsec-1,nrad)-vrad(1:nsec,nrad-1)
;Finally setup array to hold radius of each grid point. 
;Also need values for radius of ghost cells above and below rmax. 
rad2=dblarr(nrad+2)
radius=dblarr(nsec+2,nrad+2)
rad2(1:nrad)=rmed(0:nrad-1)
rad2(0)=radtmp(0)-dr
rad2(nrad+1)=radtmp(nrad)+dr
;"radius" is 2d array to hold radius of each grid point. these do not
;vary with the first index (azimuth).
for l=0, nsec+1 do radius(l,0:nrad+1)=rad2(0:nrad+1)
;;;;;;;;;;;;;;;;;;;;;;;
;VORTICITY CALCULATION;
;;;;;;;;;;;;;;;;;;;;;;;
vort(0:nsec-1,0:nrad-1)=(radius(1:nsec,2:nrad+1)*vtheta(1:nsec,2:nrad+1)-radius(1:nsec,0:nrad-1)*vtheta(1:nsec,0:nrad-1))/(2.*dr)$
-(vrad(2:nsec+1,1:nrad)-vrad(0:nsec-1,1:nrad))/(2.*dtheta)
vort(0:nsec-1,0:nrad-1)=vort(0:nsec-1,0:nrad-1)/radius(1:nsec,1:nrad)
;;;;;;Below is just a different method to calculate vort using one loop.
; for j=1, nrad do begin
; vort(0:nsec-1,j-1)=(rad2(j+1)*vtheta(1:nsec,j+1)-rad2(j-1)*vtheta(1:nsec,j-1))/(2.*dr)$
; -(vrad(2:nsec+1,j)-vrad(0:nsec-1,j))/(2.*dtheta)
; vort(0:nsec-1,j-1)=vort(0:nsec-1,j-1)/rad2(j)
; endfor
;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;PLOT AZIMUTHAL-AVERAGED SIGMA/VORTICITY VS RADIUS;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
data = dblarr(nsec,nrad)
data1d=dblarr(nrad)
data = sigma/vort
plx=info(1,k)
ply=info(2,k)
if keyword_set(plazi) then begin
phi=pltphi(plx,ply)
array=abs(azi-phi)
grid=where(array eq min(array))
data1d(0:nrad-1)=data(grid,0:nrad-1)
endif else begin
for m=0, nrad-1 do data1d(m)=mean(data(*,m))
endelse
plotrange=[0.,max(data1d)]
time=string(k*dt/p0,format='(F5.1)')
set_plot, 'ps'
device, filename=filepath(strcompress('vort1d'+ks+'.ps',/remove_all),root_dir='.',subdir=[location])$
,bits_per_pixel=8,xsize=8, ysize=6,xoffset=0,yoffset=0,/inches
plot,rmed(1:nrad-2),data1d(1:nrad-2),title=time+' orbits',xmargin=[8,6],ymargin=[3,3],xtitle=xxtitle,ytitle=yytitle $
,charsize=1.5,xrange=[min(rmed),max(rmed)],xtickinterval=0.4,xminor=4,yrange=plotrange
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;OPTION TO PLOT PLANET RADIUS;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
if keyword_set(planet) then begin
pos=0.75*plotrange(1)
plrad=sqrt(plx*plx+ply*ply)
xyouts,plrad*(1.01),pos, 'planet', charsize=1.5, charthick=1.5
vline, plrad
endif
device,/close
print, 'done '+ks
endfor
end

