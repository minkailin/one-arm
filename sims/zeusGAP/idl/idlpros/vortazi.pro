; THIS PROCEDURE CALCULATES SIGMA/OMEGA ON THE GRID POINTS USING
; CENTRAL DIFFERENCING. THIS VERSION PLOTS SIGMA/OMEGA AS WE GO AROUND
; A CIRCLE OF RADIUS RAD. DEFAULT IS TO NOT SET RAD, AND USE PLANET
; RADIUS. ALSO SET PLANET AZIMUTH=0, CENTER OF PLOT.
pro vortazi, loc=loc, start=start, finish=finish, rad=rad
common consts, pi
common arrays, data1d
pi=3.141592654
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;DATA LOCATION, GET INFO ABOUT RUN, READ PLANET ORBIT INFO. GET RADIUS;
;AND AZIMUTHAL ARRAY VALUES OF THE GRID.;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
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
azifix=dindgen(nsec)*2.*pi/nsec
;If rad is given, then use circle of radius rad. Find the grid closest to requested value.
if keyword_set(rad) then begin
rad=double(rad)
array=abs(rmed-rad)
grid=where(array eq min(array))
plpos=string(rad,format='(F4.2)')
endif
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;SETUP FOR VORTICITY CALCULATION;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;Arrays to hold data variables. Also cread 2d array to hold radii over grid.
sigma=dblarr(nsec,nrad)
vtheta= dblarr(nsec,nrad)
vrad=dblarr(nsec,nrad)
radius=dblarr(nsec,nrad)
for l=0, nsec-1 do radius(l,0:nrad-1)=rmed(0:nrad-1)
;Grid sizes
dtheta=2.*pi/nsec
dr=(radtmp(nrad)-radtmp(0))/nrad
;Array to hold results.Don't calculate vorticity for first and last
;radius (boundaries). 
vort=dblarr(nsec,nrad-2)
data1d=dblarr(nsec)
xxtitle=textoidl('\theta')
yytitle=textoidl('\Sigma/\omega')
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;SETUP DONE, START CALCULATION AND PLOTTING;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
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
;Now calculate vorticity. We need r*vtheta. For convenience define
;this array. Then begin finite differencing, Convention is (i,j) for
;(theta,radius indicies). Start with boundary values (theta=0 and theta=2pi grid points).
rvtheta=radius*vtheta
;When i=0, j=1 to nrad-2
vort(0,0:nrad-3)=(0.5*(rvtheta(0,2:nrad-1)+rvtheta(1,2:nrad-1))-0.5*(rvtheta(0,0:nrad-3)+rvtheta(1,0:nrad-3)))/(2.*dr) $
-(0.5*(vrad(1,1:nrad-2)+vrad(1,2:nrad-1))-0.5*(vrad(nsec-1,1:nrad-2)+vrad(nsec-1,2:nrad-1)))/(2.*dtheta)
;When i=nsec-1, j=1 to nrad-2
vort(nsec-1,0:nrad-3)=(0.5*(rvtheta(nsec-1,2:nrad-1)+rvtheta(0,2:nrad-1))-0.5*(rvtheta(nsec-1,0:nrad-3)+rvtheta(0,0:nrad-3)))/(2.*dr) $
-(0.5*(vrad(0,1:nrad-2)+vrad(0,2:nrad-1))-0.5*(vrad(nsec-2,1:nrad-2)+vrad(nsec-2,2:nrad-1)))/(2.*dtheta)
;Rest of grid, when i=1 to nsec-2 and j=1 to nrad-2
vort(1:nsec-2,0:nrad-3)=(0.5*(rvtheta(1:nsec-2,2:nrad-1)+rvtheta(2:nsec-1,2:nrad-1))-0.5*(rvtheta(1:nsec-2,0:nrad-3)+ $
rvtheta(2:nsec-1,0:nrad-3)))/(2.*dr) $
-(0.5*(vrad(2:nsec-1,1:nrad-2)+vrad(2:nsec-1,2:nrad-1))-0.5*(vrad(0:nsec-3,1:nrad-2)+vrad(0:nsec-3,2:nrad-1)))/(2.*dtheta)
;Still need to divide by radius at each grid to get vorticity,
vort(0:nsec-1,0:nrad-3)=vort(0:nsec-1,0:nrad-3)/radius(0:nsec-1,1:nrad-2)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;PLOT SIGMA/VORTICITY AROUND A CIRCLE OF RADIUS (DEFAULT: PLANET RADIUS);
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
data=sigma(0:nsec-1,1:nrad-2)/vort(0:nsec-1,0:nrad-3)
plx=info(1,k)
ply=info(2,k)
;Default is go around circle of planet radius. Later we can modify to
;radius of inner/outer edge of co-orbital region. This is roughly plus/minus planet Hill radius.
if not keyword_set(rad) then begin
plrad=sqrt(plx*plx+ply*ply)
array=abs(rmed-plrad)
grid=where(array eq min(array))
plpos=string(plrad,format='(F4.2)')
endif
data1d(0:nsec-1)=data(0:nsec-1,grid)
;plotrange=[min(data1d(where(data1d ge 0.))),max(data1d)]
plotrange=[0.,0.01]
time=string(k*dt/p0,format='(F5.1)')
set_plot, 'ps'
device, filename=filepath(strcompress('vortazi'+ks+'.ps',/remove_all),root_dir='.',subdir=[location])$
,bits_per_pixel=8,xsize=8, ysize=6,xoffset=0,yoffset=0,/inches,/color
; Now set planet to zero azimuth. This involves shifting azimuthal axis
; values and changing the order of plotting array values. If do not need to have
; planet at theta=0 then just use the last plot option (the else if).
phi=pltphi(plx,ply)
azi=azifix-phi
if phi gt pi then begin
excess=where(azi lt -pi)
azi(excess)=2.*pi+azi(excess)
seg1=where(azi eq min(azi))
plot,azi(seg1:nsec-1),data1d(seg1:nsec-1),title=time+' orbits, r='+plpos,xmargin=[9,5],ymargin=[3,3],xtitle=xxtitle $
,ytitle=yytitle,charsize=1.5,xrange=[-pi,pi],yrange=plotrange,xtickinterval=pi/4.,xminor=4
oplot,azi(0:seg1-1),data1d(0:seg1-1)
endif else if (0. le phi and  phi lt pi) then begin
excess=where(azi ge pi)
azi(excess)=azi(excess)-2.*pi
seg2=where(azi eq max(azi))
plot,azi(seg2+1:nsec-1),data1d(seg2+1:nsec-1),title=time+' orbits, r='+plpos,xmargin=[9,5],ymargin=[3,3],xtitle=xxtitle $
,ytitle=yytitle,charsize=1.5,xrange=[-pi,pi],yrange=plotrange,xtickinterval=pi/4.,xminor=4
oplot,azi(0:seg2),data1d(0:seg2)
endif else begin
plot,azi,data1d,title=time+' orbits, r='+plpos,xmargin=[9,5],ymargin=[3,3],xtitle=xxtitle $
,ytitle=yytitle,charsize=1.5,xrange=[-pi,pi],yrange=plotrange,xtickinterval=pi/4.,xminor=4
endelse
;Do not plot points in order of index. May be gap in plot so manually insert
;line segment. Usually unimportant.
arrow, azi(0),data1d(0),azi(nsec-1),data1d(nsec-1),hsize=0,/data
vline, 0.
device,/close
print, 'done '+ks
endfor
end

