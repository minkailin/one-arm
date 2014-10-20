; THIS PROCEDURE CALCULATES SIGMA/OMEGA ON THE GRID POINTS USING
; CENTRAL DIFFERENCING, AND DOES A CARTESIAN CONTOUR PLOT CENTERED
; ABOUT PLANET.
pro vortxy, loc=loc, start=start, finish=finish, xs=xs
common consts, pi
pi=3.141592654
factor=double((0.00028/3.)^(1./3.))
if keyword_set(xs) then begin
factor=factor*xs
endif else factor=factor*2.5
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
;plot,
azifix=dblarr(nsec)
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
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;SETUP DONE, START CALCULATION AND PLOTTING;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
xxtitle=textoidl('r-r_p')
yytitle=textoidl('(\phi-\phi_p)/\pi')
loadct,39, bottom=0
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
;(theta,radius indicies). Start with boundary values (theta=0 and
;theta=2pi grid points). We on working on the data array of size
;(nsec,nrad), but calculating vorticity everywhere except first and
;last radius value. Then put results into the smaller array, vort.
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
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;NOW DO THE POLAR CONTOUR PLOTS;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
data=sigma(0:nsec-1,1:nrad-2)/vort(0:nsec-1,0:nrad-3)
;physical size of half width. find planet radius and azimuth, and find
;range of grid points that correspond to sector
plx=info(1,k)
ply=info(2,k)
plrad=sqrt(plx*plx+ply*ply)
width=factor*plrad
arr1=abs(rmed-(plrad-width))
arr2=abs(rmed-(plrad+width))
rmns=where(arr1 eq min(arr1))
rpls=where(arr2 eq min(arr2))
phi=pltphi(plx,ply)
; arr3=abs(azi-(phi-width/plrad))
; arr4=abs(azi-(phi+width/plrad))
; phimns=where(arr3 eq min(arr3))
; phipls=where(arr4 eq min(arr4))
data2=transpose(data(0:nsec-1,rmns:rpls))
plotrange=[alog(min(data2(where(data2 gt 0.)))),alog(max(data2))]
;plotrange=[-7.5,-4.]
;Deal with negative sigma/omega by assigning them small values:
zeros=where(data2 le 0.)
if min(zeros) ne -1 then data2(zeros)=exp(-10.)
data2=alog(data2)
levels=plotrange(0)+(plotrange(1)-plotrange(0))*(dindgen(32)/32.)
time=string(k*dt/p0,format='(F5.1)')
;The actual plotting part.
set_plot, 'ps'
device, filename=filepath(strcompress('vortxy'+ks+'.ps',/remove_all),root_dir='.',subdir=[location])$
,/color, bits_per_pixel=8,xsize=16, ysize=12
rplot=rmed-plrad
azi=azifix-phi
if phi gt pi then begin
excess=where(azi lt -pi)
azi(excess)=2.*pi+azi(excess)
seg1=where(azi eq min(azi))
if seg1 ne nsec-1 then begin
contour,data2(*,seg1:nsec-1),rplot(rmns:rpls),azi(seg1:nsec-1)/pi,/fill $
,levels=levels,title=time+' orbits',ymargin=[4,2],xmargin=[8,10],color=0$
,xtitle=xxtitle,ytitle=yytitle,charsize=1.5,yrange=[-1.,1.]
endif
if seg1 ne 1 then begin
contour,data2(*,0:seg1-1),rplot(rmns:rpls),azi(0:seg1-1)/pi,/fill $
,levels=levels,/overplot,yrange=[-1.,1.],title=time+' orbits',ymargin=[4,2],xmargin=[8,10],color=0$
,xtitle=xxtitle,ytitle=yytitle,charsize=1.5
endif
endif else if (0. le phi and  phi lt pi and (pi-phi) gt dtheta ) then begin
excess=where(azi ge pi)
azi(excess)=azi(excess)-2.*pi
seg2=where(azi eq max(azi))
if seg2 ne nsec-2 then begin
contour,data2(*,seg2+1:nsec-1),rplot(rmns:rpls),azi(seg2+1:nsec-1)/pi,/fill $
,levels=levels,title=time+' orbits',ymargin=[4,2],xmargin=[8,10],color=0$
,xtitle=xxtitle,ytitle=yytitle,charsize=1.5,yrange=[-1.,1.]
endif
if seg2 ne 0 then begin
contour,data2(*,0:seg2),rplot(rmns:rpls),azi(0:seg2)/pi,/fill $
,levels=levels,/overplot,yrange=[-1.,1.],title=time+' orbits',ymargin=[4,2],xmargin=[8,10],color=0$
,xtitle=xxtitle,ytitle=yytitle,charsize=1.5
endif
endif else begin
contour,data2,rplot(rmns:rpls),azi(0:nsec-1)/pi,/fill $
,levels=levels,title=time+' orbits',ymargin=[4,2],xmargin=[8,10],color=0$
,xtitle=xxtitle,ytitle=yytitle,charsize=1.5,yrange=[-1.,1.]
endelse
;do not plot points in order of index. may be gap so manually insert
;line segment
;arrow, azi(0),data1d(0),azi(nsec-1),data1d(nsec-1),hsize=0,/data
data3=dblarr(rpls-rmns+1,2)
data3(*,0)=data2(*,nsec-1)
data3(*,1)=data2(*,0)
contour,data3(*,0:1),rplot(rmns:rpls),[azi(nsec-1),azi(0)]/pi,/fill,levels=levels,/overplot
colorbar, position=[0.82, 0.18, 0.87, 0.88],/vertical,/right,range=plotrange,format='(F5.1)',color=0,charsize=1.5
oplot,[0.,0.],[0.,0.],psym=6,symsize=1,color=-1
device,/close
print, 'done'+ks
endfor
end
