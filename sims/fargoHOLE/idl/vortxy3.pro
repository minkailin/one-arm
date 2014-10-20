; THIS PROCEDURE CALCULATES SIGMA/OMEGA ON THE GRID POINTS USING
; CENTRAL DIFFERENCING, AND DOES A POLAR CONTOUR PLOT. OPTION TO
; ROTATE DISC SUCH THAT PLANET IS AT AZIMUTH=0.
; 20/12/2008: Added option to overplot tracer fluid density. But currently disabled (commented out).
pro vortxy3, loc=loc, start=start, finish=finish,plotrange0=plotrange0,ct=ct,out=out,track=track,noplot=noplot $
            ,xrange=xrange,xtickinterval=xtickinterval,yrange=yrange, mp=mp
common consts, pi
if not keyword_set(mp) then mp=3d-4
f0=(mp/3.0)^(1.0/3.0)

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
rmed=radtmp(0:nrad-1)+radtmp(1:nrad)
rmed=rmed*0.5
;rmed = 2.0/3.0*(radtmp(1:nrad)^3.-radtmp(0:nrad-1)^3.);
;rmed(0:nrad-1) = rmed(0:nrad-1) / (radtmp(1:nrad)^2.-radtmp(0:nrad-1)^2.)

for j=in, nrad+in-1 do rad(j)=rmed(j-in)
radius=dblarr(nsec,nrad)
for l=0, nsec-1 do radius(l,0:nrad-1)=rmed(0:nrad-1)
;Cell sizes.
dtheta=2.*pi/nsec
dlogr=alog(radtmp(nrad)/radtmp(0))/nrad


;Array to hold results.Don't calculate vorticity for first and last
;radius (boundaries).
vort=dblarr(nsec,nrad-2)

data2=dblarr(nsec,nrad-2+in+1)
aziscale=azi/pi
merge=dblarr(2,finish-start+1)
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
data = vort(0:nsec-1,0:nrad-3);/sigma(0:nsec-1,1:nrad-2)


; data = transpose(data)
; sigma= transpose(sigma)
; vtheta=transpose(vtheta)
; vrad  =transpose(vrad)



; plx=info(1,k)
; ply=info(2,k)
; plrad=sqrt(plx*plx+ply*ply)
; phi=pltphi(plx,ply)
; if(phi gt !dpi) then begin
;     temp2=min(abs(azi-(phi-!dpi)),grid2)
;     dataplot(0:nrad-3,0:nsec-1-grid2) = data(0:nrad-3,grid2:nsec-1)
;     dataplot(0:nrad-3,nsec-grid2:nsec-1) = data(0:nrad-3, 0:grid2-1) 

;     sigmaT(*,0:nsec-1-grid2) = sigma(*,grid2:nsec-1)
;     sigmaT(*,nsec-grid2:nsec-1) = sigma(*, 0:grid2-1)

;     vthetaT(*,0:nsec-1-grid2) = vtheta(*,grid2:nsec-1)
;     vthetaT(*,nsec-grid2:nsec-1) = vtheta(*, 0:grid2-1)
    
;     vradT(*,0:nsec-1-grid2) = vrad(*,grid2:nsec-1)
;     vradT(*,nsec-grid2:nsec-1) = vrad(*, 0:grid2-1)
; endif

; if(phi lt !dpi) then begin
;     temp2=min(abs(azi-(phi+!dpi)),grid2)
;     dataplot(0:nrad-3,nsec-grid2:nsec-1) = data(0:nrad-3,0:grid2-1)
;     dataplot(0:nrad-3,0:nsec-1-grid2) = data(0:nrad-3, grid2:nsec-1) 

;     sigmaT(*,nsec-grid2:nsec-1) = sigma(*,0:grid2-1)
;     sigmaT(*,0:nsec-1-grid2) = sigma(*, grid2:nsec-1)

;     vthetaT(*,nsec-grid2:nsec-1) = vtheta(*,0:grid2-1)
;     vthetaT(*,0:nsec-1-grid2) = vtheta(*, grid2:nsec-1)
    
;     vradT(*,nsec-grid2:nsec-1) = vrad(*,0:grid2-1)
;     vradT(*,0:nsec-1-grid2) = vrad(*, grid2:nsec-1)
; endif

; if(phi eq !dpi) then begin
;     dataplot = data
;     sigmaT = sigma
;     vthetaT = vtheta
;     vradT = vrad
; endif


temp = min(abs(rmed - xrange(0)), x1)
temp = min(abs(rmed - xrange(1)), x2)
;temp = min(abs(azi/!dpi - 1d0 - yrange(0)), y1)
;temp = min(abs(azi/!dpi - 1d0 - yrange(1)), y2)

sub_sig = sigma(*,x1:x2)
sub_da = sub_sig
sub_vr =  vrad(*,x1:x2)
sub_vt = vtheta(*,x1:x2)

for i=0, nsec-1 do sub_da(i,*) = rmed(x1:x2)^2*dlogr*dtheta  

rplot = rmed(1:nrad-2)
azi = 2d0*!dpi*dindgen(nsec)/nsec

rplot2d = dblarr(nsec,nrad-2)
azi2d = dblarr(nsec, nrad-2)

for i=0, nsec-1 do begin
    rplot2d(i,*) = rplot(*)
    azi2d(i,*) = 2d0*!dpi*i/nsec
endfor

temp = min(abs(rplot - xrange(0)), x1)
temp = min(abs(rplot - xrange(1)), x2)
sub_vort = data(*,x1:x2)
sub_rplot2d = rplot2d(*,x1:x2)
sub_azi2d = azi2d(*,x1:x2)

temp = min(sub_vort, cent)

rcen = sub_rplot2d(cent)
tcen = sub_azi2d(cent) 

sub_rplot2d -= rcen
sub_azi2d  = rcen*(sub_azi2d - tcen)

dist = sqrt(sub_rplot2d*sub_rplot2d + sub_azi2d*sub_azi2d)

vcen = sub_vt(cent)

sub_vt -= vcen
sub_vr -= sub_vr(cent)

vsq = sub_vr^2 + sub_vt^2

filter = where( (sub_vort lt 0.0035))

print, 'vortex mass =', total(sub_sig(filter)*sub_da(filter))
print, 'vortex area= ', total(sub_da(filter))

print, 'average vsq is', mean(vsq(filter))
print, 'vortex radius', sqrt(total(sub_da(filter))/!dpi)
print, 'mean density', mean(sub_sig(filter))


sub_vort(where(sub_vort ge 0.0035)) = 0

contour, sub_vort, levels=levels

endfor

end
