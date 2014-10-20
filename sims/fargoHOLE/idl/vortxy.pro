; THIS PROCEDURE CALCULATES SIGMA/OMEGA ON THE GRID POINTS USING
; CENTRAL DIFFERENCING, AND DOES A POLAR CONTOUR PLOT. OPTION TO
; ROTATE DISC SUCH THAT PLANET IS AT AZIMUTH=0.
; 20/12/2008: Added option to overplot tracer fluid density. But currently disabled (commented out).
pro vortxy, loc=loc, start=start, finish=finish,plotrange0=plotrange0,ct=ct,out=out,track=track,noplot=noplot $
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
dr=(radtmp(nrad)-radtmp(0))/nrad
for i=0, nsec-1 do dr2d(i,*)=dr(*)
;Array to hold results.Don't calculate vorticity for first and last
;radius (boundaries).
vort=dblarr(nsec,nrad-2)
dataplot=dblarr(nrad-2,nsec)
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
;data = 1d0/data

;zeros = where( data le 0.0)
;data(zeros) = 1d-10
;data = -alog10(data)

;data=1.0/data
;data /= 1d3

;data=vort(0:nsec-1,0:nrad-3)
;data=alog10(data)
if not keyword_set(plotrange0) then begin

temp=min(abs(xrange(0)-rmed(1:nrad-2)),inner)
temp=min(abs(xrange(1)-rmed(1:nrad-2)),outer)


plotrange=[min(data(*,inner:outer)),max(data(*,inner:outer))]
endif else plotrange=plotrange0

;print, plotrange

levels=plotrange(0)+(plotrange(1)-plotrange(0))*(dindgen(32)/31.)
time=string(k*dt/p0,format='(F7.2)')

data = transpose(data)


if keyword_set(track) then begin
    data(0:subset(0)-1,*) = 10.0
    data(subset(n_elements(subset)-1):nrad-3,*) = 10.0
    min_eta = min(abs(data))
    vortices = where((abs(data) ge min_eta) and (abs(data) le min_eta+abs(min_eta)*2.75))
    ;vortices = where(data lt 0.0)
    points = n_elements(vortices)
    position=dblarr(2,points)
;scale_height = dblarr(points*(points-1)/2)
;beg = 0

    scale_height = dblarr(points,points)
    
    for l = 0, points-1 do begin
        result=array_indices(data, vortices(l))
        r = rplot(result(0))
        p = azi(result(1))
        position(0,l) = r*cos(p)
        position(1,l) = r*sin(p)
        scale_height(*,l) = r*0.05
;if(l lt points-1) then begin
;nl = points - (l+1)
;scale_height(beg: nl+beg-1) = r*0.05
;beg = nl+beg
;endif
    endfor
    
    result = distance_measure(position,/double,/matrix)
    same_vortex = where(result le 3.0*scale_height)
    result(same_vortex) = 10.0
    
    temp = min(result,min_dist)
    
    result2 = array_indices(result, min_dist)
    
    vor1 = result2(0)
    vor2 = result2(1)
    
    merge(0,k-start) = info(7,k)/p0
    merge(1,k-start) = min(result)
endif




;The actual plotting part.


if not keyword_set(noplot) then begin


    

    time=string(info(7,k)/p0,format='(F7.2)')
    plx=info(1,k)
    ply=info(2,k)
    plrad=sqrt(plx*plx+ply*ply)
    phi=pltphi(plx,ply)

    rhill = f0*plrad
    ;rmed = (rmed - plrad)/rhill

    if not keyword_set(xrange) then xrange=[min(rmed),max(rmed)]
    rplot=rmed(1:nrad-2)
    subset = where((rplot ge xrange(0)) and (rplot le xrange(1)))
    subset = where((rplot ge xrange(0)) and (rplot le xrange(1)))

    if(phi gt !dpi) then begin
        temp2=min(abs(azi-(phi-!dpi)),grid2)
        dataplot(0:nrad-3,0:nsec-1-grid2) = data(0:nrad-3,grid2:nsec-1)
        dataplot(0:nrad-3,nsec-grid2:nsec-1) = data(0:nrad-3, 0:grid2-1) 
    endif
    
    if(phi lt !dpi) then begin
        temp2=min(abs(azi-(phi+!dpi)),grid2)
        dataplot(0:nrad-3,nsec-grid2:nsec-1) = data(0:nrad-3,0:grid2-1)
        dataplot(0:nrad-3,0:nsec-1-grid2) = data(0:nrad-3, grid2:nsec-1) 
    endif

if(phi eq !dpi) then dataplot = data




set_plot, 'ps'

if not keyword_set(plotrange0) then begin

    plotrange=[min(data(subset(0):subset(n_elements(subset)-2),*)),max(data(subset(0)+1:subset(n_elements(subset)-2),*))]
    
endif else plotrange=plotrange0
levels=plotrange(0)+(plotrange(1)-plotrange(0))*(dindgen(32)/31.)

device, filename=filepath(strcompress('vortxy_'+ks+'.ps',/remove_all),root_dir='.',subdir=[location])$
,/color, bits_per_pixel=8,xsize=12, ysize=14
contour,dataplot,rmed(1:nrad-2),azi/pi - 1.0,levels=levels,title=time+' orbits', $
  xtitle='r', ytitle=textoidl('(\phi-\phi_p)/\pi'),charsize=1.5,xmargin=[7.0,6.0], ymargin=[3,2], $
  xtickinterval=xtickinterval,xrange=xrange,yrange=yrange,/fill
colorbar, position=[0.85, 0.15, 0.90, 0.90],/vertical,/right,range=plotrange,format='(f5.2)'
oplot,[0.0,0.0],[0.0,0.0]/pi,psym=7,symsize=1.5,color=!D.Table_size*0.0


;levels = dblarr(16)
;levels(0:2) = [-0.36, 0.65, 1.66]
;levels(3:15) = 2.67 + (4.02-2.67)*dindgen(13)/12.
;contour,dataplot,rmed(1:nrad-2),azi/pi - 1.0 ,levels=levels,/overplot,c_color=0

;xyouts,plrad,phi/pi,'X',charsize=1.5


; for l=0, points-1 do begin
;   x1 = position(0,l)
;   y1 = position(1,l)
;   r1 = sqrt(x1*x1+y1*y1)
;   p1 = pltphi(x1,y1)
;   oplot,[r1,r1],[p1,p1]/pi, psym=7, color=!D.Table_size*0.5

; endfor

; zeros = where(data le 0.0)
; for l = 0, n_elements(zeros)-1  do begin
; result = array_indices(data,zeros(l))
; xx = rplot(result(0))
; yy = (azi(result(1)) - phi)/pi
; r1 = sqrt(xx*xx+yy*yy)
; p1 = pltphi(xx,yy)
; oplot,[r1,r1],[p1,p1], psym=7, color=!D.Table_size*0.0
; endfor


if keyword_set(track) then begin
    
    x1 = position(0,vor1)
    y1 = position(1,vor1)
    r1 = sqrt(x1*x1+y1*y1)
    p1 = pltphi(x1,y1) - phi 
    oplot,[r1,r1],[0,0]/pi, psym=7, color=!D.Table_size*0.0
    
    x1 = position(0,vor2)
    y1 = position(1,vor2)
    r1 = sqrt(x1*x1+y1*y1)
    p1 = pltphi(x1,y1) - phi
    oplot,[r1,r1],[p1,p1]/pi, psym=7, color=!D.Table_size*0.0




endif
device,/close
endif




print, 'done '+ks

endfor

if keyword_set(track) then begin
    set_plot, 'ps'
    device, filename=filepath(strcompress('vortxy_merge.ps',/remove_all),root_dir='.',subdir=[location])$
      ,bits_per_pixel=8,xsize=8, ysize=6,xoffset=0,yoffset=0,/inches
    plot,merge(0,*),merge(1,*),xmargin=[8,4],ymargin=[3,3] $
      ,ytitle=textoidl('minimum inter-vortex distance'),xtitle='t' $
      ,charsize=1.5,xminor=4, thick=4
    device,/close
    
    openw,1,filepath(strcompress('vortxy_merge.dat',/remove_all),root_dir='.',subdir=[location])
    for i = 0, finish-start do printf,1,merge(0,i),merge(1,i)
    close,1
endif


end
