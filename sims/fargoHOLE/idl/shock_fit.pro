; THIS PROCEDURE CALCULATES |GRAD SIGMA|/SIGMA ON THE GRID POINTS USING
; CENTRAL DIFFERENCING, AND DOES A CARTESIAN CONTOUR PLOT CENTERED
; ABOUT PLANET.
pro fitting, eta=eta, yrange=yrange, order=order
common consts, pi
common grid_in, rmed, azinew, rmns, rpls
common info, plrad, data3 
common curve, coeff
;Convert to shearing box (x,y) units
x=rmed(rmns:rpls)-plrad
arr1=min(abs(azinew-yrange(0)),phimns)
arr2=min(abs(azinew-yrange(1)),phipls)
y=plrad*azinew(phimns:phipls)*pi
;Shock location
data4=data3(0:rpls-rmns,phimns:phipls)
max_grad=max(data4)
shock_loc=where(data4 ge eta*max_grad)
ind=array_indices(data4,shock_loc)
s=size(ind)
print, 'number of shock points=',s(2)
xss=dblarr(s(2))
yss=dblarr(s(2))
for i=0, s(2)-1 do begin
xss(i)=x(ind(0,i))
yss(i)=y(ind(1,i))
endfor
;Fit the curve
coeff=svdfit(xss,yss,order+1,/double)
end
pro shock_fit, loc=loc, start=start, finish=finish, xs=xs, tracer=tracer $
,plotrange=plotrange, xrange=xrange, ct=ct,yrange=yrange, eta=eta, order=order $
,noplot=noplot, out=out, mass=mass
common grid_in, rmed, azinew, rmns, rpls
common info, plrad, data3
common consts, pi
common curve, coeff
pi=3.141592654
factor0=double((mass/3.)^(1./3.))
if not keyword_set(ct) then ct=5
if not keyword_set(yrange) then yrange=[-0.1,0.1]
if keyword_set(plotrange) then plotrange0=plotrange
if n_elements(start) ne 1 then begin
cases=start
endif
if not keyword_set(finish) then begin
cases=[start]
endif else cases=start+indgen(finish-start+1)
numcases=n_elements(cases)
if not keyword_set(xs) then xs=2.0
if not keyword_set(xrange) then xrange=[-xs,xs]
if not keyword_set(eta) then eta=0.95
if not keyword_set(order) then order=2
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
radtmp=dblarr(nrad+1)
openr,1,filepath('used_rad.dat',root_dir='.',subdir=[location])
readf,1,radtmp
close,1
azifix=dindgen(nsec)*2.*pi/nsec
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;SETUP FOR GRADIENT CALCULATION;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;Arrays to hold variables.Also create 2d array to hold radii values
;over grid.
sigma=dblarr(nsec,nrad)
label=dblarr(nsec,nrad)
;vtheta= dblarr(nsec,nrad)
;vrad=dblarr(nsec,nrad)
rmed=radtmp(0:nrad-1)+radtmp(1:nrad)
rmed=rmed*0.5
radius=dblarr(nsec,nrad)
for l=0, nsec-1 do radius(l,0:nrad-1)=rmed(0:nrad-1)
;Cell sizes.
dtheta=2.*pi/nsec
dr=(radtmp(nrad)-radtmp(0))/nrad
;Array to hold results.Don't calculate for first and last
;radius (boundaries).
grad=dblarr(nsec,nrad-2)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;SETUP DONE, START CALCULATION AND PLOTTING;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
xxtitle=textoidl('(r-r_p)/d')
yytitle=textoidl('(\phi-\phi_p)/\pi')
loadct,ct, bottom=0
for k=0, numcases-1 do begin
ks=string(cases(k),format='(I03)')
;Read raw data.
openr,2,filepath(strcompress('gasdens'+string(cases(k))+'.dat',/remove_all),root_dir='.',subdir=[location]) 
readu,2,sigma
close,2
;openr,3,filepath(strcompress('gasvtheta'+string(cases(k))+'.dat',/remove_all),root_dir='.',subdir=[location]) 
;readu,3,vtheta
;close,3
;openr,4,filepath(strcompress('gasvrad'+string(cases(k))+'.dat',/remove_all),root_dir='.',subdir=[location]) 
;readu,4,vrad
;close,4
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;Now calculate grad sigma. Convention is (i,j) for
;(theta,radius indicies). Start with boundary values (theta=0 and
;theta=2pi grid points). We on working on the data array of size
;(nsec,nrad), but calculating grad everywhere except first and
;last radius value. Then put results into the smaller array, grad.
;When i=0, j=1 to nrad-2
dr=radius(0,2:nrad-1)-radius(0,0:nrad-3)
grad(0,0:nrad-3)=((sigma(1,1:nrad-2)-sigma(nsec-1,1:nrad-2))/(2.0*dtheta*radius(0,1:nrad-2)))^2.0 $
+((sigma(0,2:nrad-1)-sigma(0,0:nrad-3))/(2.0*dr))^2.0
;When i=nsec-1, j=1 to nrad-2
dr=radius(0,2:nrad-1)-radius(0,0:nrad-3)
grad(nsec-1,0:nrad-3)=((sigma(0,1:nrad-2)-sigma(nsec-2,1:nrad-2))/(2.0*dtheta*radius(nsec-1,1:nrad-2)))^2.0 $
+((sigma(nsec-1,2:nrad-1)-sigma(nsec-1,0:nrad-3))/(2.0*dr))^2.0
;Rest of grid, when i=1 to nsec-2 and j=1 to nrad-2
dr=radius(1:nsec-2,2:nrad-1)-radius(1:nsec-2,0:nrad-3)
grad(1:nsec-2,0:nrad-3)=((sigma(2:nsec-1,1:nrad-2)-sigma(0:nsec-3,1:nrad-2))/(2.0*dtheta*radius(1:nsec-2,1:nrad-2)))^2.0 $
+((sigma(1:nsec-2,2:nrad-1)-sigma(1:nsec-2,0:nrad-3))/(2.0*dr))^2.0
;Take square root to get modulus of grad sigma vector, then divide by sigma
grad(0:nsec-1,0:nrad-3)=sqrt(grad(0:nsec-1,0:nrad-3))/sigma(0:nsec-1,1:nrad-2)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;NOW DO THE POLAR CONTOUR PLOTS;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
data=grad
;physical size of half width. find planet radius and azimuth, and find
;range of grid points that correspond to sector
plx=info(1,cases(k))
ply=info(2,cases(k))
plrad=sqrt(plx*plx+ply*ply)
phi=pltphi(plx,ply)
rplot=(rmed-plrad)/(factor0*plrad)
azi=azifix-phi
arr1=min(abs(rplot-xrange(0)),rmns)
arr2=min(abs(rplot-xrange(1)),rpls)
data2=transpose(data(0:nsec-1,rmns-1:rpls-1))
data2=alog(data2)
if not keyword_set(plotrange0) then begin
plotrange=[min(data2),max(data2)]
endif else if plotrange0(0) eq 0 then begin 
plotrange=[min(data2),plotrange(1)]
endif
levels=plotrange(0)+(plotrange(1)-plotrange(0))*(dindgen(128)/128.)
time=string(info(7,cases(k))/p0,format='(F6.2)')
;Shift azimuthal index
data3=dblarr(rpls-rmns+1,nsec)
azinew=dblarr(nsec)
if phi gt pi then begin
excess=where(azi lt -pi)
azi(excess)=2.*pi+azi(excess)
seg1=where(azi eq min(azi))
if seg1 ne nsec-1 then begin
data3(*,0:nsec-seg1-1)=data2(*,seg1:nsec-1)
azinew(0:nsec-seg1-1)=azi(seg1:nsec-1)/pi
endif
if seg1 ne 1 then begin
data3(*,nsec-seg1:nsec-1)=data2(*,0:seg1-1)
azinew(nsec-seg1:nsec-1)=azi(0:seg1-1)/pi
endif
endif else if (0. le phi and  phi lt pi and (pi-phi) gt dtheta ) then begin
excess=where(azi ge pi)
azi(excess)=azi(excess)-2.*pi
seg2=where(azi eq max(azi))
if seg2 ne nsec-2 then begin
data3(*,0:nsec-seg2-2)=data2(*,seg2+1:nsec-1)
azinew(0:nsec-seg2-2)=azi(seg2+1:nsec-1)/pi
endif
if seg2 ne 0 then begin
data3(*,nsec-seg2-1:nsec-1)=data2(*,0:seg2)
azinew(nsec-seg2-1:nsec-1)=azi(0:seg2)/pi
endif
endif else begin
data3=data2
azinew(0:nsec-1)=azi(0:nsec-1)/pi
endelse
;Shock curve
fitting, eta=eta, yrange=yrange, order=order
if not keyword_set(noplot) then begin
ys=0.0
for i=0, order do begin
ys=ys+coeff(i)*(rplot*plrad*factor0)^i
endfor
ys=ys/(pi*plrad)
;Plot
set_plot, 'ps'
device, filename=filepath(strcompress('shock_fit_'+ks+'.ps',/remove_all),root_dir='.',subdir=[location])$
,/color, bits_per_pixel=8,xsize=14, ysize=12
temp=strcompress(string(plrad,format='(f5.2)'),/remove_all)
contour,data3,rplot(rmns:rpls),azinew(0:nsec-1) ,/fill $
,levels=levels,title=time+' orbits,'+textoidl(' r_p=')+temp,ymargin=[4,2],xmargin=[8,10],color=0$
,xtitle=xxtitle,ytitle=yytitle,charsize=1.5,yrange=yrange,xrange=xrange $
,xtickinterval=2.0
; for j=0, s(2)-1 do begin
; oplot, [rplot(ind(0,j)+rmns),rplot(ind(0,j)+rmns)] $
; ,[azinew(ind(1,j)+phimns),azinew(ind(1,j)+phimns)], color=0, psym=1
; endfor
oplot, rplot, ys, thick=4
colorbar, position=[0.80, 0.18, 0.85, 0.88],/vertical,/right,range=plotrange,format='(F5.1)',color=0,charsize=1.5
device,/close
endif
print, 'done'+ks
endfor
end


