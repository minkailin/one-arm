pro polar, type=type, loc=loc, start=start, finish=finish $
,ct=ct,plotrange0=plotrange0, log=log, xrange=xrange, nopert=nopert, basic=basic $
,out=out, yrange=yrange, xtickinterval=xtickinterval
if not keyword_set(ct) then ct=5
if not keyword_set(finish) then finish=start
common consts, pi, nrad, time
pi=3.141592654
;;;;;;;;;;;;;;;
;WHERE IS DATA;
;;;;;;;;;;;;;;;
location=strcompress(loc,/remove_all)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;GET DIMENSIONS AND TIME UNIT;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
dims=(read_ascii(filepath('dims.dat',root_dir='.',subdir=[location]))).(0)
nout=fix(dims(5))
nrad=fix(dims(6))
nsec=fix(dims(7))
if not keyword_set (out) then begin
info=dblarr(11,nout+1)
endif else info=dblarr(11,out+1)

;info=dblarr(11,86)
openr,3,filepath('planet0.dat',root_dir='.',subdir=[location])
readf,3,info
close,3
a0=info(1,0)
dt=info(7,1)
p0=2.*pi*(a0)^(3./2.)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;CREATE ARRAY FOR AZIMUTH AND RADIAL;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
azi=dblarr(nsec)
azi1=dblarr(nsec)
radtmp=dblarr(nrad+1)
rad=dblarr(nrad)
openr,1,filepath('used_rad.dat',root_dir='.',subdir=[location])
readf,1,radtmp
close,1
for i=0, nsec-1 do azi(i)=2.*pi*i/nsec
rad(0:nrad-1)=(radtmp(0:nrad-1)+radtmp(1:nrad))/2.
;;;;;;;;;;;;;;;;;;;;;;;;
;DO POLAR CONTOUR PLOTS;
;;;;;;;;;;;;;;;;;;;;;;;;
data = dblarr(nsec,nrad)
data2 = dblarr(nsec,nrad)
data0=dblarr(nsec,nrad)

loadct,ct, bottom=0
if not keyword_set(nopert) then begin
if not keyword_set(basic) then basic='20'
openr,2,filepath(strcompress('gas'+type+basic+'.dat',/remove_all),root_dir='.',subdir=[location])
readu,2,data0
close,2
endif

for k=start, finish do begin
ks=string(k,format='(I03)')
openr,2,filepath(strcompress('gas'+type+string(k)+'.dat',/remove_all),root_dir='.',subdir=[location]) 
readu,2,data
close,2

;openr,2,filepath(strcompress('gasdens'+string(k)+'.dat',/remove_all),root_dir='.',subdir=[location]) 
;readu,2,data2
;close,2

;data *=data2

;for i=0, nsec-1 do data(i,0:nrad-1) *= rad(0:nrad-1)

if not keyword_set(nopert) then begin

data /= data0
if not keyword_set(log) then data -= 1.0

endif

if keyword_set(log) then data=alog10(data)
;data2(*,0:in-1)=0.
;data2(*,in:nrad+in-1)=data
if not keyword_set(plotrange0) then begin 
plotrange=[min(data),max(data)]
endif else plotrange=plotrange0
levels=plotrange(0)+(plotrange(1)-plotrange(0))*(dindgen(48)/47.)
time=string(info(7,k)/p0,format='(F7.2)')
plx=info(1,k)
ply=info(2,k)

;print, 'plx, ply', plx, ply

plrad=sqrt(plx*plx+ply*ply)
phi=pltphi(plx,ply)
temp=min(abs(azi-phi),grid)
azi1(0:nsec-1)=azi(0:nsec-1)-phi*(1d0 + 0.1);-azi(grid)
set_plot, 'ps'
data2=dblarr(nsec,nrad+2)
radnew=dblarr(nrad+2)
radnew(2:nrad+1) = rad(0:nrad-1)
radnew(0:1) = [0.0,rad(0)/2.0]
data2(0:nsec-1,2:nrad+1) = data(0:nsec-1,0:nrad-1)
data2(0:nsec-1,0:1) = 10.0
device, filename=filepath(strcompress('polar_'+type+ks+'.ps',/remove_all),root_dir='.',subdir=[location])$
,/color, bits_per_pixel=8,xsize=14, ysize=12
;temp=min(abs(radnew-10.0), grid)
;data2(*,grid:nrad-1) = abs(max(data2))*2.0

data2 =congrid(data2, nsec/2, nrad/2,/interp)
azi1 = congrid(azi1,nsec/2,/interp)
radnew=congrid(radnew,nrad/2,/interp)

polar_contour,data2,azi1,radnew,/isotropic,/fill,levels=levels,title=time+' orbits',ymargin=[2,2],xmargin=[4,4],xrange=xrange,yrange=yrange, xtickinterval=xtickinterval, ytickinterval=xtickinterval
colorbar, position=[0.85, 0.07, 0.9, 0.93],/vertical,/right,range=plotrange,format='(f5.2)'
oplot,[plrad,plrad],[0.,0.],psym=6,symsize=1,color=120
device,/close
print, 'done '+ks
endfor
end
