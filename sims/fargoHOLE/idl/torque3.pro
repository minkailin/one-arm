pro torque3, loc=loc, start=start, finish=finish, out=out, mp=mp, xtickinterval=xtickinterval $
            , scale=scale, xrange=xrange, yrange=yrange, plotrange0=plotrange0
if not keyword_set(finish) then finish=start
common consts, pi
pi=!dpi
if not keyword_set(mp) then mp=3d-4
f0=(mp/3.0)^(1.0/3.0)
h=0.05
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
period=2.*pi*(a0)^(3./2.)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;CREATE ARRAY FOR AZIMUTH AND RADIAL;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
azi=dblarr(nsec)
radtmp=dblarr(nrad+1)
rad=dblarr(nrad)
openr,1,filepath('used_rad.dat',root_dir='.',subdir=[location])
readf,1,radtmp
close,1
for i=0, nsec-1 do azi(i)=2.*pi*i/nsec
azi1 = azi/!dpi - 1.0
temp = min(abs(azi1),p0)

rad(0:nrad-1)=(radtmp(0:nrad-1)+radtmp(1:nrad))/2.0
dlogr = alog(radtmp(nrad)/radtmp(0))/nrad
dphi = 2.0*!dpi/nsec
;;;;;;;;;;;;;
;TORQUE CALC;
;;;;;;;;;;;;;
density = dblarr(nsec,nrad)
density_1d = dblarr(nrad)
torque_2d = dblarr(nsec,nrad)
torque_1d = dblarr(nrad)
dataplot = dblarr(nrad,nsec)
dataplot2 = dblarr(nrad,nsec)
torque_time = dblarr(2,finish-start+1)

for k=start, finish do begin
    ks=string(k,format='(I03)')
    openr,2,filepath(strcompress('gasdens'+string(k)+'.dat',/remove_all),root_dir='.',subdir=[location]) 
    readu,2,density
    close,2
    
    for i = 0, nrad -1 do density_1d(i) = mean(density(*,i))

    plx=info(1,k)
    ply=info(2,k)
    phi=pltphi(plx,ply)
    plrad = sqrt(plx*plx + ply*ply)
    temp = min(abs(rad-plrad),r0)
    
    rhill = 1.0*plrad*f0
    xs = 2.5*rhill
    temp = min(abs(rad-(plrad+xs)),rppxs)
    temp = min(abs(rad-(plrad-xs)),rpmxs)
    
    mhill = 0.0
    mhs = 0.0
    torque_acc = 0.0
    totalmass=0.0    
    for i=0, nsec-1 do begin
        for j=0, nrad-1 do begin
            dm = rad(j)*rad(j)*dphi*dlogr*density(i,j)
            totalmass += dm
            if((j le rppxs) and (j ge rpmxs)) then mhs += dm
            
            x = rad(j)*cos(azi(i))
            y = rad(j)*sin(azi(i))
            d = sqrt((x-plx)^2.0 + (y-ply)^2.0)
            
            if (d gt rhill) then begin
                dm /= d^3.0
                fx = (x-plx)*dm
                fy = (y-ply)*dm
                torque = plx*fy - ply*fx
                torque_acc += torque
                torque_2d(i,j) = torque
            endif else begin
                torque_2d(i,j) = 0.0
                mhill += dm
            endelse
        endfor
    endfor
    
    mpprime = mp + mhill
    meff = mp + mhs

print, 'total mass', totalmass

for i=0, nrad-1 do torque_1d(i) = mean(torque_2d(*,i))


 if(torque_acc gt 0.0) then begin
     temp = min(abs(azi1-(f0)),p1)
     sigma_s = mean(dataplot(r0:rppxs,p0:p1))
     sgn = 1.0
 endif
 
 if(torque_acc lt 0.0) then begin
     temp = min(abs(azi1-(-f0)),p1)
     sigma_s = mean(dataplot(rpmxs:r0,p1:p0))
     sgn = -1.0
 endif

 if(torque_acc eq 0.0) then sigma_s = 0.0

 torque_ratio = 4.0*!dpi*xs*plrad/meff
 torque_ratio *= sgn*sigma_s/sqrt(1.0+mpprime)

 torque_time(0,k-start) = info(7,k)/period
 torque_time(1,k-start) = torque_ratio

dataplot2 /= dataplot

pos = where(dataplot2 ge 0d0)
neg = where(dataplot2 le 0d0)

if not keyword_set(scale) then scale = 1d0
dataplot2 *= scale
;dataplot2 = alog(abs(dataplot2))

loadct,5, bottom=0

newazi = azi/pi - 1.0
temp = min(abs(newazi - yrange(0)), grid1)
temp = min(abs(newazi - yrange(1)), grid2)
temp = min(abs(rad - xrange(0)), grid3)
temp = min(abs(rad - xrange(1)), grid4)


print, 'done '+ks

endfor

openw,1,filepath(strcompress('torque_1d.dat',/remove_all),root_dir='.',subdir=[location])
for i = 0, nrad-1 do printf,1,rad(i),torque_1d(i)
close,1

openw,1,filepath(strcompress('torque_time.dat',/remove_all),root_dir='.',subdir=[location])
for i = 0, finish-start do printf,1,torque_time(0,i),torque_time(1,i)
close,1

set_plot, 'ps'
device, filename=filepath(strcompress('torque_time.ps',/remove_all),root_dir='.',subdir=[location])$
,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches
plot,torque_time(0,*),torque_time(1,*),xmargin=[8,4],ymargin=[3,3] $
,ytitle=textoidl('\Gamma_3/\Gamma_T'),xtitle='t/orbits' $
,charsize=1.5,xminor=4, thick=4
device,/close

end
