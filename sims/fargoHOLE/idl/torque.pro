pro torque, loc=loc, start=start, finish=finish, out=out, mp=mp;, xtickinterval=xtickinterval $
            ;, scale=scale, xrange=xrange, yrange=yrange, plotrange0=plotrange0
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
azi1 = azi/pi - 1d0
temp=min(abs(azi1), p0)

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

torque_time = dblarr(4,finish-start+1)

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
    totalmass= 0.0
;    mhs = 0.0
    torque_inc = 0.0
    torque_exl = 0.0  
    for i=0, nsec-1 do begin
        for j=0, nrad-1 do begin
            dm = rad(j)*rad(j)*dphi*dlogr*density(i,j)
            totalmass += dm
           
            x = rad(j)*cos(azi(i))
            y = rad(j)*sin(azi(i))
            d = sqrt((x-plx)^2.0 + (y-ply)^2.0)
            
                dm2 = dm/d^3.0
                fx = (x-plx)*dm2
                fy = (y-ply)*dm2
                torque = plx*fy - ply*fx ;torque per planet mass
                                ;torque_inc += torque 
                torque_2d(i,j) = torque
                
                                ;if (d gt rhill) then begin
                torque_exl += torque

		if (d le 1.0*rhill) then begin 
                	torque_2d(i,j) = 0.0	
                	mhill += dm  
                endif
        endfor
    endfor
    
;    mpprime = mp + mhill
;    meff = mp + mhs

print, 'total mass', totalmass

for i=0, nrad-1 do begin
;    nonzero = where(torque_2d(*,i) ne 0.0)
    torque_1d(i) = total(torque_2d(*,i))/(rad(i)*dlogr)
endfor


; re-shuffle grid
    
; dataT = transpose(density)
; dataT2 = transpose(torque_2d)
    
;  if(phi gt !dpi) then begin
;      temp2=min(abs(azi-(phi-!dpi)),grid2)
;      dataplot(0:nrad-1,0:nsec-1-grid2) = dataT(0:nrad-1,grid2:nsec-1)
;      dataplot(0:nrad-1,nsec-grid2:nsec-1) = dataT(0:nrad-1, 0:grid2-1)

;      dataplot2(0:nrad-1,0:nsec-1-grid2) = dataT2(0:nrad-1,grid2:nsec-1)
;      dataplot2(0:nrad-1,nsec-grid2:nsec-1) = dataT2(0:nrad-1, 0:grid2-1)
;  endif
 
;  if(phi lt !dpi) then begin
;      temp2=min(abs(azi-(phi+!dpi)),grid2)
;      dataplot(0:nrad-1,nsec-grid2:nsec-1) = dataT(0:nrad-1,0:grid2-1)
;      dataplot(0:nrad-1,0:nsec-1-grid2) = dataT(0:nrad-1, grid2:nsec-1)
     
;      dataplot2(0:nrad-1,nsec-grid2:nsec-1) = dataT2(0:nrad-1,0:grid2-1)
;      dataplot2(0:nrad-1,0:nsec-1-grid2) = dataT2(0:nrad-1, grid2:nsec-1)
;  endif
 
;  if(phi eq !dpi) then begin
;      dataplot = dataT
;      dataplot2 = dataT2
;  endif

 if(torque_inc gt 0.0) then begin
;      temp = min(abs(azi1-(f0)),p1)
;      sigma_s = mean(dataplot(r0:rppxs,p0:p1))

     sigma_s = density_1d(rppxs)
     sgn = 1.0
 endif
 if(torque_inc lt 0.0) then begin
;     temp = min(abs(azi1-(-f0)),p1)
;     sigma_s = mean(dataplot(rpmxs:r0,p1:p0))
     sigma_s = density_1d(rpmxs)
     sgn = -1.0
 endif
 if(torque_inc eq 0.0) then begin
     sigma_s = 0.0
     sgn = 0.0
 endif

 torque_ratio_inc = 4.0*!dpi*xs*plrad*sigma_s*sgn/mp


 if(torque_exl gt 0.0) then begin
     sigma_s = density_1d(rppxs)
     sgn = 1.0
 endif
 if(torque_exl lt 0.0) then begin
     sigma_s = density_1d(rpmxs)
     sgn = -1.0
 endif
 if(torque_exl eq 0.0) then sigma_s = 0.0
 torque_ratio_exl = 4.0*!dpi*xs*plrad*sigma_s*sgn/(mp+mhill)


; torque_ratio = 4.0*!dpi*xs*plrad/meff
; torque_ratio *= sgn*sigma_s/sqrt(1.0+mpprime)

 torque_time(0,k-start) = info(7,k)/period
 torque_time(1,k-start) = torque_ratio_inc
 torque_time(2,k-start) = torque_ratio_exl
 torque_time(3,k-start) = total(torque_2d)

;dataplot2 /= dataplot

;pos = where(dataplot2 ge 0d0)
;neg = where(dataplot2 le 0d0)

;if not keyword_set(scale) then scale = 1d0
;dataplot2 *= scale
;dataplot2 = alog(abs(dataplot2))

;loadct,5, bottom=0

;newazi = azi/pi - 1.0
;temp = min(abs(newazi - yrange(0)), grid1)
;temp = min(abs(newazi - yrange(1)), grid2)
;temp = min(abs(rad - xrange(0)), grid3)
;temp = min(abs(rad - xrange(1)), grid4)

;if not keyword_set(plotrange0) then begin 
;    plotrange=[min(dataplot2( grid3:grid4, grid1:grid2)),max(dataplot2( grid3:grid4, grid1:grid2))]
;endif else plotrange=plotrange0
;levels=plotrange(0)+(plotrange(1)-plotrange(0))*(dindgen(32)/31.)
;time=string(info(7,k)/period,format='(F7.2)')

;set_plot, 'ps'
;device, filename=filepath(strcompress('torque2d_'+ks+'.ps',/remove_all),root_dir='.',subdir=[location])$
;  ,/color, bits_per_pixel=8,xsize=12, ysize=14
;contour,dataplot2,rad,newazi ,/fill,levels=levels,title=time+' orbits', $
;  xtitle='r', ytitle=textoidl('(\phi-\phi_p)/\pi'),charsize=1.5,xmargin=[7.0,6.0], ymargin=[3,2], $
;  xtickinterval=xtickinterval,xrange=xrange, yrange=yrange, ytickinterval=.1
;colorbar, position=[0.85, 0.15, 0.90, 0.90],/vertical,/right,range=plotrange,format='(f5.2)'
;oplot,[plrad,plrad],[0.0,0.0]/pi,psym=7,symsize=1.5,color=!D.Table_size*0.0
;device,/close

print, 'done '+ks

openw,1,filepath(strcompress('torque_1d'+ks+'.dat',/remove_all),root_dir='.',subdir=[location])
for i = 0, nrad-1 do printf,1,rad(i),torque_1d(i)
close,1

openw,1,filepath(strcompress('torque_time.dat',/remove_all),root_dir='.',subdir=[location])
for i = 0, finish-start do printf,1,torque_time(0,i),torque_time(1,i),torque_time(2,i), torque_time(3,i)
close,1

endfor

;set_plot, 'ps'
;device, filename=filepath(strcompress('torque_time.ps',/remove_all),root_dir='.',subdir=[location])$
;,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches
;plot,torque_time(0,*),torque_time(1,*),xmargin=[8,4],ymargin=[3,3] $
;,ytitle=textoidl('\Gamma_3/\Gamma_T'),xtitle='t/orbits' $
;,charsize=1.5,xminor=4, thick=4
;device,/close

end
