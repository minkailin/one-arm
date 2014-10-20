pro toomreQ, loc=loc, start=start, finish=finish,plotrange0=plotrange0,ct=ct $
             ,xrange=xrange, xtickinterval=xtickinterval, width=width, onedim=onedim $
             ,yrange=yrange, smallh=smallh, kep=kep
common consts, pi
pi = !dpi
h = smallh
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
nlines = file_lines(filepath('planet0.dat',root_dir='.',subdir=[location]))
info=dblarr(11,nlines)
openr,3,filepath('planet0.dat',root_dir='.',subdir=[location])
readf,3,info
close,3
a0=info(1,0)
dt=info(7,1)
p0=2.*pi*(a0)^(3./2.)
;Create 1d arrays for radial and azimuthal values for contour 
;plot.
azi=dblarr(nsec)
radtmp=dblarr(nrad+1)
openr,1,filepath('used_rad.dat',root_dir='.',subdir=[location])
readf,1,radtmp
close,1
azi=dindgen(nsec)*2.*pi/nsec

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;SETUP FOR KAPPASQ CALCULATION;
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

;Array to hold results.Don't calculate vorticity for first and last
;radius (boundaries).
kappasq=dblarr(nsec,nrad)
data = dblarr(nsec,nrad)
dataplot=dblarr(nrad,nsec)
sigmaplot=dblarr(nrad,nsec)
data_1d=dblarr(nsec)
sigma_1d=dblarr(nsec)
azi1=dblarr(nsec)
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
    
    rvtheta = radius*vtheta
    for i = 0, nsec-1 do begin
        kappasq(i, 0:nrad-1) = deriv(rmed(0:nrad-1), rvtheta(i,0:nrad-1)^2.)/(rmed(0:nrad-1)^3d0)
    endfor

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;NOW DO THE POLAR CONTOUR PLOTS;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
       if not keyword_set(kep) then begin
        data=h*sqrt(kappasq)/(pi*sigma*sqrt(radius))
       endif else begin
;        data = h*(vtheta/radius)/(pi*sigma*sqrt(radius))
         data = h*(radius^(-1.5))/(pi*sigma*sqrt(radius))
        endelse

    dataT=transpose(data)
    sigmaT=transpose(sigma)

    if not keyword_set(plotrange0) then begin 
        plotrange=[min(dataT(1:nrad-2,*)),max(dataT(1:nrad-2,*))]
    endif else plotrange=plotrange0
    levels=plotrange(0)+(plotrange(1)-plotrange(0))*(dindgen(32)/31.)
    time=string(info(7,k)/p0,format='(F7.2)')
    plx=info(1,k)
    ply=info(2,k)
    plrad=sqrt(plx*plx+ply*ply)
    phi=pltphi(plx,ply)
    temp=min(abs(azi-phi),grid)
;azi1(0:nsec-1)=azi(0:nsec-1)-azi(grid)
    
        dataplot = dataT
        sigmaplot = sigmaT
    
    if not keyword_set(onedim) then begin
        set_plot, 'ps'
        device, filename=filepath(strcompress('Qxy_'+ks+'.ps',/remove_all),root_dir='.',subdir=[location])$
          ,/color, bits_per_pixel=8,xsize=12, ysize=14
        contour,dataplot,rmed,azi/(2.0*pi),/fill,levels=levels,title=time+' orbits', $
          xtitle=textoidl('r/r_0'), ytitle=textoidl('\phi/2\pi'),charsize=1.5,xmargin=[7.0,6.0], ymargin=[3,2] $
          ,xrange=xrange, xtickinterval=xtickinterval, yrange=yrange, xstyle=1
        colorbar, position=[0.85, 0.15, 0.90, 0.90],/vertical,/right,range=plotrange,format='(f4.2)'
;        oplot,[plrad,plrad],[0.0,0.0]/pi,psym=7,symsize=1.5,color=!D.Table_size*0.5
;xyouts,plrad,phi/pi,'X',charsize=1.5
        device,/close
    endif else begin
        ;rhill = f0*plrad
        ;temp=min(abs(rmed-(plrad+width(0)*rhill)), r0)
        ;temp=min(abs(rmed-(plrad+width(1)*rhill)), r1)

        temp=min(abs(rmed-width(0)), r0)
        temp=min(abs(rmed-width(1)), r1)

        for j=0, nsec-1 do begin
            data_1d(j) = mean(dataplot(r0:r1,j))
            sigma_1d(j)= mean(sigmaplot(r0:r1,j))
        endfor
        sigma_1d *=1d3
stop
        set_plot, 'ps'
        device, filename=filepath(strcompress('Q1d_'+ks+'.ps',/remove_all),root_dir='.',subdir=[location])$
          ,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches
        plot,azi/(2.*pi), sqrt(data_1d),title=time+' orbits',xmargin=[6,6],ymargin=[4,2] $
          ,xtitle=textoidl('\phi/2\pi') $
          ,charsize=1.5,/nodata,ystyle=4
        axis, /save, yaxis=0, ytitle='Q, solid',charsize=1.5, yrange=[min(data_1d),max(data_1d)],ytickinterval=1
        oplot, azi/(2.*pi), data_1d, thick=4, linestyle=0
        axis, /save, yaxis=1, ytitle= textoidl('\Sigma\times10^3, dotted'),charsize=1.5, yrange=[min(sigma_1d),max(sigma_1d)]
        oplot, azi/(2.*pi), sigma_1d, thick=4, linestyle=1
       device,/close
   endelse
   print, 'done '+ks
endfor
end
