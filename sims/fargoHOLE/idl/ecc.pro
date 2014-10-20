pro ecal, k=k
common consts, pi, nrad, time
common shared, location, sigma, vrad, vtheta, nsec, vthetac, vradc, radius, eccentricity, azi

;get raw data
ks=string(k,format='(I03)')
openr,2,filepath(strcompress('gasdens'+string(k)+'.dat',/remove_all),root_dir='.',subdir=[location]) 
readu,2,sigma
close,2
openr,2,filepath(strcompress('gasvrad'+string(k)+'.dat',/remove_all),root_dir='.',subdir=[location]) 
readu,2,vrad
close,2
openr,2,filepath(strcompress('gasvtheta'+string(k)+'.dat',/remove_all),root_dir='.',subdir=[location]) 
readu,2,vtheta
close,2

;interpolate to get velocities at grid centres
vthetac(0,1:nrad-1) = 0.5*(vtheta(0,1:nrad-1) + vtheta(nsec-1,1:nrad-1))

vthetac(1:nsec-1,1:nrad-1) = 0.5*(vtheta(1:nsec-1,1:nrad-1) + vtheta(0:nsec-2,1:nrad-1))
vradc(1:nsec-1,1:nrad-1)   = 0.5*(vrad(1:nsec-1,1:nrad-1) + vrad(1:nsec-1,0:nrad-2))

;j=0 is not valid!
vsq = vradc*vradc + vthetac*vthetac
rvrad = radius*vradc
esq = vsq*vsq*radius*radius - vsq*rvrad*rvrad-2.0*vsq*radius+2.0*rvrad*rvrad/radius + 1.0
esq = 1.0+radius*vthetac*vthetac*(radius*vsq - 2.0)

eccentricity=sqrt(esq)

end

pro ecc, loc=loc, start=start, finish=finish,ct=ct,plotrange0=plotrange0, width=width, pol=pol, out=out,xrange=xrange, yrange=yrange
if not keyword_set(ct) then ct=5
if not keyword_set(finish) then finish=start
common consts, pi, nrad, time
common shared, location, sigma, vrad, vtheta, nsec, vthetac, vradc, radius, eccentricity, azi
pi=3.141592654
mp = 3d-4
f0=(mp/3.0)^(1./3.)
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
if not keyword_set(out) then begin
    info=dblarr(11,nout+1)
endif else info = dblarr(11, out+1)

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
radius=dblarr(nsec,nrad)
openr,1,filepath('used_rad.dat',root_dir='.',subdir=[location])
readf,1,radtmp
close,1

dlogr = alog(radtmp(nrad)/radtmp(0))/nrad
dtheta = 2.0*pi/nsec

for i=0, nsec-1 do azi(i)=dtheta*i
for j=0, nrad-1 do rad(j)=(radtmp(j)+radtmp(j+1))/2.
for i=0, nsec-1 do radius(i,*)=rad(*)

;;;;;;;;;;;;;;;;;;;;;;;;
;DO POLAR CONTOUR PLOTS;
;;;;;;;;;;;;;;;;;;;;;;;;
sigma = dblarr(nsec,nrad)
vrad  = dblarr(nsec,nrad)
vtheta= dblarr(nsec,nrad)
vradc = dblarr(nsec,nrad)
vthetac=dblarr(nsec,nrad)
eccentricity=dblarr(nsec,nrad)
rphi = dblarr(nsec)
std = dblarr(3)
ecc_time = dblarr(4,finish-start +1)

;basic state

openr,2,filepath(strcompress('gasdens0.dat',/remove_all),root_dir='.',subdir=[location]) 
readu,2,sigma
close,2

sigma0=sigma

for k=start, finish do begin
    ks=string(k,format='(I03)')
    time=string(info(7,k)/p0,format='(F7.2)')

    plx=info(1,k)
    ply=info(2,k)
    plrad=sqrt(plx*plx+ply*ply)
    phi=pltphi(plx,ply)
    azi1(0:nsec-1)=azi(0:nsec-1)-phi
;rhill = plrad*f0
;xs = 2.5*rhill
    temp = min(abs(rad - width(0)),x1)
    temp = min(abs(rad - width(1)),x2)
    
    openr,2,filepath(strcompress('gasdens'+string(k)+'.dat',/remove_all),root_dir='.',subdir=[location]) 
    readu,2,sigma
    close,2

    ecal, k=k
    ;subarray = eccentricity(*,x1:x2)
    test = (sigma(*,x1:x2) - sigma0(*,x1:x2))/sigma0(*,x1:x2)
    testrad = rad(x1:x2)
   
    for i=0, nsec-1 do begin
        filter = min(abs(test(i,*)), edge)
        rphi(i) = testrad(edge)
    endfor
    rphi = smooth(rphi,10)
    xellip = rphi*cos(azi)
    yellip = rphi*sin(azi)
  
    params = mpfitellipse(xellip, yellip, start_params, /tilt,/quiet)

    a = params(0)
    b = params(1)
    e = sqrt(1.0 - (b/a)^2)
    phi_e = params(4)

    theta = azi + phi_e - pi/2.

    rfit = ((cos(theta)/b)^2 + (sin(theta)/a)^2)^(-0.5)
    x = (-a + 2*a*dindgen(101)/100.)
    y = tan(pi-phi_e-phi)*x

    ecc_time(0,k-start) = info(7,k)/p0
    ecc_time(1,k-start) = e
    ecc_time(2,k-start) =  phi_e

    rhill = plrad*f0
    xs = 2.5*rhill
    temp = min(abs(rad - (plrad-xs)),x1)
    temp = min(abs(rad - (plrad+xs)),x2)
    test = (sigma(*,x1:x2) - sigma0(*,x1:x2))/sigma0(*,x1:x2)
    ecc_time(3,k-start) = mean(test)

    if keyword_set(pol) then begin
        loadct,ct, bottom=0

        if not keyword_set(plotrange0) then begin 
            plotrange=[min(eccentricity),max(eccentricity)]
        endif else plotrange=plotrange0
        levels=plotrange(0)+(plotrange(1)-plotrange(0))*(dindgen(32)/31.)
        
        set_plot, 'ps'
        device, filename=filepath(strcompress('ecc_'+ks+'.ps',/remove_all),root_dir='.',subdir=[location])$
          ,/color, bits_per_pixel=8,xsize=16, ysize=12
        polar_contour,eccentricity(*,1:nrad-1),azi1(*),rad(1:nrad-1),/isotropic,/fill,levels=levels,title=time+' orbits' $
          ,ymargin=[2.5,2.5],xrange=yrange, yrange=yrange $
          , xtickinterval=2., ytickinterval=2.,/dither     
        colorbar, position=[0.8, 0.1, 0.85, 0.9],/vertical,/right,range=plotrange,format='(e10.2)'
;        oplot,[plrad,plrad],[0.,0.],psym=6,symsize=1,color=120
        device,/close
    endif 

     set_plot, 'ps'
     device, filename=filepath(strcompress('ellipse_'+ks+'.ps',/remove_all),root_dir='.',subdir=[location])$
       , xsize=6, ysize=6, xoffset=0, yoffset=0, /inches,/color,bits_per_pixel=8
     plot, rphi, azi, xmargin=[4,4], ymargin=[4,4], ytitle=yytitle $
       , charsize=1.5,thick=4,linestyle=0, title=time+' orbits', xrange=xrange, yrange=yrange $
       , /polar,/isotropic
     ;oplot, b*cos(azi), a*sin(azi), linestyle=2
     oplot, rfit, azi1, linestyle=2, /polar
     oplot, x, y, linestyle=1
     device,/close
    



;     array=dblarr(nrad-1)
;     for j=0, nrad-2 do array(j) = mean(eccentricity(*,j+1))
    
;     set_plot, 'ps'
;     device, filename=filepath(strcompress('ecc1d_'+ks+'.ps',/remove_all),root_dir='.',subdir=[location])$
;       , xsize=8, ysize=4.5, xoffset=0, yoffset=0, /inches,/color,bits_per_pixel=8
;     plot, rad(1:nrad-1), array, xmargin=[6,2], ymargin=[3,2], ytitle=yytitle $
;       , charsize=1.5,thick=4,linestyle=0, xtitle='r', title=time+' orbits', xrange=xrange, yrange=yrange
;     device,/close
    print, 'done '+ks
endfor

 set_plot, 'ps'
 device, filename=filepath(strcompress('ecc_time.ps',/remove_all),root_dir='.',subdir=[location])$
   , xsize=8, ysize=4.5, xoffset=0, yoffset=0, /inches,/color,bits_per_pixel=8
 plot, ecc_time(0,*), smooth(ecc_time(2,*),3)/pi, xmargin=[6,2], ymargin=[3,2], ytitle=yytitle $
   , charsize=1.5,thick=4,linestyle=0, xtitle='t/orbits'
 ;oplot, ecc_time(0,*), ecc_time(3,*)/pi, linestyle=1
 device,/close

 openw,1,filepath(strcompress('ecc_time.dat',/remove_all),root_dir='.',subdir=[location])
 for j=0, finish-start do printf,1,ecc_time(0,j),ecc_time(1,j),ecc_time(2,j), ecc_time(3,j),format='(4(e12.5,1x))'
 close,1

end
