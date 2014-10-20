pro alpha, loc=loc, start=start, finish=finish, out=out, xrange=xrange, mp=mp $
           ,width=width, yrange=yrange

pi=3.141592654
if not keyword_set(mp) then mp=3d-4
f0 = (mp/3.0)^(1.0/3.0)
h = 0.05

if not keyword_set(finish) then finish = start

;;;;;;;;;;;;;;;
;WHERE IS DATA;
;;;;;;;;;;;;;;;
location=strcompress(loc,/remove_all)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;GET DIMENSIONS AND TIME UNIT;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
dims=(read_ascii(filepath('dims.dat',root_dir='.',subdir=[location]))).(0)
nrad=fix(dims(6))
nsec=fix(dims(7))
nout=fix(dims(5))
if not keyword_set(out) then begin
    info=dblarr(11,nout+1)
endif else info = dblarr(11,out+1)

openr,3,filepath('planet0.dat',root_dir='.',subdir=[location])
readf,3,info
close,3
a0=info(1,0)
dt=info(7,1)
p0=2.*pi*(a0)^(3./2.)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;CREATE ARRAY FOR RADIUS AND AZIMUTHAL VALUES;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
radtmp=dblarr(nrad+1)
rad=dblarr(nrad)
openr,1,filepath('used_rad.dat',root_dir='.',subdir=[location])
readf,1,radtmp
close,1
rmed = (radtmp(0:nrad-1) + radtmp(1:nrad))/2.0

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;ARRAYS TO HOLD DATA AND DATA FOR PLOTTING;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
vrad   = dblarr(nsec,nrad)
vtheta = dblarr(nsec,nrad)
vtheta_e = dblarr(nsec,nrad)

sigma0=dblarr(nsec,nrad)
sigma = dblarr(nsec,nrad)
sigma_1d = dblarr(nrad)
sigma_1d_e = dblarr(nrad)
omega = dblarr(nsec, nrad)
omega_e=dblarr(nsec, nrad)

vrad_1d=dblarr(nrad)

; temperature = dblarr(nsec,nrad)
; temperature_e = dblarr(nsec,nrad)

vtheta_1d = dblarr(nrad)
alpha = dblarr(nsec,nrad)
alpha_1d  = dblarr(nrad)
alpha_time=dblarr(3,finish-start+1)

;;;;;;;;;;;;;;;;;;;;;;;;;;;
;SETUP DONE BEGIN PLOTTING;
;;;;;;;;;;;;;;;;;;;;;;;;;;;

openr,2,filepath(strcompress('gasdens0.dat',/remove_all),root_dir='.',subdir=[location]) 
readu,2,sigma0
close,2

for k=start, finish do begin
    ks=string(k,format='(I03)')

    openr,2,filepath(strcompress('gasvrad'+string(k)+'.dat',/remove_all),root_dir='.',subdir=[location]) 
    readu,2,vrad
    close,2
    vrad_old = vrad

    openr,2,filepath(strcompress('gasvtheta'+string(k)+'.dat',/remove_all),root_dir='.',subdir=[location]) 
    readu,2,vtheta
    close,2
    
    openr,2,filepath(strcompress('gasdens'+string(k)+'.dat',/remove_all),root_dir='.',subdir=[location]) 
    readu,2,sigma
    close,2

; openr,2,filepath(strcompress('gasTemperature'+string(k)+'.dat',/remove_all),root_dir='.',subdir=[location])
; readu,2,temperature
; close,2
    
    
    vtheta_e(0:nsec-2,1:nrad-1) = (vtheta(0:nsec-2,1:nrad-1) + vtheta(1:nsec-1,1:nrad-1) + vtheta(0:nsec-2, 0:nrad-2) $ 
                                   + vtheta(1:nsec-1,0:nrad-2))/4.0
    vtheta_e(nsec-1,1:nrad-1) =  (vtheta(nsec-1,1:nrad-1) + vtheta(0,1:nrad-1) + vtheta(nsec-1, 0:nrad-2) + vtheta(0,0:nrad-2))/4.0
    
; temperature_e(0:nsec-2,1:nrad-1) = (temperature(0:nsec-2,1:nrad-1)+ temperature(0:nsec-2,0:nrad-2))/2.0
; temperature_e(nsec-1,1:nrad-1) = (temperature(nsec-1,1:nrad-1)+ temperature(0,0:nrad-2))/2.0
    
    for i = 1, nrad-1 do begin
        mean_vtheta = mean(vtheta_e(0:nsec-1,i))
        mean_vrad = mean(vrad(0:nsec-1,i))
        vtheta_e(0:nsec-1,i) -= mean_vtheta
        vrad(0:nsec-1,i)     -= mean_vrad
    endfor    
    alpha(0:nsec-1,1:nrad-1) = vtheta_e(0:nsec-1,1:nrad-1)*vrad(0:nsec-1,1:nrad-1)    
    for i = 1, nrad-1 do begin
        alpha(*,i) /= h*h/radtmp(i)
        alpha_1d(i) = mean(alpha(*,i))
    endfor

    dsigma = (sigma - sigma0)/sigma0

    alpha_time(0,k-start) = info(7,k)/p0


    if not keyword_set(width) then begin
        rpmxs = 0
        rppxs = nrad-1
    endif else begin
        temp = min(abs(rmed - width(0)),rpmxs)
        temp = min(abs(rmed - width(1)),rppxs)
    endelse
    alpha_time(2,k-start) = mean(alpha(*,rpmxs:rppxs)) ; over region near planet
  ;  subgrid = dsigma(*, rpmxs:rppxs)
  ;  overdense = where(subgrid gt 0d0)
    alpha_time(1,k-start) = mean(alpha)


;      set_plot, 'ps'
;      time=string(alpha_time(0,k-start),format='(F7.2)')
;      device, filename=filepath(strcompress('alpha_'+ks+'.ps',/remove_all),root_dir='.',subdir=[location])$
;        ,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches
;      plot,rmed, alpha_1d*1d4,title=time+' orbits',xmargin=[6,2],ymargin=[3,2] $
;        ,ytitle=textoidl('<\alpha>_\phi\times10^{4}'),xtitle='r' $
;        ,charsize=1.5,xrange=xrange,xminor=4, thick=4, yrange=yrange
;      device,/close
    
    openw,1,filepath(strcompress('alpha_'+ks+'.dat',/remove_all),root_dir='.',subdir=[location])
    for i = 1, nrad-1 do printf,1, radtmp(i), alpha_1d(i)
    close,1

    ;openw,1,filepath(strcompress('shear_'+ks+'.dat',/remove_all),root_dir='.',subdir=[location])
    ;for i = 0, nrad-1 do printf,1, rmed(i), shear(i)
    ;close,1
    
    print, 'done '+ks
endfor

; set_plot, 'ps'
; device, filename=filepath(strcompress('alpha_time.ps',/remove_all),root_dir='.',subdir=[location])$
;   ,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches
; plot,alpha_time(0,*),alpha_time(1,*),xmargin=[8,4],ymargin=[3,3] $
;   ,ytitle=textoidl('<\alpha>\times10^{3}'),xtitle='t' $
;   ,charsize=1.5,xminor=4, thick=4
; device,/close 
  

openw,1,filepath(strcompress('alpha_time.dat',/remove_all),root_dir='.',subdir=[location])
for i = 0, finish-start do printf,1,alpha_time(0,i),alpha_time(1,i), alpha_time(2,i)
close,1

end


;  alpha from 1d alpha-disc model
;     if keyword_set(visc1d) then begin
;         for i = 0, nrad-1 do begin
;             omega(0:nsec-1,i) = vtheta(0:nsec-1,i)/rmed(i) 
;             sigma_1d(i) = mean(sigma(*,i))
;             vrad_1d(i) = mean(vrad_old(*,i))
;         endfor
;         for i = 1, nrad-1 do begin
;             omega_e(0:nsec-1,i) = (omega(0:nsec-1,i) + omega(0:nsec-1,i-1))/2.0   
;             sigma_1d_e(i) = (sigma_1d(i) + sigma_1d(i-1))/2.0 
;         endfor

;         hprime = deriv(radtmp(1:nrad-1), omega_e(1:nrad-1)*radtmp(1:nrad-1)*radtmp(1:nrad-1))
;         g = radtmp(1:nrad-1)*sigma_1d_e(1:nrad-1)*vrad_1d(1:nrad-1)*hprime

;         omegaprime = deriv(rmed, omega)
;         alpha_visc1d(rpmxs) = (alpha_1d(rpmxs)/1d3)*h*h*sigma_1d(rpmxs)*rmed(rpmxs)^(3.5)*omegaprime(rpmxs)
;         for i= rpmxs+1, nrad-1 do begin
;             alpha_visc1d(i) = g(i-1)*(rmed(i) - rmed(i-1)) + alpha_visc1d(i-1)
;         endfor

;         f = h*h*rmed^(3.5)*sigma_1d*omegaprime
;         alpha_visc1d /= f
        
;         print, alpha_visc1d

;         for i=0, nsec-1 do begin
;             alpha_visc2d(i, 1:nrad-1) = omega_e(i,1:nrad-1)*vrad_old(i, 1:nrad-1)/deriv(radtmp(1:nrad-1), omega_e(i, 1:nrad-1))
;         endfor
;         for i=0, nrad-1 do alpha_visc1d(i) = mean(alpha_visc2d(*,i))
        
;         alpha_visc1d(1:nrad-1)/=h*h*sqrt(radtmp(1:nrad-1))

;         print, alpha_visc1d

;     endif
