function func, phi, parameters
m = parameters(0)
r = parameters(1)
rprime = parameters(2)
h = parameters(3)

soft = 0.3*rprime*h
result = cos(m*phi)/sqrt(r*r + rprime*rprime - 2.0*r*rprime*cos(phi)+ soft*soft)

return, result
end

function shift, data, plx, ply
common shared, nrad, nsec, azi

dataplot = dblarr(nrad, nsec)

dataT = transpose(data)
phi=pltphi(plx,ply)
temp=min(abs(azi-phi),grid)
    
if(phi gt !dpi) then begin
    temp2=min(abs(azi-(phi-!dpi)),grid2)
    dataplot(0:nrad-1,0:nsec-1-grid2) = dataT(0:nrad-1,grid2:nsec-1)
    dataplot(0:nrad-1,nsec-grid2:nsec-1) = dataT(0:nrad-1, 0:grid2-1) 
endif
if(phi lt !dpi) then begin
    temp2=min(abs(azi-(phi+!dpi)),grid2)
    dataplot(0:nrad-1,nsec-grid2:nsec-1) = dataT(0:nrad-1,0:grid2-1)
    dataplot(0:nrad-1,0:nsec-1-grid2) = dataT(0:nrad-1, grid2:nsec-1) 
endif
if(phi eq !dpi) then dataplot = dataT

return, dataplot
end

pro stress3, loc=loc, start=start, finish=finish, out=out, xrange=xrange, modes=modes
common shared, nrad, nsec, azi

scale = 4
pi=3.141592654
G = 1.0
mp = 3d-4
f0 = (mp/3.0)^(1.0/3.0)
h = 0.05
m = modes

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
radtmp = dblarr(nrad+1)
rad = dblarr(nrad)
azi = dblarr(nsec)

openr,1,filepath('used_rad.dat',root_dir='.',subdir=[location])
readf,1,radtmp
close,1

rad(0:nrad-1) = (radtmp(0:nrad-1) + radtmp(1:nrad))/2.0
dlogr = alog(radtmp(1)/radtmp(0))/double(nrad/scale)
for i=0, nsec-1 do azi(i)=2.*pi*i/nsec

csq = h*h/rad

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; ARRAYS              ;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
kernel_small = dblarr(nrad/scale, nrad/scale)
sigma  =  dblarr(nsec,nrad)
vrad  =  dblarr(nsec,nrad)
vtheta  =  dblarr(nsec,nrad)

vrad_fft = dcomplexarr(nsec,nrad)
vrad_1d = dcomplexarr(nrad)

vtheta_fft = dcomplexarr(nsec,nrad)
vtheta_1d = dcomplexarr(nrad)
omega0 = dblarr(nrad)
kappasq0 = dblarr(nrad)
sigma0_1d = dblarr(nrad)

sigma_fft = dcomplexarr(nsec,nrad)
sigma_1d = dcomplexarr(nrad)

potential = dcomplexarr(nrad)
dpotential = dcomplexarr(nrad)
flux = dblarr(nrad)
flux_reynolds = dblarr(nrad)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; FILL KERNEL        ;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

 rad_small = rebin(rad,nrad/scale)

 for i=0, nrad/scale-1 do begin
     print, i
     for j=0, nrad/scale-1 do begin
         rprime = rad_small(j)
         r = rad_small(i)
         parameters = [m, r, rprime, h]
         kernel_small(i,j) = qpint1d('func', 0d0, 2d0*pi, parameters, epsabs=1d-16, epsrel=1d-16)
          kernel_small(i,j)*= -G*dlogr*rprime*rprime
      endfor
  endfor

  kernel = rebin(kernel_small,nrad,nrad)

; for i=0, nrad-1 do begin
;     print, i
;     for j=0, nrad-1 do begin
;         rprime = rad(j)
;         r = rad(i)
;         parameters = [m, r, rprime, h]
;         kernel(i,j) = qpint1d('func', 0d0, 2d0*pi, parameters, epsabs=1d-16, epsrel=1d-16)
;         kernel(i,j)*= -G*dlogr*rprime*rprime
;     endfor
; endfor

print, 'got the kernel'

plx=info(1,20)
ply=info(2,20)

openr,2,filepath(strcompress('gasdens20.dat',/remove_all),root_dir='.',subdir=[location]) 
readu,2,sigma
close,2 
sigma0 = shift(sigma, plx, ply)

openr,2,filepath(strcompress('gasvrad20.dat',/remove_all),root_dir='.',subdir=[location]) 
readu,2,vrad
close,2
vrad0 = shift(vrad, plx, ply)

openr,2,filepath(strcompress('gasvtheta20.dat',/remove_all),root_dir='.',subdir=[location]) 
readu,2,vtheta
close,2
vtheta0 = shift(vtheta, plx, ply)

for i=0, nrad-1 do begin
    omega0(i) = mean(vtheta0(*,i))/rad(i)
    sigma0_1d(i) = mean(sigma0(*,i))
endfor
kappasq0 =  deriv(rad, rad^4.*omega0^2.)/(rad^3.0)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;GET THE SELF-GRAVTY ACCELERATIONS;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
for l=start, finish do begin
    ks=string(l,format='(I03)')

    plx=info(1,l)
    ply=info(2,l)

    openr,2,filepath(strcompress('gasdens'+string(l)+'.dat',/remove_all),root_dir='.',subdir=[location]) 
    readu,2,sigma
    close,2 
    sigma = shift(sigma, plx, ply) - sigma0
;    for i=0, nrad-1 do sigma(*,i) -= sigma0_1d(i)
    sigma_fft = 2.0*fft(sigma, dimension=1,/double)
    sigma_1d(*) = sigma_fft(modes,*)

    openr,2,filepath(strcompress('gasvrad'+string(l)+'.dat',/remove_all),root_dir='.',subdir=[location]) 
    readu,2,vrad
    close,2 
    vrad = shift(vrad, plx, ply) - vrad0
;    for i=0, nrad-1 do vrad(*,i) -= mean(vrad0(*,i))
    vrad_fft = 2.0*fft(vrad, dimension=1,/double)
    vrad_1d(*) = vrad_fft(modes,*)

    openr,2,filepath(strcompress('gasvtheta'+string(l)+'.dat',/remove_all),root_dir='.',subdir=[location]) 
    readu,2,vtheta
    close,2 
    vtheta = shift(vtheta, plx, ply) - vtheta0
;    for i=0, nrad-1 do vtheta(*,i) -= mean(vtheta0(*,i))
    vtheta_fft = 2.0*fft(vtheta, dimension=1,/double)
    vtheta_1d(*) = vtheta_fft(modes,*)

    ;SG flux
    potential = kernel#sigma_1d
    dpotential = dcomplex(deriv(rad, real_part(potential)), deriv(rad, imaginary(potential)))
    ;print, 'sigma', mean(abs(sigma_1d)),  ' potential', mean(abs(potential)), ' dpotential', mean(abs(dpotential))

    flux_sg = 2.0*h*rad*(modes/2.0)*imaginary(dpotential*conj(potential))/(4.0*pi*G)
    
    ;flux from lagragian displacement

     dsigma = sigma_1d/sigma0_1d
     diff_sigma = dcomplex(deriv(rad, real_part(dsigma)), deriv(rad, imaginary(dsigma)))
    
     w_nsg = csq*dsigma
     wprime_nsg = csq*diff_sigma

     w = w_nsg + potential
     wprime = wprime_nsg + dpotential
    
     zeta = vrad_1d/vtheta_1d

     top = zeta*wprime*kappasq0/(2.0*omega0) + 2.0*dcomplex(0d0, 1d0)*modes*omega0*w/rad
     bot = dcomplex(0d0, 1d0)*wprime + zeta*modes*w/rad
     sbar = -top/bot

     displacement = -dcomplex(0d0, 1d0)*vrad_1d/sbar

     flux_disp = -(modes/2d0)*sigma*imaginary(conj(displacement)*w)

     ; reynolds flux

     flux_rey = 0.5*rad*sigma*real_part(vtheta_1d*conj(vrad_1d))


    openw,1,filepath(strcompress('stress3_'+ks+'.dat',/remove_all),root_dir='.',subdir=[location])
    for i = 1, nrad-1 do printf,1, rad(i), flux_sg(i), flux_disp(i), flux_rey(i)
    close,1   
 
    print, 'done '+ks

    set_plot, 'ps'
    time=string(info(7,l)/p0,format='(F7.2)')
    device, filename=filepath(strcompress('stress3_sg_'+ks+'.ps',/remove_all),root_dir='.',subdir=[location])$
      ,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches
    plot,rad, flux_sg/abs(flux_sg(nrad/2)), title=time+' orbits',xmargin=[6,1],ymargin=[3,2] $
      ,ytitle=textoidl('Ang. mom. flux'),xtitle='r' $
      ,charsize=1.5,xrange=xrange, thick=4, xtickinterval=1
    oplot, rad, flux_disp/abs(flux_disp(nrad/2)), thick=4, linestyle=1
    oplot, rad, flux_rey/abs(flux_rey(nrad/2)), thick=4, linestyle=2
    device,/close
    
;     set_plot, 'ps'
;     time=string(info(7,l)/p0,format='(F7.2)')
;     device, filename=filepath(strcompress('stress3_rey_'+ks+'.ps',/remove_all),root_dir='.',subdir=[location])$
;       ,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches
;     plot,rad, flux_reynolds_mode,title=time+' orbits',xmargin=[6,1],ymargin=[3,2] $
;       ,ytitle=textoidl('Ang. mom. flux'),xtitle='r' $
;       ,charsize=1.5,xrange=xrange, thick=4, xtickinterval=1
;     ;oplot, rad, flux_reynolds_mode, linestyle=1, thick=4
;     device,/close
endfor

; set_plot, 'ps'

; device, filename=filepath(strcompress('alpha_time.ps',/remove_all),root_dir='.',subdir=[location])$
; ,bits_per_pixel=8,xsize=8, ysize=6,xoffset=0,yoffset=0,/inches
; plot,alpha_time(0,*),alpha_time(1,*),xmargin=[8,4],ymargin=[3,3] $
; ,ytitle=textoidl('<\alpha>\times10^{3}'),xtitle='t' $
; ,charsize=1.5,xminor=4, thick=4, xrange=[20,200]
; device,/close

; openw,1,filepath(strcompress('alphaSG_time.dat',/remove_all),root_dir='.',subdir=[location])
; for i = 0, finish-start do printf,1,alpha_time(0,i),alpha_time(1,i), alpha_time(2,i)
; close,1

end
