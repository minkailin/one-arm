

pro disk_modes, loc=loc, start=start, finish=finish, legend=legend, xrange=xrange, yrange=yrange $
                ,plotmodes=plotmodes, xtickinterval=xtickinterval, job=job $
                ,mp=mp, basic=basic
  common consts, pi
  pi=!dpi

  f0=(mp/3.0)^(1.0/3.0)
  if not keyword_set(finish) then finish=start
  
  modes = 10
  
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;DATA LOCATION,GET RUN INFO, READ PLANET ORBIT INFO.;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  location=strcompress(loc,/remove_all)
  dims=(read_ascii(filepath('dims.dat',root_dir='.',subdir=[location]))).(0)
  nout=fix(dims(5))
  nrad=fix(dims(6))
  r1 = nrad - 1
  nsec=fix(dims(7))
  if not keyword_set(out) then begin
     info=dblarr(11,nout+1)
  endif else info=dblarr(11,out+1)
  nlines = file_lines(filepath('planet0.dat',root_dir='.',subdir=[location]))
  info=dblarr(11,nlines)
  
  openr,3,filepath('planet0.dat',root_dir='.',subdir=[location])
  readf,3,info
  close,3
  a0=info(1,0)
  dt=info(7,1)
  p0=2.*pi*(a0)^(3./2.)
  
  azi=dblarr(nsec)
  radtmp=dblarr(nrad+1)
  openr,1,filepath('used_rad.dat',root_dir='.',subdir=[location])
  readf,1,radtmp
  close,1
  rmed=radtmp(0:nrad-1)+radtmp(1:nrad)
  rmed=rmed*0.5
  
  azi=dindgen(nsec)*2.*pi/nsec
  radius=dblarr(nsec,nrad)
  for l=0, nsec-1 do radius(l,0:nrad-1)=rmed(0:nrad-1)
  
;;;;;;;
;SETUP;
;;;;;;;
  
;Arrays to hold raw data, and center-valued velocities.
  sigma  = dblarr(nsec,nrad)
  vtheta = dblarr(nsec,nrad)
  vrad   = dblarr(nsec,nrad)
  
  sigma0  = dblarr(nsec,nrad)
  vtheta0 = dblarr(nsec,nrad)
  vrad0   = dblarr(nsec,nrad)
  
  sigma_1d   = dblarr(nrad) 
  sigma_1d0  = dblarr(nrad)
  vtheta_1d  = dblarr(nrad)
  omega0_1d  = dblarr(nrad) 
  
  vthetac = dblarr(nsec,nrad)
  vradc   = dblarr(nsec,nrad)
        
  fft_sigma      = dcomplexarr(nsec,nrad)
  fft_sigmavrad  = dcomplexarr(nsec,nrad)
  fft_sigmavtheta= dcomplexarr(nsec,nrad)
  
  am    = dcomplexarr(modes+1, nrad)
  amdot = dcomplexarr(modes+1, nrad)
  
  integrated_amp = dcomplexarr(modes+1)
  integrated_freq = dcomplexarr(modes+1)
  growthrate     = dblarr(modes)
  corotation     = dblarr(modes)
  
;Basic state angular velocity
  openr,3,filepath(strcompress('gasvtheta'+string(basic)+'.dat',/remove_all),root_dir='.',subdir=[location]) 
  readu,3,vtheta0
  close,3
  for i=0, nrad-1 do omega0_1d(i) = mean(vtheta0(*,i))/rmed(i)
  
  case job of
     'calc': begin

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;GET DATA AND INTERPOLATE TO GET CELL-CENTER VELOCITIES;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
        openr,2,filepath(strcompress('gasdens0.dat',/remove_all),root_dir='.',subdir=[location]) 
        readu,2,sigma
        close,2
        for i = 0, nrad -1 do begin 
           sigma_1d0(i) = mean(sigma(*,i))
        endfor
        
        openw,10,filepath(strcompress('mode_amplitudes.dat',/remove_all),root_dir='.',subdir=[location])
        openw,20,filepath(strcompress('mode_growthrate.dat',/remove_all),root_dir='.',subdir=[location])
        openw,30,filepath(strcompress('mode_corotation.dat',/remove_all),root_dir='.',subdir=[location])
        
        for k=start, finish do begin
           time = info(7,k)/p0
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
           
           for i = 0, nsec-2 do vthetac(i,*) = (vtheta(i,*) + vtheta(i+1,*))/2.0
           vthetac(nsec-1,*) = (vtheta(nsec-1,*) + vtheta(0,*))/2.0
           
           for j = 0, nrad-2 do vradc(*,j) = (vrad(*,j) + vrad(*,j+1))/2.0
           vradc(*,nrad-1) = vrad(*,nrad-1)
           
           sigmavrad   = radius*sigma*vradc
           
           for i = 0, nsec-1 do begin
              sigmavrad(i,*) = deriv(rmed(*), sigmavrad(i,*)) ;d/dr(r sigma u_r)
           endfor
           sigmavtheta = sigma*vthetac
           
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;FFT w.r.t. phi of sigma, sigmavrad and sigmavtheta;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
           
           for i=0, nrad-1 do begin
              fft_sigma(*,i)       = fft(sigma(*,i),/double) 
              fft_sigmavrad(*,i)   = fft(sigmavrad(*,i),/double) 
              fft_sigmavtheta(*,i) = fft(sigmavtheta(*,i),/double)
           endfor
           
;Extract the m=modes components
           
           am = fft_sigma(0:modes,*)
 
;Get amdot 
           for i = 0, nrad-1 do begin
              amdot(0:modes, i) = -(1d0/rmed(i))*fft_sigmavrad(0:modes,i) $
                                  - dcomplex(0d0, 1d0)*dindgen(modes+1)*fft_sigmavtheta(0:modes,i)/rmed(i)
           endfor
           
           complex_freq = dcomplex(0d0, 1d0)*amdot/am
           
;Integrate the modes, outer disc
           
           plx = info(1, k)
           ply = info(2, k)
           plrad = sqrt(plx*plx + ply*ply)
           rh    = f0*plrad
           temp = min(abs(rmed - plrad), r0)
               
           drad = rmed(r1) - rmed(r0)
           
           for i=0, modes do begin
              integrated_amp(i) = dcomplex(int_tabulated(rmed(r0:r1),real_part(am(i,r0:r1)),/double) $
                                           ,int_tabulated(rmed(r0:r1),imaginary(am(i,r0:r1)),/double))
              integrated_freq(i) = dcomplex(int_tabulated(rmed(r0:r1),real_part(complex_freq(i,r0:r1)),/double) $
                                           ,int_tabulated(rmed(r0:r1),imaginary(complex_freq(i,r0:r1)),/double))
           endfor
           integrated_amp  /= drad
           integrated_freq /= drad
           
;Mode amplitudes
           printf,10, time, abs(integrated_amp(1:modes)/integrated_amp(0)) $ 
                  ,format=strcompress('('+string(modes+1)+'(e10.4,2x))',/remove_all) 

;Mode co-rotation radius and growth rates
           for i=1, modes do begin
              
              real_freq = real_part(integrated_freq(i)))/double(i)
              temp = min(abs(omega0_1d - real_freq, grid)
              rc = rmed(grid)
              
              if (omega0_1d(grid) lt  real_freq) then begin
                 rc = mean(rmed(grid:grid-1))
                 omega_c = mean(oemga0_1d(grid:grid-1))
              endif
              
              if (omega0_1d(grid) gt  real_freq) then begin
                 rc = mean(rmed(grid:grid+1))
                 omega_c = mean(omega0_1d(grid:grid+1))
              endif

              corotation(i-1) = (rc - plrad)/rh
              growthrate(i-1) = imaginary(integrated_freq(i))/omega_c
              
           endfor
           
           printf,20, time, growthrate(0:modes-1) $ 
                  ,format=strcompress('('+string(modes+1)+'(e10.4,2x))',/remove_all) 
           printf,30, time, corotation(0:modes-1) $ 
                  ,format=strcompress('('+string(modes+1)+'(e10.4,2x))',/remove_all) 
           
           close, 10
           close, 20
           close, 30
        end
     end
end

