pro nonaxi,  r0=r0, loc=loc, start=start, finish=finish, xrange=xrange, yrange=yrange, rrange=rrange, azimode=azimode, smallh=smallh  

!p.font=0
 
nmodes = 8

if not keyword_set(r0) then r0=1d0

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
 
  nlines = file_lines(filepath('planet0.dat',root_dir='.',subdir=[loc]))
  info=dblarr(11,nlines)
  openr,3,filepath('planet0.dat',root_dir='.',subdir=[location])
  readf,3,info
  close,3

  a0=info(1,0)
  dt=info(7,1)
  p0=2.*!dpi*(a0)^(3./2.)
;;;;;;;;;;;;;;;;;;
;ARRAY FOR RADIUS:
;;;;;;;;;;;;;;;;;;
  radtmp=dblarr(nrad+1)
  radius=dblarr(nrad)

  openr,1,filepath('used_rad.dat',root_dir='.',subdir=[location])
  readf,1,radtmp
  close,1
  radius(0:nrad-1)= 0.5*(radtmp(0:nrad-1)+radtmp(1:nrad))

  temp = min(abs(radius - rrange(0)),x1)
  temp = min(abs(radius - rrange(1)),x2)
  dr   = radius(x2) - radius(x1)
  
  temp = min(abs(radius - r0),r1)
  
;;;;;;;;;;
;GET DATA;
;;;;;;;;;;
  data   = dblarr(nsec,nrad)
  eigen  = dcomplexarr(nsec,nrad)
  omega0 = dblarr(nrad)

  vphi1d  = dblarr(nrad)
  kappa1d = dblarr(nrad)
  omega1d = dblarr(nrad)
  sigma1d = dblarr(nrad)

  omega1d_time = dblarr(finish-start+1,nrad)

  time      = dblarr(finish-start+1)
  time_norm = dblarr(finish-start+1) 

  mode_amp  = dblarr(nmodes,finish-start+1)
  mode_real = dblarr(nmodes,finish-start+1, nrad)  
  mode_imag = dblarr(nmodes,finish-start+1, nrad)

  mode_corot    = dblarr(nmodes,finish-start+1)
  mode_corot_ex = dblarr(nmodes,finish-start+1)

  mode_growth   = dblarr(nmodes,finish-start+1) 
  mode_growth_ex= dblarr(nmodes,finish-start+1)

  mode_cmplx_amp= dcomplexarr(nmodes,finish-start+1)

;use axisymmetric component of t=0 as normalization
  openr,2,filepath(strcompress('gasdens0.dat',/remove_all),root_dir='.',subdir=[location]) 
  readu,2,data
  close,2
  data_fft = fft(data, -1, dimension=1,/double)
  real = real_part(data_fft(0,*))
  imag = imaginary(data_fft(0,*))
  norm = dcomplex(int_tabulated(radius(x1:x2),real(x1:x2)), int_tabulated(radius(x1:x2),imag(x1:x2)))
  czero = abs(norm)

;initial angular speed
  openr,2,filepath(strcompress('gasvtheta0.dat',/remove_all),root_dir='.',subdir=[location]) 
  readu,2,data
  close,2
  for i=0, nrad-1 do vphi1d(i)  = mean(data(*,i))
  omega0  = vphi1d(r1)/radius(r1); reference frequency

;sound speed
  cs = smallh/sqrt(radius)

  for k=start, finish do begin
;time array
     time(k-start)      = info(7,k)
     time_norm(k-start) = time(k-start)/p0

;read data
     if (k lt 1000) then begin
        ks=string(k,format='(I03)')  
     endif else ks=string(k,format='(I04)')

     openr,2,filepath(strcompress('gasdens'+string(k)+'.dat',/remove_all),root_dir='.',subdir=[location]) 
     readu,2,data
     close,2
     sigma=data
     for i=0, nrad-1 do sigma1d(i)  = mean(data(*,i))     

     openr,2,filepath(strcompress('gasvrad'+string(k)+'.dat',/remove_all),root_dir='.',subdir=[location]) 
     readu,2,data
     close,2
     vrad=data
     for j = 0, nrad-2 do vrad(*,j) = (data(*,j) + data(*,j+1))/2.0
     vrad(*,nrad-1) = data(*,nrad-1)
     
     openr,2,filepath(strcompress('gasvtheta'+string(k)+'.dat',/remove_all),root_dir='.',subdir=[location]) 
     readu,2,data
     close,2
     vphi=data
     for i = 0, nsec-2 do vphi(i,*) = (data(i,*) + data(i+1,*))/2.0
     vphi(nsec-1,*) = (data(nsec-1,*) + data(0,*))/2.0

;use current angular velocity profile for getting co-rotation , lindblad resonances, and wkb propagation zones     
     for i=0, nrad-1 do vphi1d(i)  = mean(data(*,i))
     omega1d = vphi1d/radius
     omega1d_time(k-start,*) = omega1d     

     kappa1d = deriv(radius, radius*radius*vphi1d^2)/radius^3
     kappa1d = sqrt(kappa1d)

     toomreq = kappa1d*cs/(!dpi*sigma1d) 
       
;fft      
     sigma_fft    = fft(sigma, -1, dimension=1,/double)
     
     rsigmavr = sigma*vrad
     for i=0, nsec-1 do begin
        temp = rsigmavr(i,*) 
        rsigmavr(i,*) = -deriv(radius, radius*temp)/radius
     endfor
     rsigmavr_fft =  fft(rsigmavr, -1, dimension=1,/double)

     sigmavt_fft = fft(sigma*vphi, -1, dimension=1,/double)
     for i=0,nsec-1 do sigmavt_fft(i,*) *= -dcomplex(0,1d0)*i/radius(*)

     sigmam_dot = rsigmavr_fft + sigmavt_fft

     eigen = dcomplex(0,1d0)* sigmam_dot/sigma_fft
     real_freq = real_part(eigen)
     growth    = imaginary(eigen)
     
     corot = real_freq
     for m=1, nmodes do begin
        corot(m,*) = real_freq(m,*)/m
        corot(m,*) = corot(m,*)^(-2./3.)
     endfor
     
     set_plot, 'ps'
     device, filename=filepath(strcompress('nonaxi_realfreq'+ks+'.ps',/remove_all),root_dir='.',subdir=[location]) $ 
             ,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches
     plot, radius/r0, corot(1,*)/r0,xmargin=[8,2],ymargin=[3.2,1.8], ytitle=textoidl('r_c/r_0'), ystyle=0,xrange=xrange,xstyle=1, yrange=[0,10] $;yrange=[radius(0),radius(nrad-1)]/r0  $
           ,charsize=1.5, thick=4, xtitle=textoidl('r/r_0'), title = strcompress(string(info(7,k)/p0,format='(F7.2)')+textoidl('P_0'),/remove_all)
;     for m=2, nmodes do begin
;        oplot, radius/r0, corot(m,*)/r0, thick=4, linestyle=m-1
;     endfor
     device,/close

     set_plot, 'ps'
     device, filename=filepath(strcompress('nonaxi_growth'+ks+'.ps',/remove_all),root_dir='.',subdir=[location]) $ 
             ,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches
     plot, radius/r0, 1d2*growth(1,*)/omega0,xmargin=[8,2],ymargin=[3.2,1.8], ytitle=textoidl('\gamma/\Omega_0'), ystyle=0,xrange=xrange,xstyle=1, yrange=[0,1]  $
           ,charsize=1.5, thick=4, xtitle=textoidl('r/r_0'), title = strcompress(string(info(7,k)/p0,format='(F7.2)')+textoidl('P_0'),/remove_all)
;     for m=2, nmodes do begin
;        oplot, radius/r0, growth(m,*)/omega0, thick=4, linestyle=m-1
;     endfor
     device,/close


     for m=1, nmodes do begin
;amplitude of density fourier components
        real = real_part(sigma_fft(m,*))
        imag = imaginary(sigma_fft(m,*))
        
        mode_real(m-1,k-start,*) = real 
        mode_imag(m-1,k-start,*) = imag

        amp = dcomplex(int_tabulated(radius(x1:x2),real(x1:x2)), int_tabulated(radius(x1:x2),imag(x1:x2)))
        mode_amp(m-1, k-start) = abs(amp)/czero
        
        mode_cmplx_amp(m-1,k-start) = amp/czero
        
;average real frequency, find co-rotation radius
        omega = int_tabulated(radius(x1:x2),real_freq(m,x1:x2))
        omega/= dr

        temp = min(abs(m*omega1d - omega),grid)
        if( (grid eq nrad-1) or (grid eq 0) ) then begin
           rc = (omega/m)^(-2d0/3d0)
        endif else rc = radius(grid)
        mode_corot(m-1, k-start) = rc/r0 
        
;average growth rate 
        gmma  = int_tabulated(radius(x1:x2),growth(m,x1:x2))
        gmma /= dr
        mode_growth(m-1, k-start) = gmma/omega0 
        
     endfor
     
     print, 'time, amp, rc, growth ', time_norm(k-start), mode_corot(azimode-1, k-start), mode_growth(azimode-1, k-start) 
  endfor 
  
;explicit rate calculation for azimode  
  
  growth_ex = deriv(time, alog(abs(mode_cmplx_amp(azimode-1,*))) )

  rate = dcomplexarr(nmodes, finish-start+1, nrad)  

  for i=0, nrad-1 do begin 
  rate(azimode-1,*,i) = dcomplex( deriv(time, mode_real(azimode-1,*,i  )), $
                                  deriv(time, mode_imag(azimode-1,*,i) ) )

  rate(azimode-1,*,i) /= dcomplex( mode_real(azimode-1,*,i) ,             $
                                   mode_imag(azimode-1,*,i) ) 
  endfor 

  rate    *= dcomplex(0, 1d0)

  rate_re = real_part(rate)

  for i=0, finish-start do begin
     result = int_tabulated(radius(x1:x2),rate_re(azimode-1,i,x1:x2))
     omega  = result/dr

     temp = min(abs(azimode*omega1d_time(i,*) - omega),grid)
     if( (grid eq nrad-1) or (grid eq 0) ) then begin
        rc = (omega/azimode)^(-2d0/3d0)
     endif else rc = radius(grid)
     
     print, time_norm(i), rc/r0, growth_ex(i)/omega0, omega 
  endfor
  
  set_plot, 'ps'
  device, filename=filepath('fftdensity_.ps',root_dir='.',subdir=[location]) $
          , xsize=8, ysize=4.5, xoffset=0, yoffset=0 $
          , /inches,/color,bits_per_pixel=8
  plot, time_norm(*), alog10(abs(mode_cmplx_amp(azimode-1,*))), xmargin=[6,2], ymargin=[3.5,0.5] , ytitle=textoidl('log_{10}(C_m/C_0)') $
        , charsize=1.5,thick=4,linestyle=0, xtitle=textoidl('t/P_0'), xrange=xrange, yrange=yrange, xstyle=1
  device,/close
end 
