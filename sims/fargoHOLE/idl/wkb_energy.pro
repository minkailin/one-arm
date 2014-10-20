function epicyclesq, rad, vtheta
  kappasq = deriv(rad, rad*rad*vtheta*vtheta)/(rad^3.0)
  return, kappasq
end

function soundspeed, pressure, sigma
  return, sqrt(pressure/sigma)
end

function toomreQ, rad, vtheta, sigma, pressure
  cs    = soundspeed(pressure, sigma)
  kappa = sqrt(epicyclesq(rad,vtheta))
  return, cs*kappa/(!dpi*sigma)
end

function get_grid, loc, start, r0
  
  dims=(read_ascii(filepath('dims.dat',root_dir='.',subdir=[loc]))).(0)
  nout=fix(dims(5))
  nrad=fix(dims(6))
  nsec=fix(dims(7))
  
  radtmp=dblarr(nrad+1)
  rad   =dblarr(nrad)
  rplot =dblarr(nrad)
  
  nlines = file_lines(filepath('planet0.dat',root_dir='.',subdir=[loc]))
  info=dblarr(11,nlines)
  openr,3,filepath('planet0.dat',root_dir='.',subdir=[loc])
  readf,3,info
  close,3
 
  openr,1,filepath('used_rad.dat',root_dir='.',subdir=[loc])
  readf,1,radtmp
  close,1
  rad(0:nrad-1)=(radtmp(0:nrad-1)+radtmp(1:nrad))/2.
  
  rplot = rad/r0
   
  return, rplot 
end

pro wkb_energy, loc=loc, start=start, finish=finish, legend=legend,label=label, xrange=xrange, yrange=yrange, ytickinterval=ytickinterval, r0=r0, xtickinterval=xtickinterval, rrange=rrange, smallh=smallh, azimodes=azimodes, rc=rc, scale=scale

  !p.font = 0  
  nmodes = n_elements(azimodes)
  location=strcompress(loc,/remove_all)
  if not keyword_set(scale) then scale=1

  if not keyword_set(r0) then r0 = 1.0 
  rplot = get_grid(loc, start, r0)
  rad   = rplot*r0
  xtitle=textoidl('r/r_0') 
  temp = min(abs(rplot - rrange(0)), r1)
  temp = min(abs(rplot - rrange(1)), r2)
  
  dims=(read_ascii(filepath('dims.dat',root_dir='.',subdir=[location]))).(0)
  nrad=fix(dims(6))
  nsec=fix(dims(7))

  azi         = 2d0*!dpi*dindgen(nsec)/(nsec - 1d0)
  
  data      = dblarr(nsec,nrad)
  sigma1d   = dblarr(nrad) 
  vtheta1d  = dblarr(nrad)

  time        = dblarr(finish-start + 1)
  mode_energy = dblarr(nmodes, finish-start + 1)
  
  nlines = file_lines(filepath('planet0.dat',root_dir='.',subdir=[loc]))
  info=dblarr(11,nlines)
  openr,3,filepath('planet0.dat',root_dir='.',subdir=[location])
  readf,3,info
  close,3  
  p0=2d0*!dpi*(r0)^(3./2.)
  
;  pattern = rc^(-1.5)

  csq = smallh^2/rad 

  for k=start, finish do begin

     if (k lt 1000) then begin
        ks=string(k,format='(I03)')
     endif else ks=string(k,format='(I04)')
     
     time(k-start)   = info(7,k)/p0
     
     openr,2,filepath(strcompress('gasdens'+string(k)+'.dat',/remove_all),root_dir='.',subdir=[location]) 
     readu,2,data
     close,2
     sigma = data 
     for i=0, nrad-1 do sigma1d(i) = mean(data(*,i))
     fft_sigma = fft(data, -1, dimension=1, /double)

     openr,2,filepath(strcompress('gasvtheta'+string(k)+'.dat',/remove_all),root_dir='.',subdir=[location]) 
     readu,2,data
     close,2
     vtheta = data
     for i = 0, nsec-2 do vtheta(i,*) = (data(i,*) + data(i+1,*))/2.0
     vtheta(nsec-1,*) = (data(nsec-1,*) + data(0,*))/2.0
     for i=0, nrad-1 do vtheta1d(i) = mean(data(*,i))
     omega  = vtheta1d/rad
     kappa2 = epicyclesq(rad, vtheta1d)
     
     openr,2,filepath(strcompress('gasvrad'+string(k)+'.dat',/remove_all),root_dir='.',subdir=[location]) 
     readu,2,data
     close,2
     vrad=data
     for j = 0, nrad-2 do vrad(*,j) = (data(*,j) + data(*,j+1))/2.0
     vrad(*,nrad-1) = data(*,nrad-1)
               
     sigmavr_fft =  fft(sigma*vrad, -1, dimension=1,/double)
     for i=0,nsec-1 do sigmavr_fft(i,*) /= rad(*)

     sigmavt_fft = fft(sigma*vtheta, -1, dimension=1,/double)
     for i=0,nsec-1 do sigmavt_fft(i,*) *= dcomplex(0,1d0)*i/rad(*)

     for i = 0, nmodes-1 do begin
        m       = azimodes(i)
        sigma_m = fft_sigma(m,*)
    
        re  = real_part(sigma_m)
        dre = deriv(rad, re)
        
        im  = imaginary(sigma_m)
        dim = deriv(rad, im)
       
        kr  = dim*re - dre*im
        kr /= re^2 + im^2 + 1d-16

        if not keyword_set(rc) then begin
           sigmam = dcomplex(int_tabulated(rad(r1:r2), re(r1:r2)) ,$
                             int_tabulated(rad(r1:r2), im(r1:r2)))
           
           sigmam_dot = sigmavr_fft(m,r2) - sigmavr_fft(m,r1)
           
           re = real_part(sigmavr_fft(m,*))
           im = imaginary(sigmavr_fft(m,*))
           
           sigmam_dot+= dcomplex(int_tabulated(rad(r1:r2), re(r1:r2)) ,$
                                 int_tabulated(rad(r1:r2), im(r1:r2)))
           
           re = real_part(sigmavt_fft(m,*))
           im = imaginary(sigmavt_fft(m,*))
           
           sigmam_dot+= dcomplex(int_tabulated(rad(r1:r2), re(r1:r2)) ,$
                                 int_tabulated(rad(r1:r2), im(r1:r2)))
           
           sigmam_dot *= -1d0 
           
           pattern = dcomplex(0d0, 1d0)*sigmam_dot/sigmam
           pattern = real_part(pattern)/m 
        endif else begin
           pattern = rc(i)^(-1.5)
        endelse

;        pattern = int_tabulated(rad(r1:r2),real_freq(m,r1:r2))/(rad(r2) - rad(r1))
;        pattern/= m 
;        pattern  = real_freq(m,*)/m 

        en_den = 2d0*m^2*!dpi^2*abs(sigma_m)^2*sigma1d*(pattern - omega)
        en_den /= (kappa2 - m^2*(pattern - omega)^2 + kr^2*csq)^2  
        
        mode_energy(i,k-start) = int_tabulated(rad(r1:r2), en_den(r1:r2)*2d0*!dpi*rad(r1:r2))
     endfor
     print, k
  endfor

  mode_energy(*,*) *= scale

  set_plot, 'ps'
  device, filename=filepath(strcompress('wkb_energy.ps',/remove_all),root_dir='.',subdir=[location]) $ 
          ,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches
  plot, time, mode_energy(0,*),xmargin=[8,2],ymargin=[3.2,1.8], ytitle=textoidl('J_m') $
        ,charsize=1.5, thick=4, xtitle=textoidl('t/P_0'),yrange=yrange,ystyle=1,xstyle=1
  
  for i = 1, nmodes-1 do begin
     oplot, time, mode_energy(i,*), thick=4, linestyle=i
  endfor
  
  if keyword_set(legend) then begin
     x0=legend(0)
     x1=legend(1)
     y0=legend(2)
     dy=legend(3)
     for j=0, n_elements(label)-1 do begin
        oplot, [x0,x1], [y0,y0]-dy*j, thick=4, linestyle=j
        xyouts, x1, y0-dy*j,textoidl(label(j)),charsize=1.5
     endfor
  endif

  device,/close
end

