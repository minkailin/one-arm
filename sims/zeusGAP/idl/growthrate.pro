function hdf, n, name
  
  fileid = hdf_sd_start(name,/read)
  sds = hdf_sd_select(fileid,n)
  hdf_sd_getdata,sds,data
  
  return, data
end

pro growthrate, loc=loc, basic=basic, start=start, finish=finish, azimode=azimode, ravg=ravg
  !p.font=0
  ii     = dcomplex(0d0, 1d0) 
  nmodes = 10

  if not keyword_set(finish) then finish=start 
; estimating growth rates of a mode from a single snapshot

;get the basic info
  location =strcompress(loc,/remove_all)
  if not keyword_set(basic) then begin
     filename = strcompress('hdfaa.000',/remove_all)
  endif else begin
     name = string(basic,format='(I03)')
     filename = strcompress('hdfaa.'+name,/remove_all)
  endelse
  fileloc  = filepath(filename,root_dir='.',subdir=[location])

  rad  = hdf(2, fileloc)
  theta= hdf(1, fileloc)
  phi  = hdf(0, fileloc)
  
  nrad  = n_elements(rad)
  ntheta= n_elements(theta)
  nphi  = n_elements(phi)
  
  v0   =  hdf(11, fileloc)
  data_axisymmetric = dblarr(nrad, ntheta)  
  
  for j=0, ntheta-1 do begin
     for i=0, nrad-1 do begin
        data_axisymmetric(i,j) = mean(v0(i,j,*))
     endfor
  endfor
  
 

  for k=0, nphi-1 do begin
     v0(*,*,k) = data_axisymmetric(*,*)
  endfor
 
; assume uniform theta and phi spacing 
  dtheta = theta(1) - theta(0)
  dphi = phi(1) - phi(0)
  
; planet info
  nlines = file_lines(filepath('planetxy_hdf.dat',root_dir='.',subdir=[location]))
  planetinfo = dblarr(7,nlines)
  openr, 1, filepath('planetxy_hdf.dat',root_dir='.',subdir=[location])
  readf,1,planetinfo
  close,1
  torb = 2d0*!dpi*planetinfo(1,0)^(3d0/2d0)
  
;normalization for fourier density amplitudes 
  n=0
  plx=planetinfo(1,n)
  ply=planetinfo(2,n)
  plz=planetinfo(3,n)
  phiplt = pltphi(plx,ply)
  plrad = sqrt(plx^2 + ply^2 + plz^2)
  
  if keyword_set(mp) then begin
     rhill = plrad*(mp/3d0)^(1d0/3d0)
     rplot = (rad - plrad)/rhill
  endif else begin
     rplot = rad/plrad
  endelse
  temp = min(abs(rplot-1d0), r0)
  temp = min(abs(rplot-ravg(0)), r1)
  temp = min(abs(rplot-ravg(1)), r2)
  nrnew = r2-r1+1 

  den        = hdf(19, fileloc)
  den_fft    = fft(den(r1:r2,*,*), -1, dimension=3, /double)
  dfft_norm  = real_part(den_fft(*,*,0)) ;axis. component is real
  norm1      = dblarr(nrnew)

  for i=0, nrnew-1 do begin
     norm1(i) = int_tabulated(theta,  dfft_norm(i,*))
  endfor
  norm1 /= theta(ntheta-1) - theta(0)
  
  norm = int_tabulated(rad(r1:r2), norm1)
  norm/= rad(r2) - rad(r1)

  time     = dblarr(finish-start + 1)
  mode_amp = dblarr(finish-start + 1)
 
  for n=start, finish do begin
     time(n-start) = planetinfo(0,n)
     ks   = string(n,format='(I03)')
     filename = strcompress('hdfaa.'+ks,/remove_all)
     fileloc  = filepath(filename,root_dir='.',subdir=[location])
     
     
;density field
     den   = hdf(19, fileloc)
;momentum 
     vrad  = hdf(3, fileloc)*den
     vtheta= hdf(7, fileloc)*den
     vphi  = hdf(11, fileloc)*den
      
     plx=planetinfo(1,n)
     ply=planetinfo(2,n)
     plz=planetinfo(3,n)
     phiplt = pltphi(plx,ply)
     plrad = sqrt(plx^2 + ply^2 + plz^2)
     
     if keyword_set(mp) then begin
        rhill = plrad*(mp/3d0)^(1d0/3d0)
        rplot = (rad - plrad)/rhill
     endif else begin
        rplot = rad/plrad
     endelse
     temp = min(abs(rplot-1d0), r0)
     temp = min(abs(rplot-ravg(0)), r1)
     temp = min(abs(rplot-ravg(1)), r2)
     nrnew = r2-r1+1    
     radsmall    = dblarr(nrnew)
     
     dfft_2d     = dcomplexarr(r2-r1+1, ntheta)
     vradfft_2d  = dcomplexarr(r2-r1+1, ntheta)
     vthetafft_2d= dcomplexarr(r2-r1+1, ntheta)
     vphifft_2d  = dcomplexarr(r2-r1+1, ntheta)
     
     dvrad_dr           = dcomplexarr(r2-r1+1, ntheta)
     dvtheta_dtheta     = dcomplexarr(r2-r1+1, ntheta)
     
     
     den_fft    = fft(den(r1:r2,*,*), -1, dimension=3, /double)
     vrad_fft   = fft(vrad(r1:r2,*,*), -1, dimension=3, /double)
     vtheta_fft = fft(vtheta(r1:r2,*,*), -1, dimension=3, /double)
     vphi_fft   = fft(vphi(r1:r2,*,*), -1, dimension=3, /double)
     
     dfft_2d     = den_fft(*,*,azimode)
    
     vradfft_2d  = vrad_fft(*,*,azimode)
     vthetafft_2d= vtheta_fft(*,*,azimode)
     vphifft_2d  = vphi_fft(*,*,azimode)
     
;radial derivative of radial momentum
     radsmall(*) = rad(r1:r2)
     for j=0, ntheta-1 do begin
        vradfft_2d(*,j) *= radsmall(*)^2
        dvrad_dr(*,j)  = dcomplex( deriv(radsmall, real_part(vradfft_2d(*,j))), $
                                   deriv(radsmall, imaginary(vradfft_2d(*,j))) )    
        dvrad_dr(*,j) /= radsmall(*)^2
     endfor
;theta derivative of polar momentum
     for i=0, nrnew-1 do begin
        vthetafft_2d(i,*)     *= sin(theta(*))
        dvtheta_dtheta(i,*)  = dcomplex( deriv(theta, real_part(vthetafft_2d(i,*))), $
                                         deriv(theta, imaginary(vthetafft_2d(i,*))) )
        dvtheta_dtheta(i,*) /= radsmall(i)*sin(theta(*))
     endfor
;phi derivative of azimuthal momentum, just multiply by im
     vphifft_2d *= ii*azimode
     for j=0, ntheta-1 do begin
        vphifft_2d(*,j) /= radsmall*sin(theta(j))
     endfor
     
     rhodot = -(dvrad_dr + dvtheta_dtheta + vphifft_2d)
     complex_rate = rhodot/dfft_2d
     
     real_rate = real_part(complex_rate)
     imag_rate = -imaginary(complex_rate)    
     
     avg1 = dblarr(nrnew)
     for i=0, nrnew-1 do avg1(i) = int_tabulated(theta, real_rate(i,*))
     avg1 /= theta(ntheta-1) - theta(0)
     
     avg_growth = int_tabulated(radsmall, avg1)
     avg_growth/= rad(r2) - rad(r1)
     avg_growth/= (v0(r0,ntheta-1,0)/plrad)
     
     for i=0, nrnew-1 do avg1(i) = int_tabulated(theta, imag_rate(i,*))
     avg1 /= theta(ntheta-1) - theta(0)
     
     avg_freq = int_tabulated(radsmall, avg1)
     avg_freq/= rad(r2) - rad(r1)

     temp = min(abs(azimode*v0(*,ntheta-1,0)/rad(*) - avg_freq), grid)
     if( (grid eq nrad-1) or (grid eq 0) ) then begin
        rc = (avg_freq/azimode)^(-2d0/3d0)
     endif else rc = rad(grid)
     
     print, time(n-start)/torb, rc, avg_growth

;     avg_patt   = mean(imag_rate(1:nrnew-2,1:ntheta-2))/(azimode*v0(r0,ntheta-1,0)/plrad)
;     print, '<growth rate>/omega0=', avg_growth
;     print, '<pattern speed>/m.omega0=', avg_patt
;     print, 'average mode amp', mean(abs(dfft_2d))
     

;mode amplitudes, and average them
     real = real_part(dfft_2d)
     imag = imaginary(dfft_2d)
     
     norm1_re      = dblarr(nrnew)
     norm1_im      = dblarr(nrnew)

     for i=0, nrnew-1 do begin
        norm1_re(i) = int_tabulated(theta,  real(i,*))
        norm1_im(i) = int_tabulated(theta,  imag(i,*))
     endfor
     norm1_re /= theta(ntheta-1) - theta(0)
     norm1_im /= theta(ntheta-1) - theta(0)
    
     norm_re  = int_tabulated(rad(r1:r2), norm1_re)
     norm_im  = int_tabulated(rad(r1:r2), norm1_im)

     norm_re /= rad(r2) - rad(r1)
     norm_im /= rad(r2) - rad(r1)

     amp = sqrt(norm_re^2 + norm_im^2)

     time(n-start)     = planetinfo(0,n)
     mode_amp(n-start) = amp/norm
  endfor

;explicit rate cal.
  growth_ex = deriv(time, alog(abs(mode_amp)))
  for i=0, finish-start do begin
     print, time(i)/torb, growth_ex(i)
  endfor


  set_plot, 'ps'
  device, filename=filepath('growthrate_.ps',root_dir='.',subdir=[location]) $
          , xsize=8, ysize=4.5, xoffset=0, yoffset=0 $
          , /inches,/color,bits_per_pixel=8
  plot, time/torb, alog10(abs(mode_amp)), xmargin=[6,2], ymargin=[3.5,0.5] , ytitle=textoidl('log_{10}(C_m/C_0)') $
        , charsize=1.5,thick=4,linestyle=0, xtitle=textoidl('t/P_0'), xrange=xrange, yrange=yrange, xstyle=1
  device,/close
  
end
