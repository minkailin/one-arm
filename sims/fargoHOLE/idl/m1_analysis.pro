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

function get_data, loc, start, mode
  
  dims=(read_ascii(filepath('dims.dat',root_dir='.',subdir=[loc]))).(0)
  nout=fix(dims(5))
  nrad=fix(dims(6))
  nsec=fix(dims(7))
  
  radtmp=dblarr(nrad+1)
  rad   =dblarr(nrad)
  openr,1,filepath('used_rad.dat',root_dir='.',subdir=[loc])
  readf,1,radtmp
  close,1
  rad(0:nrad-1)=(radtmp(0:nrad-1)+radtmp(1:nrad))/2d0
  
  data  = dblarr(nsec,nrad)
  
  ;; initial density field 
  ;; openr,2,filepath(strcompress('gasdens'+'0'+'.dat',/remove_all),root_dir='.',subdir=[loc])
  ;; readu,2,data
  ;; close,2
  ;; sigma0 = data
  
  ;; fft_sigma0 = fft(sigma0, -1, dimension=1, /double)
  ;; data0       = abs(fft_sigma0(0,*))
  
  ;current density field 
  openr,3,filepath(strcompress('gasdens'+string(start)+'.dat',/remove_all),root_dir='.',subdir=[loc])
  readu,3,data
  close,3
  sigma = data 
 
  fft_sigma = fft(sigma, -1, dimension=1, /double)
  data1d     = fft_sigma(mode,*);/data0

  return, data1d
end

pro m1_analysis, loc=loc, start=start, finish=finish, legend=legend,label=label, xrange=xrange, yrange=yrange, ytickinterval=ytickinterval, r0=r0, xtickinterval=xtickinterval, rrange=rrange, gmma=gmma, smallh=smallh, basic=basic, azimode=azimode, rc=rc  

  !p.font = 0
  if not keyword_set(azimode) then azimode=1
  azistring = strcompress(string(azimode),/remove_all)
 
  location=strcompress(loc,/remove_all)
  
  if not keyword_set(r0) then r0 = 1.0 
  rplot = get_grid(loc, start, r0)
  rad   = rplot*r0
  xtitle=textoidl('r/r_0') 
  temp = min(abs(rplot - rrange(0)), r1)
  temp = min(abs(rplot - rrange(1)), r2)
  
  dims=(read_ascii(filepath('dims.dat',root_dir='.',subdir=[location]))).(0)
  nrad=fix(dims(6))
  nsec=fix(dims(7))
  
  nlines = file_lines(filepath('planet0.dat',root_dir='.',subdir=[loc]))
  info=dblarr(11,nlines)
  openr,3,filepath('planet0.dat',root_dir='.',subdir=[location])
  readf,3,info
  close,3
  
  p0=2d0*!dpi*(r0)^(3./2.)

  data    = dblarr(nsec,nrad)
  sigma   = dblarr(nrad) 
  vtheta  = dblarr(nrad)
  kappa2  = dblarr(nrad)
  pressure= dblarr(nrad)
  kr_plus = dblarr(nrad)
  kr_minus= dblarr(nrad)


  
  m1_density  = dblarr(nsec, nrad)
  azi         = 2d0*!dpi*dindgen(nsec)/(nsec - 1d0)
  azi        /= azimode   

  time        = dblarr(finish-start + 1)
  outnum      = dblarr(finish-start + 1)
  mode_amp    = dblarr(finish-start + 1)
  mode_phimax = dblarr(finish-start + 1)
  


  if not keyword_set(basic) then basic = start
  
;use axisymmetric component of t=0 as normalization
  openr,2,filepath(strcompress('gasdens'+string(basic)+'.dat',/remove_all),root_dir='.',subdir=[location])
  readu,2,data
  close,2
  sigma0_fft = fft(data, -1, dimension=1,/double)
  
  openr,2,filepath(strcompress('gasvtheta'+string(basic)+'.dat',/remove_all),root_dir='.',subdir=[location])
  readu,2,data
  close,2
  vtheta0_fft = fft(data, -1, dimension=1,/double)
  
  angmom0 = int_tabulated(rad, 2d0*!dpi*rad*rad*real_part(vtheta0_fft(0,*))*real_part(sigma0_fft(0,*)))
  
  for k=start, finish do begin
     if (k lt 1000) then begin
        ks=string(k,format='(I03)')
     endif else ks=string(k,format='(I04)')
     
     outnum(k-start) = k 
     time(k-start)   = info(7,k)
     
     fourier = get_data(loc, k, azimode)
     
     for j=0, nsec-1 do begin
        for i=r1, r2 do begin
           m1_density(j,i) = 2d0*real_part( fourier(i)*( cos(azimode*azi(j)) + dcomplex(0d0,1d0)*sin(azimode*azi(j))) )  
        endfor
     endfor
     
     mode_amp(k-start)    = max(alog(abs(fourier(r1:r2))))
     
     temp        = max(m1_density, grid)
     result      = array_indices(m1_density, grid)
     mode_phimax(k-start) = azi(result(0))
     
        
     re  = real_part(fourier)
     dre = deriv(rad, re)
     
     im  = imaginary(fourier)
     dim = deriv(rad, im)

     kr  = dim*re - dre*im
     kr /= re^2 + im^2 + 1d-16
     
     openr,2,filepath(strcompress('gasdens'+string(k)+'.dat',/remove_all),root_dir='.',subdir=[location]) 
     readu,2,data
     close,2
     for i=0, nrad-1 do sigma(i) = mean(data(*,i))

  
     openr,2,filepath(strcompress('gasvtheta'+string(k)+'.dat',/remove_all),root_dir='.',subdir=[location]) 
     readu,2,data
     close,2
     for i=0, nrad-1 do vtheta(i) = mean(data(*,i))
     omega = vtheta/rad
     
     kappa2 = epicyclesq(rad, vtheta)
  
     
     if keyword_set(gmma) then begin
        openr,3,filepath(strcompress('gasTemperature'+string(k)+'.dat',/remove_all),root_dir='.',subdir=[loc])
        readu,3,data
        close,3
        for i=0, nrad-1 do pressure(i) = mean(data(*,i))*sigma(i)
     endif else begin
        pressure = (smallh*smallh/rad)*sigma
     endelse
     
     Q = toomreQ(rad, vtheta, sigma, pressure)


     kr_norm = kr/(kappa2/(2d0*!dpi*sigma))
  

     if keyword_set(rc) then begin
        omega_mode=azimode*rc^(-1.5) ;approx rotation curve as keplerian for 
        nu = omega_mode-azimode*omega
        nu/= sqrt(kappa2)
        for i=0, nrad-1 do begin
           qb = 1d0-Q(i)^2*(1.0-nu(i)^2)
           if(qb ge 0d0) then begin
              kr_plus(i)  = (2d0/Q(i))*(1d0 + sqrt(qb))
              kr_minus(i) = (2d0/Q(i))*(1d0 - sqrt(qb))
           endif else begin
              kr_plus(i) = 0d0
              kr_minus(i)= 0d0
           endelse
        endfor

        set_plot, 'ps'
        device, filename=filepath(strcompress('m1_analysis_kr'+ks+'.ps',/remove_all),root_dir='.',subdir=[location]) $ 
                ,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches
        plot, rplot, kr_plus,xmargin=[8,2],ymargin=[3.2,1.8], ytitle=textoidl('|k(r)|/k_T(r)'), ystyle=0,xrange=xrange,yrange=yrange,xstyle=1  $
              ,charsize=1.5, thick=4, xtitle=textoidl('r/r_0'), title = strcompress(string(info(7,k)/p0,format='(F7.2)')+textoidl('P_0'),/remove_all) +', m='+azistring
        oplot, rplot, kr_minus, thick=4, linestyle=1
        oplot, rplot, kr_norm, thick=4,linestyle=2
        device,/close

        ;theoretical angular mometum
        

        jdens=2d0*azimode*!dpi^2*sigma*(omega_mode - azimode*omega)*abs(fourier)^2
        jdens/=(kappa2-(omega_mode - azimode*omega)^2+kr^2*pressure/sigma)^2

        ang= int_tabulated(rad(r1:r2), 2d0*!dpi*rad(r1:r2)*jdens(r1:r2))

 
     endif
     
     print, k;, ang/angmom0
  endfor
  
                                ;plot w.r.t. output number to select
                                ;appropriate range of outputs (linear growth)
  set_plot, 'ps'
  device, filename=filepath(strcompress('m1_analysis_amp.ps',/remove_all),root_dir='.',subdir=[location]) $ 
          ,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches
  plot, outnum, mode_amp,xmargin=[8,2],ymargin=[3.2,1.8], ytitle=textoidl('max(ln|\Sigma_m/\Sigma_0|)'), ystyle=0  $
        ,charsize=1.5, thick=4, xtitle='Output number', title = 'm='+azistring 
  device,/close
  
  set_plot, 'ps'
  device, filename=filepath(strcompress('m1_analysis_phimax.ps',/remove_all),root_dir='.',subdir=[location]) $ 
          ,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches
  plot, outnum, mode_phimax/(2d0*!dpi/azimode),xmargin=[8,2],ymargin=[3.2,1.8], ytitle=textoidl('\phi_{max}/(2\pi/m)'), ystyle=0  $
        ,charsize=1.5, thick=4, xtitle='Output number' , title = 'm='+azistring
  device,/close
  
;linear fits to amplitude growth and 
  Result = LINFIT(time, mode_amp, /DOUBLE) 
  rate   = result(1)

  if not keyword_set(rc) then begin  
  Result = LINFIT(time, mode_phimax, /DOUBLE) 
  pattern = result(1)
  rc = pattern^(-2./3.) ; approx rotation curve as kepler
  endif else begin
  pattern = rc^(-1.5)
  endelse

  print, 'rate, rc', rate, rc

;dispersion relation plots
  real_freq = pattern*azimode
  nu2 = (real_freq - azimode*omega)^2/kappa2
  
  set_plot, 'ps'
  device, filename=filepath(strcompress('m1_analysis_LR.ps',/remove_all),root_dir='.',subdir=[location]) $ 
          ,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches
  plot, rplot, nu2 - 1,xmargin=[8,2],ymargin=[3.2,1.8], ytitle=textoidl('\nu^2 - 1'), ystyle=0  $
        ,charsize=1.5, thick=4, xtitle=textoidl('r/r_0'), xstyle=1, title = 'm='+azistring, yrange=yrange
  oplot, [0,1]*max(rad), [0,0], linestyle=1
  device,/close
  
  set_plot, 'ps'
  device, filename=filepath(strcompress('m1_analysis_Qbar.ps',/remove_all),root_dir='.',subdir=[location]) $ 
          ,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches
  plot, rplot, nu2 - 1 + 1d0/Q^2,xmargin=[8,2],ymargin=[3.2,1.8], ytitle=textoidl('\nu^2 - 1 + Q^{-2}'), ystyle=0  $
        ,charsize=1.5, thick=4, xtitle=textoidl('r/r_0'), xstyle=1, title = 'm='+azistring, yrange=yrange
  oplot, [0,1]*max(rad), [0,0], linestyle=1
  device,/close

  ;; cs      = soundspeed(pressure,sigma)
  ;; cg_norm = kr*cs - sqrt(kappa2)/Q ;this kr is that of output k=finish
  ;; cg_norm/= azimode*(pattern - omega)
  
  ;; set_plot, 'ps'
  ;; device, filename=filepath(strcompress('m1_analysis_cg'+ks+'.ps',/remove_all),root_dir='.',subdir=[location]) $ 
  ;;         ,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches
  ;; plot, rplot, cg_norm,xmargin=[8,2],ymargin=[3.2,1.8], ytitle=textoidl('c_g/c_s'), ystyle=0  $
  ;;       ,charsize=1.5, thick=4, xtitle=textoidl('r/r_0'), xstyle=1,xrange=xrange, title = 'm='+azistring 
  ;; oplot, [0,1]*max(rad), [0,0], linestyle=1
  ;; device,/close
  


end

