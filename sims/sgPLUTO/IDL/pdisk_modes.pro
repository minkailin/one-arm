function radavg, rad, data, r1, r2, t1, t2
  data_re = real_part(data)
  data_im = imaginary(data)
  
  res_re = mean(data_re(r1:r2, t1:t2))
  res_im = mean(data_im(r1:r2, t1:t2))
  denom  = 1.0
  
  return, dcomplex(res_re, res_im)/denom
end

pro pdisk_modes, start=start, finish=finish, r0=r0, o8=o8, o10=o10, omit=omit, zslice=zslice, $
                 legend=legend, label=label, xtickinterval=xtickinterval, ytickinterval=ytickinterval $
                 ,xrange=xrange, yrange=yrange, range=range, zeronorm=zeronorm 
  
  COMMON PLUTO_GRID,  nx1,nx2,nx3,x1,x2,x3,$
     dx1,dx2,dx3, geometry, xpos, ypos,$
     AMRLevel, AMRBoxes         ;  ** Chombo data structure **
                                ;  ** loaded when HDF5LOAD is called **
  
  COMMON PLUTO_VAR, NVAR, rho, vx1, vx2, vx3, $
     bx1, bx2, bx3, $
     Ax1, Ax2, Ax3, $
     bx1s, bx2s, bx3s,pot,$
                                ; ----------------------------------------
     v1, v2, v3, $              ; Kept for backward
     b1, b2, b3, $              ; compatibility with 
     A1, A2, A3, $              ; PLUTO 3 data name
     b1s, b2s, b3s, $           ;
     pr,            $           ;
                                ; -----------------------------------------
     prs, psi_glm, fneut, fion, fh2, $
     tmp, tr1, tr2, tr3, tr4,$
     fhI, fheI, fheII,$
     fcI, fcII, fcIII, fcIV, fcV,$
     fnI, fnII, fnIII, fnIV, fnV,$
     foI, foII, foIII, foIV, foV,$
     fneI, fneII, fneIII, fneIV, fneV,$
     fsI, fsII, fsIII, fsIV, fsV, vars, vname
  COMMON PLUTO_RUN,  t, dt, nlast, first_call
  
  !p.font = 0

  nmodes = 6
  nmodes_p1 = nmodes+1
  fmtstring = '('+string(nmodes_p1,format='(I2)')+'(e15.8,x))'
  
  if not keyword_set(finish) then finish=start 
  if not keyword_set(r0) then r0=1.0
  
  pload, 0
  
  nrad  = nx1
  ntheta= nx2
  nphi  = nx3
  
  rad   = x1
  theta = x2
  phi   = x3
  time  = t
  
  dr     = rad(1) - rad(0)      ; assume uniform
  dtheta = theta(1) - theta(0)  ; assume uniform
  
;define averaging zone 
  temp = min(abs(rad - r0), rzero)
if not keyword_set(range) then begin
  temp = min(abs(rad -0.8*r0), r1)
  temp = min(abs(rad -1.2*r0), r2)
endif else begin
  temp = min(abs(rad -range(0)*r0), r1)
  temp = min(abs(rad -range(1)*r0), r2)
endelse  

  if not keyword_set (omit) then omit = 5
  t2   = ntheta-1 - fix(omit)
  if keyword_set(zslice) then begin
     t1   = t2
  endif else begin
     t1   =  omit
  endelse

  ;; t1 = ntheta-1
  ;; t2 = ntheta-1
 
  omega_zero = vx3(rzero,ntheta-1,0)/rad(rzero)
  omega1d    = vx3(*,ntheta-1,0)/rad(*)

;work arrays
  fft_dens      = dcomplexarr(nrad, ntheta, nphi)
  fft_dens_dot  = dcomplexarr(nrad, ntheta, nphi)
  fft_mom_rad   = dcomplexarr(nrad, ntheta, nphi)
  fft_mom_theta = dcomplexarr(nrad, ntheta, nphi)
  fft_mom_phi   = dcomplexarr(nrad, ntheta, nphi)
  
  work_theta_1d = dcomplexarr(ntheta)
  work_rad_1d   = dcomplexarr(nrad)
  sigma         = dcomplexarr(nrad, ntheta)
  
;result arrays
  amp_2d       = dblarr(nrad, ntheta, nmodes+1, finish-start+1)
  avg_amp      = dblarr(nmodes+1,finish-start+1) ;average of (amplitude of(complex amp))
  
  growth_rates = dblarr(nmodes+1,finish-start+1)
  mode_freqs   = dblarr(nmodes+1,finish-start+1) 
  mode_corot   = dblarr(nmodes+1,finish-start+1)  
; amplitudes normalized wrt initial magnitude of m=0 mode
  fft_dens = fft(rho, -1, dimension=3, /double)
  res      = radavg(rad, abs(fft_dens(*,*,0)), r1, r2, t1, t2)
  amp_zero = real_part(res)
  
;time loop
  for n=start, finish do begin
     
     avg_amp(  0, n-start)    = t(n)
     growth_rates(0, n-start) = t(n)
     mode_freqs(  0, n-start) = t(n)
     mode_corot(  0, n-start) = t(n)
;read data
     pload, n,/silent

;fft the density field in azimth
     fft_dens = fft(rho, -1, dimension=3, /double)
     
;get normalized mode amplitudes
     if keyword_set(zeronorm) then begin;then we are normalizing wrt current time 
     res      = radavg(rad, abs(fft_dens(*,*,0)), r1, r2, t1, t2)
     amp_zero = real_part(res)
     omega_zero = vx3(rzero,ntheta-1,0)/rad(rzero)
     omega1d = vx3(*,ntheta-1,0)/rad(*)
     endif
     for m=1, nmodes do begin
        amp_2d(*,*,m,n-start) = abs(fft_dens(*,*,m))
        avg_amp(m, n-start) = real_part(radavg(rad, amp_2d(*,*,m,n-start), r1, r2, t1, t2))/amp_zero
     endfor

;compute complex frequencies associated with m-th mode

     rho_vr     = rho*vx1
     rho_vtheta = rho*vx2
     rho_vphi   = rho*vx3
     
     fft_mom_rad   = fft(rho_vr, -1, dimension=3, /double)
     fft_mom_theta = fft(rho_vtheta, -1, dimension=3, /double)
     fft_mom_phi   = fft(rho_vphi, -1, dimension=3, /double)
     
     for m=1, nmodes do begin
        
;radial contribution to diverence
        for j=t1, t2 do begin
           work_rad_1d = fft_mom_rad(*,j,m)
           if ( (not keyword_set(o8)) and (not keyword_set(o10)) ) then begin
              fft_mom_rad(*,j,m) = dcomplex(deriv(rad,real_part(work_rad_1d)), $
                                            deriv(rad,imaginary(work_rad_1d)) )  
           endif
           for i=r1, r2 do begin
;8th order FD
              if keyword_set(o8) then begin
                 fft_mom_rad(i,j,m) = (1.0/280.)*work_rad_1d(i-4) - (4.0/105.0)*work_rad_1d(i-3) $
                                      +(1.0/5.0)*work_rad_1d(i-2) - (4.0/5.0)*work_rad_1d(i-1) $
                                      +(4.0/5.0)*work_rad_1d(i+1) - (1.0/5.0)*work_rad_1d(i+2) $
                                      +(4.0/105.0)*work_rad_1d(i+3) -(1.0/280.)*work_rad_1d(i+4)
                 fft_mom_rad(i,j,m)/= dr
              endif
              if keyword_set(o10) then begin
;10th order FD
                 fft_mom_rad(i,j,m) = -2.0*work_rad_1d(i-5) + 25.0*work_rad_1d(i-4) -150.0*work_rad_1d(i-3) $
                                      + 600.0*work_rad_1d(i-2) - 2100.0*work_rad_1d(i-1) $
                                      +2100.0*work_rad_1d(i+1) - 600.0*work_rad_1d(i+2) $
                                      +150.0*work_rad_1d(i+3) -25.0*work_rad_1d(i+4) +2.0*work_rad_1d(i+5)
                 fft_mom_rad(i,j,m)/= 2520.0*dr
              endif
           endfor
           fft_mom_rad(*,j,m) += 2.0*work_rad_1d/rad 
        endfor
        
;theta contribution to divergence        
        for i=r1, r2 do begin
           work_theta_1d = fft_mom_theta(i,*,m)
           if ( (not keyword_set(o8)) and (not keyword_set(o10)) ) then begin
              fft_mom_theta(i,*,m) = dcomplex(deriv(theta, real_part(work_theta_1d)), $
                                              deriv(theta, imaginary(work_theta_1d) ))
           endif
           for j=t1, t2 do begin
              if keyword_set(o8) then begin
;8th order FD
                 fft_mom_theta(i,j,m) = (1.0/280.)*work_theta_1d(j-4) - (4.0/105.0)*work_theta_1d(j-3) $
                                        +(1.0/5.0)*work_theta_1d(j-2) - (4.0/5.0)*work_theta_1d(j-1) $
                                        +(4.0/5.0)*work_theta_1d(j+1) - (1.0/5.0)*work_theta_1d(j+2) $
                                        +(4.0/105.0)*work_theta_1d(j+3) -(1.0/280.)*work_theta_1d(j+4)
                 fft_mom_theta(i,j,m)/= dtheta
              endif
              if keyword_set(o10) then begin
;10th order FD
                 fft_mom_theta(i,j,m) = -2.0*work_theta_1d(j-5) + 25.0*work_theta_1d(j-4) -150.0*work_theta_1d(j-3) $
                                        + 600.0*work_theta_1d(j-2) - 2100.0*work_theta_1d(j-1) $
                                        +2100.0*work_theta_1d(j+1) - 600.0*work_theta_1d(j+2) $
                                        +150.0*work_theta_1d(j+3) -25.0*work_theta_1d(j+4) +2.0*work_theta_1d(j+5)
                 fft_mom_theta(i,j,m)/= 2520.0*dtheta
              endif
           endfor 
           fft_mom_theta(i,*,m) += cos(theta)*work_theta_1d/sin(theta)
           fft_mom_theta(i,*,m) /=rad(i)
        endfor
        
;azimuthal contribution to divergence        
        for j=t1, t2 do begin
           fft_mom_phi(r1:r2,j,m) *= dcomplex(0.0,1.0)*m/(rad(r1:r2)*sin(theta(j)))
        endfor
        
;construct divergence
        fft_dens_dot(r1:r2,t1:t2,m) = - (fft_mom_rad(r1:r2,t1:t2,m) + fft_mom_theta(r1:r2,t1:t2,m) + fft_mom_phi(r1:r2,t1:t2,m))
        
;get spatially-dependent complex frequency then average it over space       
        if(min(abs(fft_dens(r1:r2,t1:t2,m))) eq 0.0) then begin
           freq = 0.0
        endif else begin
           sigma(r1:r2,t1:t2) = dcomplex(0.0,1.0)*fft_dens_dot(r1:r2,t1:t2,m)/fft_dens(r1:r2,t1:t2,m)
           freq = radavg(rad, sigma(*,*), r1, r2, t1, t2)
        endelse
        growth_rates(m, n-start) = imaginary(freq)
        growth_rates(m, n-start)/= omega_zero
        
        mode_freqs(m,   n-start)  = real_part(freq)
        mode_freqs(m,   n-start) /= m*omega_zero 

        temp = min(abs(real_part(freq) - m*omega1d(r1:r2)), grid)
        mode_corot(m, n-start) = (real_part(freq)/m)^(-2d0/3d0);rad(r1+grid)

     endfor
  endfor
  
  openw,1, 'pdisk_modes_amp.dat'
  for i=0, finish-start do begin
     printf, 1, avg_amp(0:nmodes,i), format=fmtstring
  endfor
  close,1

;if (n_elements(avg_amp(0, *)) ge 10) then begin
;;find largest mode at t=10 (assume this is 10th output) 
;  temp = max(growth_rates(3:nmodes, 10),mmax);first slot is time
;  mmax += 3
;  print, 'fastest growing mode at t=10 is m=', mmax, ' with amplitude', avg_amp(mmax, 10)
;endif;

;plot the mode amp as function of time 
  avg_amp(0,*) /= (2.0*!dpi/omega_zero)
  if not keyword_set(zeronorm) then begin
  ytitle = strcompress(textoidl('ln(a_m)'),/remove_all)
  endif else begin
  ytitle = strcompress(textoidl('ln(a_m/a_0)'),/remove_all)
  endelse  
 
  fname = 'pdisk_modes_amp.ps'
  set_plot, 'ps'
  avg_amp(1:nmodes,1:finish-start) = alog(avg_amp(1:nmodes,1:finish-start))
  device, filename=fname $
          ,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches
  plot, avg_amp(0,1:finish-start), avg_amp(1,1:finish-start),xmargin=[7.5,1.5],ymargin=[3.5,0.5], ystyle=0  $
        ,charsize=1.5, thick=4, xrange=xrange, xtitle=textoidl('t/P_0'), yrange=yrange, $ ;[min(avg_amp(1:nmodes,*)),max(avg_amp(1:nmodes,*))],$ 
        linestyle = 0, ytitle = ytitle, xtickinterval=xtickinterval, ytickinterval=ytickinterval
  for k=2, nmodes do begin
  oplot, avg_amp(0,1:finish-start), avg_amp(k,1:finish-start), thick=4, linestyle=k-1
  endfor
  if keyword_set(legend) then begin
   x0=legend(0)
   x1=legend(1)
   y0=legend(2)
   dy=legend(3)
   for k=0, n_elements(label)-1 do begin
      oplot, [x0,x1], [y0,y0]-dy*k, thick=4, linestyle=k
      xyouts, x1, y0-dy*k,textoidl(label(k)),charsize=1.5
   endfor
  endif
  device,/close

  openw,1, 'pdisk_modes_freqs.dat'
  for i=0, finish-start do begin
     printf, 1, mode_freqs(0:nmodes,i), format=fmtstring
  endfor
  close,1
  

  openw,1, 'pdisk_modes_corot.dat'
  for i=0, finish-start do begin
     printf, 1, mode_corot(0:nmodes,i), format=fmtstring
  endfor
  close,1



;; if (n_elements(avg_amp(0, *)) ge 10) then begin
;;  ;freq of max mode, averaged over 5 orbits
;;   arr =  mode_freqs(mmax,5:10)
;;   om  = mean(arr)
;;   print, 'time average of freq, dev =', om, STDDEV(arr,/DOUBLE)/om
;; endif  

   openw,1, 'pdisk_modes_rates.dat'
   for i=0, finish-start do begin
      printf, 1, growth_rates(0:nmodes,i), format=fmtstring
   endfor
   close,1

;; if (n_elements(avg_amp(0, *)) ge 10) then begin
;; ;rate of max mode, averaged over 5 orbits
;;   arr =  growth_rates(mmax,5:10)
;;   qm  = mean(arr)
;;   err = STDDEV(arr,/DOUBLE)/qm
;;   print, 'time average of rate, dev =', qm, err
;; endif
    

;explicit time differentiation of complex amplitudes to get rates
;get spatially-dependent rate, then average 
  tt = growth_rates(0,*)
  rates = dblarr(nrad, ntheta,finish-start+1)
  
  for m=1, nmodes do begin
     growth_rates(m,0) = 0.0
     for j=t1, t2 do begin
        for i=r1, r2 do begin
           rates(i,j,1:finish-start) = deriv(tt(1:finish-start), alog(amp_2d(i,j,m,1:finish-start)))
        endfor
     endfor
     
     for k=1, finish-start do begin
        growth_rates(m,k) = real_part(radavg(rad, rates(*,*,k), r1, r2, t1, t2))
     endfor
     
     growth_rates(m,1:finish-start)/=omega_zero
  endfor
  
  openw,1, 'pdisk_modes_rates_ex.dat'
  for i=0, finish-start do begin
     printf, 1, growth_rates(0:nmodes,i), format=fmtstring
  endfor
  close,1

;; if (n_elements(avg_amp(0, *)) ge 10) then begin
;; ;rate of max mode, averaged over 10 orbits
;;   arr =  growth_rates(mmax,5:10)
;;   qm_ex  = mean(arr)
;;   err_ex = STDDEV(arr,/DOUBLE)/qm_ex
;;   print, 'time average of rate_ex, dev =', qm_ex, err_ex
;; endif

;; if (n_elements(avg_amp(0, *)) ge 10) then begin
;; ;relative difference in growth rates between implicit and explicit computation
;;   dqm = abs(qm - qm_ex)
;;   if(err_ex lt err) then begin
;;      dqm/=qm_ex
;;   endif else begin
;;      dqm/=qm
;;   endelse
;;   print, 'diff. `between the qm', dqm
;; endif


;find largest mode at t=tend
  temp = max(avg_amp(1:nmodes, finish-start ),mmax_end);first slot is time
  mmax_end += 1
  print, 'max(am) at tend is m=', mmax_end, exp(avg_amp(mmax_end, finish-start))


  
end
