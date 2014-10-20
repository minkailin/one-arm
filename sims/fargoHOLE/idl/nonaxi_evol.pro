pro nonaxi_evol,  r0=r0, loc=loc, mmax=mmax, basic=basic, azimodes=azimodes, start=start, finish=finish, $
                  xrange=xrange, yrange=yrange, rrange=rrange, average=average, max=max, legend=legend, $
                  label=label, smallh=smallh 
  
  !p.font=0
  
  if not keyword_set(smallh) then smallh=0.05
  if not keyword_set(basic) then basic = 0 
  if not keyword_set(azimodes) then begin
     if not keyword_set(mmax) then mmax=4
     azimodes = dindgen(mmax) + 1d0
  endif
  
  nmodes = n_elements(azimodes) 
  
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

  phi  = 2d0*!dpi*dindgen(nsec)/(nsec - 1d0) 
  dphi = phi(1) - phi(0)

  radius_2D = dblarr(nsec,nrad)
  dVol_1D = dblarr(nrad)
  dVol_2D = dblarr(nsec,nrad)

  for i=0, nrad-1 do begin
     dVol_1D(i) = 2d0*!dpi*radius(i)*(radtmp(i+1)-radtmp(i))
  endfor
  for j=0, nsec-1 do begin
     for i=0, nrad-1 do begin
        dVol_2D(j,i) = radius(i)*(radtmp(i+1)-radtmp(i))*dphi
        radius_2D(j,i) = radius(i)
     endfor
  endfor
  
  if keyword_set(rrange) then begin
  temp = min(abs(radius - rrange(0)),r1)
  temp = min(abs(radius - rrange(1)),r2)
  endif else begin
  r1 = 0
  r2 = nrad-1
  endelse
;;;;;;;;;;
;GET DATA;
;;;;;;;;;;
  data   = dblarr(nsec,nrad)
  data1d = dblarr(nrad) 

  time_norm = dblarr(finish-start+1) 

  mode_amp  = dblarr(nmodes,finish-start+1)
  mode_ang  = dblarr(nmodes,finish-start+1)

  min_q = dblarr(finish-start+1)
  ring_mass= dblarr(finish-start+1)

  tot_ang = dblarr(finish-start+1)
  
;use axisymmetric component of t=0 as normalization
  openr,2,filepath(strcompress('gasdens'+string(basic)+'.dat',/remove_all),root_dir='.',subdir=[location]) 
  readu,2,data
  close,2
  sigma0 = data
  sigma0_fft = fft(sigma0, -1, dimension=1,/double)
  czero = total(real_part(sigma0_fft(0,*))*dVol_1D)
 

  openr,2,filepath(strcompress('gasvtheta'+string(basic)+'.dat',/remove_all),root_dir='.',subdir=[location])
  readu,2,data
  close,2
  vtheta0 = data 
  vtheta0_fft = fft(vtheta0, -1, dimension=1,/double)
  angmom0 = total(real_part(vtheta0_fft(0,r1:r2))*real_part(sigma0_fft(0,r1:r2))*radius(r1:r2)*dVol_1D(r1:r2))

;  vphi0 = data
;  for i = 0, nsec-2 do vphi0(i,*) = (data(i,*) + data(i+1,*))/2.0
;  vphi0(nsec-1,*) = (data(nsec-1,*) + data(0,*))/2.0 
;  sigma0_vtheta0 = sigma0*vphi0
;  angmom0 = total(sigma0_vtheta0*radius_2D*dVol_2D)


  for k=start, finish do begin
;time array
     time_norm(k-start) = info(7,k)/p0
     
;read data
     if (k lt 1000) then begin
        ks=string(k,format='(I03)')  
     endif else ks=string(k,format='(I04)')
     
     openr,2,filepath(strcompress('gasdens'+string(k)+'.dat',/remove_all),root_dir='.',subdir=[location]) 
     readu,2,data
     close,2
     sigma=data
     
     openr,2,filepath(strcompress('gasvtheta'+string(k)+'.dat',/remove_all),root_dir='.',subdir=[location]) 
     readu,2,data
     close,2
     vtheta=data                ;staggered 
     vphi = data                ;centered with surface density 
     for i = 0, nsec-2 do vphi(i,*) = (data(i,*) + data(i+1,*))/2.0
     vphi(nsec-1,*) = (data(nsec-1,*) + data(0,*))/2.0 
     
     
;fft      
     sigma_fft    = fft(sigma, -1, dimension=1,/double)
     vtheta_fft   = fft(vtheta, -1, dimension=1,/double)
     
     for j=0, nmodes-1 do begin
        m = azimodes(j)
;amplitude of density fourier components
        re = total(real_part(sigma_fft(m,*))*dVol_1D)
        im = total(imaginary(sigma_fft(m,*))*dVol_1D)
        amp= dcomplex(re,im)

        if not keyword_set(max) then begin
           mode_amp(j, k-start) = alog10( abs(amp)/czero + 1d-16)
        endif else begin
           mode_amp(j, k-start) = max(alog10( abs(sigma_fft(m,r1:r2)/sigma0_fft(0, r1:r2)) + 1d-16 ))
        endelse        
        
        jtot = total(radius(r1:r2)*real_part(sigma_fft(m,r1:r2)*conj(vtheta_fft(m,r1:r2)))*dVol_1D(r1:r2))
        if (m ne 0d0) then jtot *= 2d0 
        mode_ang(j, k-start) = jtot
     endfor

;compute total angular momentum 
     sigma_vtheta = sigma*vphi 
     tot_ang(k-start) = total(sigma_vtheta*radius_2D*dVol_2D)
;total mass
     mass = total(real_part(sigma_fft(0,*))*dVol_1D)

  print, k , tot_ang(k-start)/angmom0 - 1d0, total(mode_ang(*,k-start))/angmom0 - 1d0 , mass, format='(I03,x,3(e22.15,x))' 

;toomre q
  kappa2 = deriv(radius, radius^2d0*real_part(vtheta_fft(0,*))^2d0)/radius^3d0
  cs = smallh/sqrt(radius)      ;locally isothermal only 
  toomreq = cs*sqrt(kappa2)/(!dpi*real_part(sigma_fft(0,*)))
  min_q(k-start) = min(toomreq(r1:r2))
;mass between r1 and r2
  ring_mass(k-start) = total(real_part(sigma_fft(0,r1:r2))*dVol_1D(r1:r2))
  endfor
  

  if not keyword_set(max) then begin
  ytitle=textoidl('log_{10}(C_m/C_{0,t=0})')
  endif else begin
  ytitle=textoidl('max(log_{10}|\Sigma_m/\Sigma_{0,t=0}|)') 
  endelse


  set_plot, 'ps'
  device, filename=filepath('nonaxi_evol.ps',root_dir='.',subdir=[location]) $
          , xsize=8, ysize=4.5, xoffset=0, yoffset=0 $
          , /inches,/color,bits_per_pixel=8
  if not keyword_set(yrange) then begin 
  yrange1=[min(mode_amp),max(mode_amp)]
  endif else yrange1=yrange 

  plot, time_norm(*), mode_amp(0,*), xmargin=[6.5,1.5], ymargin=[3.5,0.5] , ytitle=ytitle, ystyle=1 $
        , charsize=1.5,thick=4,linestyle=0, xtitle=textoidl('t/P_0'), xrange=xrange, yrange=yrange1, xstyle=1
  for j=1, nmodes-1 do begin
  oplot, time_norm(*), mode_amp(j,*), thick=4, linestyle=j
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

  tot = mode_ang(0,*)
  for i=0, finish-start do tot(i) = total(mode_ang(*,i))

  for j=0, nmodes-1 do begin
  m=azimodes(j)
  if (m lt 1) then mode_ang(j,*) = mode_ang(j,*) - angmom0;the reference state only has ang mom in m=0
  endfor
  mode_ang /= angmom0

  if not keyword_set(yrange) then begin 
  yrange1=[min(mode_ang),max(mode_ang)]
  endif else yrange1=yrange

  set_plot, 'ps'
  device, filename=filepath('nonaxi_evol_ang.ps',root_dir='.',subdir=[location]) $
          , xsize=8, ysize=4.5, xoffset=0, yoffset=0 $
          , /inches,/color,bits_per_pixel=8
  plot, time_norm(*), mode_ang(0,*), xmargin=[8.5,1.5], ymargin=[3.5,0.5] , ytitle=textoidl('\DeltaJ_m/J_{ref}'), ystyle=1 $
        , charsize=1.5,thick=4,linestyle=0, xtitle=textoidl('t/P_0'), xrange=xrange, yrange=yrange1, xstyle=1
  for j=1, nmodes-1 do begin
  oplot, time_norm(*), mode_ang(j,*), thick=4, linestyle=j
  endfor
  oplot, time_norm(*), tot/angmom0-1., thick=1, linestyle=0
;  oplot, time_norm(*), tot_ang/angmom0-1., thick=1, linestyle=0
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

  openw,1,filepath('nonaxi_evol_ang.dat',root_dir='.',subdir=[location])
  for i=0, finish-start do begin
  printf, 1, time_norm(i), mode_ang(0:nmodes-1,i), tot(i), format='(4(e22.15,x))'
  endfor
  close,1  

  tot_ang /= angmom0
  tot_ang -= 1d0 
  set_plot, 'ps'
  device, filename=filepath('nonaxi_evol_totj.ps',root_dir='.',subdir=[location]) $
          , xsize=8, ysize=4.5, xoffset=0, yoffset=0 $
          , /inches,/color,bits_per_pixel=8
  plot, time_norm(*), tot_ang(*), xmargin=[8.5,1.5], ymargin=[3.5,0.5] , ytitle=textoidl('\Delta J_{tot}/J_{ref}'), ystyle=1 $
        , charsize=1.5,thick=4,linestyle=0, xtitle=textoidl('t/P_0'), xrange=xrange, yrange=[min(tot_ang),max(tot_ang)], xstyle=1
  device,/close

  set_plot, 'ps'
  device, filename=filepath('nonaxi_evol_minq.ps',root_dir='.',subdir=[location]) $
          , xsize=8, ysize=4.5, xoffset=0, yoffset=0 $
          , /inches,/color,bits_per_pixel=8
  plot, time_norm(*),min_q, xmargin=[8.5,1.5], ymargin=[3.5,0.5] , ytitle=textoidl('min(Q)'), ystyle=1 $
        , charsize=1.5,thick=4,linestyle=0, xtitle=textoidl('t/P_0'), xrange=xrange, yrange=[min(min_q),max(min_q)], xstyle=1
  device,/close


  set_plot, 'ps'
  device, filename=filepath('nonaxi_evol_ringmass.ps',root_dir='.',subdir=[location]) $
          , xsize=8, ysize=4.5, xoffset=0, yoffset=0 $
          , /inches,/color,bits_per_pixel=8
  plot, time_norm(*),ring_mass, xmargin=[8.5,1.5], ymargin=[3.5,0.5],$
        ytitle=textoidl('min(Q)'), ystyle=1 $
        , charsize=1.5,thick=4,linestyle=0, xtitle=textoidl('t/P_0'), xrange=xrange, $
        yrange=[min(ring_mass),max(ring_mass)], xstyle=1
  device,/close
end 
