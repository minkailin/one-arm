pro compare_angmom,  rrange=rrange, smallh=smallh, r0=r0, loc=loc, basic=basic, azimodes=azimodes, start=start, finish=finish, xrange=xrange, yrange=yrange, legend=legend, label=label, scale=scale

!p.font=0

if not keyword_set(smallh) then smallh=0.05 
if not keyword_set(finish) then finish=start
if not keyword_set(basic) then basic = 0

nmodes = n_elements(azimodes) 

if not keyword_set(r0) then r0=1d0

if not keyword_set(scale) then begin
scale = dblarr(nmodes)
scale(*) = 1d0
endif 

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
  
  if not keyword_set(rrange) then begin
  r1=0
  r2=nrad-1
  endif
  temp = min(abs(radius-rrange(0)),r1)
  temp = min(abs(radius-rrange(1)),r2)

  phi = 2d0*!dpi*dindgen(nsec)/(nsec - 1d0) 
;;;;;;;;;;
;GET DATA;
;;;;;;;;;;
  data   = dblarr(nsec,nrad)
  data1d = dblarr(nmodes,nrad) 

  time_norm = dblarr(finish-start+1)
  dtotj0_sigdvphi = time_norm
  dtotj0_dsigvphi = time_norm
  dtotj0_dsigdvphi= time_norm

  angmom_source_dens=dblarr(nrad)
  angmom_source     = dblarr(finish-start+1)

  angmom1_source     = dblarr(finish-start+1)
  angmom1_source0     = dblarr(finish-start+1)

;use axisymmetric component of t=0 as normalization
  openr,2,filepath(strcompress('gasdens'+string(basic)+'.dat',/remove_all),root_dir='.',subdir=[location]) 
  readu,2,data
  close,2
  sigma0_fft = fft(data, -1, dimension=1,/double)
  sig0 = real_part(sigma0_fft(0,*))

  openr,2,filepath(strcompress('gasvtheta'+string(basic)+'.dat',/remove_all),root_dir='.',subdir=[location])
  readu,2,data
  close,2
  vtheta0_fft = fft(data, -1, dimension=1,/double)
  kappa2 = deriv(radius, (radius*real_part(vtheta0_fft(0,*)))^2)/radius^3
  cs     = smallh/sqrt(radius)

  vphi0 = real_part(vtheta0_fft(0,*))

  x1=0
  x2=nrad-1
  j0_init = 2d0*!dpi*radius(*)*real_part(vtheta0_fft(0,*))*real_part(sigma0_fft(0,*))
  angmom0 = int_tabulated(radius(x1:x2), 2d0*!dpi*radius(x1:x2)*radius(x1:x2)*real_part(vtheta0_fft(0,x1:x2))*real_part(sigma0_fft(0,x1:x2)))
  

  for k=start, finish do begin
    time_norm(k-start) = info(7,k)/p0
    time=string(time_norm(k-start),format='(F7.2)')
    title = 't='+time+textoidl('P_0, ')

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
     vtheta=data 

     openr,2,filepath(strcompress('gasvrad'+string(k)+'.dat',/remove_all),root_dir='.',subdir=[location])
     readu,2,data
     close,2
     vrad=data
     for j = 0, nrad-2 do vrad(*,j) = (data(*,j) + data(*,j+1))/2.0
     vrad(*,nrad-1) = data(*,nrad-1)

     
;fft      
     sigma_fft    = fft(sigma, -1, dimension=1,/double)
     vtheta_fft   = fft(vtheta, -1, dimension=1,/double)
     vrad_fft     = fft(vrad, -1, dimension=1,/double)     

     dsig0 = real_part(sigma_fft(0,*)) - sig0
     dvphi0= real_part(vtheta_fft(0,*))- vphi0


     dj0_sigdvphi = 2d0*!dpi*radius*sig0*dvphi0
     dj0_dsigvphi = 2d0*!dpi*radius*dsig0*vphi0
     dj0_dsigdvphi= 2d0*!dpi*radius*dsig0*dvphi0

     dtotj0_sigdvphi(k-start) = int_tabulated(radius, dj0_sigdvphi*radius)/angmom0
     dtotj0_dsigvphi(k-start) = int_tabulated(radius, dj0_dsigvphi*radius)/angmom0
     dtotj0_dsigdvphi(k-start )=int_tabulated(radius, dj0_dsigdvphi*radius)/angmom0


     for j=0, nmodes-1 do begin
        m = azimodes(j)

     jm = 4d0*!dpi*radius*real_part(sigma_fft(m,*)*conj(vtheta_fft(m,*)))
     if(m lt 1d0) then begin
        jm /= 2d0 
        jm -= j0_init
     endif
 
     data1d(j,*) = jm/angmom0  
    endfor
 
  temp = sigma_fft(0,*)*vrad_fft(0,*) + sigma_fft(1,*)*conj(vrad_fft(1,*)) + conj(sigma_fft(1,*))*vrad_fft(1,*) ;+ sigma_fft(2,*)*conj(vrad_fft(2,*)) + conj(sigma_fft(2,*))*vrad_fft(2,*)
  angmom_source_dens = real_part(temp)
;  angmom_source(k-start)=int_tabulated(radius(r1:r2), radius(r1:r2)^(0.5)*deriv(radius(r1:r2),radius(r1:r2)*angmom_source_dens(r1:r2))) 

  angmom_source(k-start) =  radius(r2)^1.5*angmom_source_dens(r2) - radius(r1)^1.5*angmom_source_dens(r1)
  angmom_source(k-start)+= -0.5*int_tabulated(radius(r1:r2), sqrt(radius(r1:r2))*angmom_source_dens(r1:r2))


  temp = sigma_fft(0,*)*vtheta_fft(1,*)+sigma_fft(1,*)*vtheta_fft(0,*)
  angmom1_source(k-start) =1d5*abs( dcomplex(int_tabulated(radius,real_part(temp)),$
                                int_tabulated(radius,imaginary(temp))))

  temp = sigma0_fft(0,*)*vtheta_fft(1,*)+sigma_fft(1,*)*vtheta0_fft(0,*)
  angmom1_source0(k-start) =1d5*abs( dcomplex(int_tabulated(radius,real_part(temp)),$
                                int_tabulated(radius,imaginary(temp))))

  if not keyword_set(yrange) then begin
  yrange0=[min(data1d),max(data1d)]
  endif else yrange0=yrange 

  set_plot, 'ps'
  device, filename=filepath('compare_angmom_'+ks+'.ps',root_dir='.',subdir=[location]) $
          , xsize=8, ysize=4.5, xoffset=0, yoffset=0 $
          , /inches,/color,bits_per_pixel=8
  plot, radius/r0, data1d(0,*)*scale(0), xmargin=[6.5,1.5], ymargin=[3.5,1.7] , ytitle=textoidl('\Deltaj_m/J_{ref}'), ystyle=1 $
        , charsize=1.5,thick=4,linestyle=0, xtitle=textoidl('r/r_0'), xrange=xrange, yrange=yrange0, xstyle=1, title=title
  for j=1, nmodes-1 do begin
  oplot, radius/r0, data1d(j,*)*scale(j), thick=4, linestyle=j
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

  set_plot, 'ps'
  device, filename=filepath('compare_angmom0_'+ks+'.ps',root_dir='.',subdir=[location]) $
          , xsize=8, ysize=4.5, xoffset=0, yoffset=0 $
          , /inches,/color,bits_per_pixel=8
  plot, radius/r0, 2d0*!dpi*radius*sig0*dvphi0/angmom0, xmargin=[6.5,1.5], ymargin=[3.5,1.7] , ytitle=textoidl('\Deltaj_0/J_{ref}'), ystyle=1 $
        , charsize=1.5,thick=4,linestyle=0, xtitle=textoidl('r/r_0'), xrange=xrange, yrange=yrange0, xstyle=1, title=title
  oplot, radius/r0, 2d0*!dpi*radius*dsig0*vphi0/angmom0, thick=4, linestyle=1
  oplot, radius/r0, 2d0*!dpi*radius*dsig0*dvphi0/angmom0, thick=4, linestyle=2

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


;  dsig0 = real_part(sigma_fft(0,*))-real_part(sigma0_fft(0,*))
;  kr    = deriv(radius, dsig0)/(dsig0 + 1d-16) 
;  om2   = kappa2 + cs^2*abs(kr)^2 - 2.0*!dpi*abs(kr)*real_part(sigma0_fft(0,*))
;  om2  /= (real_part(vtheta0_fft(0,*))/radius)^2  
;
;  set_plot, 'ps'
;  device, filename=filepath('compare_angmom_om2_'+ks+'.ps',root_dir='.',subdir=[location]) $
;          , xsize=8, ysize=4.5, xoffset=0, yoffset=0 $
;          , /inches,/color,bits_per_pixel=8
;  plot, radius/r0, om2, xmargin=[6.5,1.5], ymargin=[3.5,1.7] , ytitle=textoidl('\omega_0^2/\Omega^2'), ystyle=1 $
;        , charsize=1.5,thick=4,linestyle=0, xtitle=textoidl('r/r_0'), xrange=xrange, yrange=yrange, xstyle=1, title=title
;
;   if keyword_set(legend) then begin
;   x0=legend(0)
;   x1=legend(1)
;   y0=legend(2)
;   dy=legend(3)
;        for j=0, n_elements(label)-1 do begin
;    oplot, [x0,x1], [y0,y0]-dy*j, thick=4, linestyle=j
;    xyouts, x1, y0-dy*j,textoidl(label(j)),charsize=1.5
;    endfor
;    endif
;  device,/close

  print, ks 
  endfor

set_plot, 'ps'
  device, filename=filepath('compare_angmom0.ps',root_dir='.',subdir=[location]) $
          , xsize=8, ysize=4.5, xoffset=0, yoffset=0 $
          , /inches,/color,bits_per_pixel=8
  plot, time_norm,dtotj0_sigdvphi ,xmargin=[6.5,1.5], ymargin=[3.5,0.5] , ytitle=textoidl('\DeltaJ_0/J_{ref}'), ystyle=1 $
        , charsize=1.5,thick=4,linestyle=0, xtitle=textoidl('t/P_0'), xrange=xrange, $
        yrange=[min([min(dtotj0_sigdvphi),min(dtotj0_dsigvphi),min(dtotj0_dsigdvphi)]),max([max(dtotj0_sigdvphi),max(dtotj0_dsigvphi),max(dtotj0_dsigdvphi)])], xstyle=1
  oplot, time_norm, dtotj0_dsigvphi, thick=4, linestyle=1
  oplot, time_norm, dtotj0_dsigdvphi, thick=4, linestyle=2

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
  
  set_plot, 'ps'
  device, filename=filepath('compare_angmom0_source.ps',root_dir='.',subdir=[location]) $
          , xsize=8, ysize=4.5, xoffset=0, yoffset=0 $
          , /inches,/color,bits_per_pixel=8
  plot, time_norm, -angmom_source*1d8,xmargin=[6.5,1.5], ymargin=[3.5,0.5] , ytitle=textoidl('J_0 source'), ystyle=1 $
        , charsize=1.5,thick=4,linestyle=0, xtitle=textoidl('t/P_0'), xrange=xrange, xstyle=1
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

  set_plot, 'ps'
  device, filename=filepath('compare_angmom1_source.ps',root_dir='.',subdir=[location]) $
          , xsize=8, ysize=4.5, xoffset=0, yoffset=0 $
          , /inches,/color,bits_per_pixel=8
  plot, time_norm, angmom1_source,xmargin=[6.5,1.5], ymargin=[3.5,0.5] , ytitle=textoidl('m_1 source'), ystyle=1 $
        , charsize=1.5,thick=4,linestyle=0, xtitle=textoidl('t/P_0'), xrange=xrange, xstyle=1
  oplot, time_norm, angmom1_source0, thick=4, linestyle=1
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
