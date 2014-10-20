function get_vorticity, vrad, vtheta, vphi, rad, theta, phi
  nrad = n_elements(rad)
  ntheta = n_elements(theta)
  nphi = n_elements(phi)
  
  dvphi_sintheta = vphi
  dvtheta = vtheta
  
;get omega_r*cos(theta)
  for i=0, nrad-1 do begin
     for k=0, nphi-1 do begin
        dvphi_sintheta(i,*,k) = deriv(theta(*), vphi(i,*,k)*sin(theta(*)))
        dvphi_sintheta(i,*,k) /= rad(i)*tan(theta(*))
     endfor
     
     for j=0, ntheta-1 do begin
        dvtheta(i,j,*) = deriv(phi(*), vtheta(i,j,*))
        dvtheta(i,j,*) /= rad(i)*tan(theta(j))
     endfor
  endfor
  omegar = dvphi_sintheta - dvtheta
  
;get omega_theta*sin(theta) !factors involving sin(theta) in the expression for the curl is absorbed with pre-factor sin(theta) here
  dvrad = vrad
  dr_vphi = vphi
  
  for j=0, ntheta-1 do begin
     for i=0, nrad-1 do begin
        dvrad(i,j,*) = deriv(phi(*), vrad(i,j,*))/rad(i)
     endfor
     
     for k=0, nphi-1 do begin
        dr_vphi(*,j,k) = deriv(rad(*), rad(*)*vphi(*,j,k))
        dr_vphi(*,j,k) *= sin(theta(j))/rad(*)
     endfor
  endfor
  omegatheta = dvrad - dr_vphi
  
  omegaz = omegar - omegatheta  ; geometric factors already incoroporated. 
  
  return, omegaz
end

pro pdisk_vort_rz2, start=start, finish=finish,xrange=xrange, yrange=yrange, $
                   ct=ct, plotrange=plotrange, r0=r0, phi0=phi0, $
                   ytickinterval=ytickinterval, smallh=smallh, mp=mp, xtickinterval=xtickinterval, nonaxi=nonaxi, $
                   basic=basic
  
  COMMON PLUTO_GRID,  nx1,nx2,nx3,x1,x2,x3,$
     dx1,dx2,dx3, geometry, xpos, ypos,$
     AMRLevel, AMRBoxes         ;  ** Chombo data structure **
                                ;  ** loaded when HDF5LOAD is called **
  
  COMMON PLUTO_VAR, NVAR, rho, vx1, vx2, vx3, $
     bx1, bx2, bx3, $
     Ax1, Ax2, Ax3, $
     bx1s, bx2s, bx3s,$
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
  
  if not keyword_set(ct) then ct =5   
  if not keyword_set(finish) then finish=start 
  if not keyword_set(zslice) then zslice= 0.99
  if not keyword_set(r0) then r0=1.0
  if not keyword_set(basic) then basic = 0
  
  pload, basic
  d0 = rho
  
  if keyword_set(basic) then begin
     for j=0, nx2-1 do begin
        for i=0, nx1-1 do begin
           d0(i,j,*) = mean(d0(i,j,*))
        endfor
     endfor
  endif
  
  data     = dblarr(nx1, nx2)
  dataRZ   = dblarr(nx1, nx2)
  daxi     = dblarr(nx1, nx2, nx3)
  data_xy  = dblarr(nx1, nx3)
  vortz_2d = dblarr(nx1, nx2)
  

  
;assume uniform spacing in r and theta
  dr  = x1(1) - x1(0)
  dth = x2(1) - x2(0)
  
  zmin = 0.0
  zmax = x1(nx1-1)*tan(!dpi/2.0 - x2(0))
  dz   = (zmax - zmin)/(nx2-1.0)
  zaxis= zmin + dindgen(nx2)*dz
  
  temp = min(abs(x1 - xrange(0)*r0), r1)
  temp = min(abs(x1 - xrange(1)*r0), r2)
  Rmin = x1(r1)
  Rmax = x1(r2)
  dbR   = (Rmax - Rmin)/(nx1 - 1.0)
  Raxis= Rmin + dindgen(nx1)*dbR
  
  if keyword_set(mp) then begin
     rh    = (mp/3.0)^(1.0/3.0)*r0
     rad  = (x1-r0)/rh
     xtitle=textoidl('(R-r_0)/r_h')
  endif else begin
     rad   = x1/r0
     xtitle=textoidl('R/r_0')
  endelse
  
  ytitle = textoidl('z/H(r_0)')
  
  azi   = x3/(2d0*!dpi)
  time  = t/(2d0*!dpi*r0^1.5)
   
  loadct,ct,/silent
  for i=start, finish do begin
     pload, i,/silent
     
     d    = rho
     vortz = get_vorticity(vx1, vx2, vx3, x1, x2, x3)
     
     for jj=0, nx2-1 do begin
        for ii=0, nx1-1 do begin
           avg = mean(d(ii,jj,*))
           daxi(ii,jj,*) = avg
           
           avg = mean(vx3(ii,jj,*))
           vx3(ii,jj,*) = avg/(x1(ii)*sin(x2(jj)))
           
           vortz_2d(ii,jj) = mean(vortz(ii,jj,*))
        endfor
     endfor
     
     if keyword_set(nonaxi) then begin
        for k=0, nx3-1 do vortz(*,*,k) -= vortz_2d(*,*)
        d0 = daxi
     endif else begin
        d0 = dzero
     endelse

     vortz /= 2d0*vx3
     
     if (phi0 lt 0.0) then begin
        d /= d0
        data_xy(*,*) = d(*,nx2-1,*)
        temp = max(data_xy, maxloc)
        ind  = array_indices(data_xy, maxloc)
        nphi = ind(1)
        azi_string = strcompress(string(azi(nphi), format='(f4.2)'),/remove_all)
     endif else begin
        nphi = phi0*nx3
        azi_string = strcompress(string(azi(nphi), format='(f4.2)'),/remove_all)
     endelse
     
     data(*,*) = vortz(*,*,nphi)
     
     for jj=0, nx2-1 do begin
        z = zaxis(jj)
        for ii=0, nx1-1 do begin
           R = Raxis(ii)
           
           r_t = sqrt(R*R + z*z)
           th_t = atan(R, z)
           
           if(th_t ge x2(0)) then begin
              
              temp = min(abs(x1 - r_t),   x0)
              temp = min(abs(x2 - th_t),y0)
              ip = x0 + (r_t - x1(x0))/(dr/2.0) + 0.0
              jp = y0 + (th_t - x2(y0))/(dth/2.0) + 0.0
              dataRZ(ii,jj) = bilinear(data, ip, jp)

              ;; v=data(*,jj)
              ;; u=[r_t]
              ;; result1 = interpol(v, x1, u,/quadratic)

              ;; v=data(ii,*)
              ;; u=[th_t]
              ;; result2 = interpol(v, x2, u,/quadratic)
              
              ;; dataRZ(ii,jj) = mean([result1,result2])

           endif else dataRZ(ii,jj) = 10.0
           
        endfor
     endfor
     
     if not keyword_set(plotrange) then begin
        mask = where(dataRZ ne 10.0)
        plotrange0 = [min(dataRZ(mask)),max(dataRZ(mask))]
     endif else begin
        plotrange0 = plotrange
     endelse
     levels = plotrange0(0) + (plotrange0(1)-plotrange0(0))*dindgen(16)/15d0
     
     
     name2 = string(i, format='(I03)')
     title1 = textoidl('t=')+strcompress(string(time(i), format='(f5.1)'),/remove_all)+textoidl('P_0')
     title  = title1+ textoidl(', \phi/2\pi=')+azi_string
     set_plot,'ps'
     device, filename=strcompress('pdisk_vort_rz2_'+name2+'.ps',/remove_all) $
             , bits_per_pixel=8, xsize=18, ysize=9,/color
     
     contour,dataRZ,Raxis/r0,zaxis/(r0*smallh),levels=levels, xstyle=1,/fill $
             , xtitle=xtitle, ytitle=ytitle,charsize=1.0,xmargin=[7.0,10.0], ymargin=[4,2], $
             xtickinterval=xtickinterval,xrange=xrange, yrange=yrange, title=title ;, xminor=5
     colorbar, position=[0.88, 0.16, 0.93, 0.92],/vertical,/right,range=plotrange0,format='(f5.2)'
     device,/close
  endfor
end
