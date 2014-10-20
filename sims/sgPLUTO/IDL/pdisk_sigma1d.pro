pro pdisk_sigma1d, loc=loc, start=start, finish=finish, xrange=xrange, yrange=yrange, $
                   r0=r0, legend=legend, label=label, nz=nz, nrad=nrad, basic=basic $
                   ,ytickinterval=ytickinterval, smallh=smallh, mp=mp, xtickinterval=xtickinterval,nopert=nopert
  
  COMMON PLUTO_GRID,  nx1,nx2,nx3,x1,x2,x3,$
     dx1,dx2,dx3, geometry, xpos, ypos,$
     AMRLevel, AMRBoxes         ;  ** Chombo data structure **
                                ;  ** loaded when HDF5LOAD is called **
  
  COMMON PLUTO_VAR, NVAR, rho, vx1, vx2, vx3, $
     bx1, bx2, bx3, $
     Ax1, Ax2, Ax3, $
     bx1s, bx2s, bx3s, pot, $
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
  
  
  if not keyword_set(finish) then finish=start 
  if not keyword_set(r0) then r0=1.0
  if not keyword_set(nopert) then begin
     ytitle = textoidl('<\Sigma>_\phi/\Sigma(t=0) - 1')
  endif else begin
     ytitle = textoidl('<\Sigma>_\phi/<\Sigma>_\phi(r_0)')
  endelse
  
  if not keyword_set(basic) then begin
  pload, 0 
  endif else pload, basic  

  nrad  = nx1
  ntheta= nx2
  nphi  = nx3

  if not keyword_set(nz) then nz   = ntheta
  if not keyword_set(nrad) then nrad = nx1

  rad   = x1
  theta = x2
  phi   = x3

  theta1 = !dpi/2d0 - theta
  tol = 1d-6*max([max(rad),max(theta1)])

  if not keyword_set(xrange) then xrange=[min(rad), max(rad)]
 
;theta is uniform spacing
  dth = theta(1) - theta(0)
  
  time  = t/(2d0*!dpi*r0^1.5)
  
  den_2d = dblarr(nrad, ntheta)
  
  zmin = 0.0
  zmax = rad(nx1-1)*tan(!dpi/2.0 - theta(0))
  dz   = (zmax - zmin)/(nz-1.0)
  zaxis= zmin + dindgen(nz)*dz
  
  temp = min(abs(rad - xrange(0)*r0), r1)
  temp = min(abs(rad - xrange(1)*r0), r2)
  Rmin = rad(r1)
  Rmax = min([rad(r2), x1(nx1-1)*sin(theta(0))])
  dbR   = (Rmax - Rmin)/(nrad - 1.0)
  Raxis= Rmin + dindgen(nrad)*dbR
  temp = min(abs(Raxis - r0), rgrid)
  
  data2d = dblarr(nrad, nz)
  data1d = dblarr(nrad)
  
  for jj=0, ntheta-1 do begin
     for ii=0, nrad-1 do begin 
        den_2d(ii,jj) = mean(rho(ii,jj,*))      
     endfor
  endfor

  for jj=0, nz-1 do begin
     z = zaxis(jj)
     for ii=0, nrad-1 do begin
        R = Raxis(ii)       
        r_t = sqrt(R*R + z*z)
        th_t = atan(R, z)        
        if(th_t ge theta(0)) then begin
           temp = min(abs(rad - r_t),   x0)
           temp = min(abs(theta - th_t),y0)

           if x0 eq nx1-1 then begin
              dr   = rad(x0) - rad(x0-1)
           endif else if x0 eq 0 then begin
              dr   = rad(x0+1) - rad(x0)
           endif else begin
           dr   = abs(rad(x0-1) - rad(x0+1))/2d0
           endelse 

           ip = x0 + (r_t - rad(x0))/(dr/2.0) + 0.0
           jp = y0 + (th_t - theta(y0))/(dth/2.0) + 0.0
           data2d(ii,jj) = bilinear(den_2d, ip, jp)
        endif else data2d(ii,jj) = 0.0
     endfor
  endfor


;  temp = polar_surface1(den_2d, rad, theta1, /grid, spacing=[dbR, dz], bounds=[Rmin,zmin,Rmax,zmax], tol=tol,/quintic)
;  data2d = congrid(temp, nrad, nz,/interp)

  result = data2d   
  for ii=0, nrad-1 do begin
     data1d(ii) = int_tabulated(zaxis, result(ii,*))
  endfor
  
  data1d_0 = data1d
  
  
  for i=start, finish do begin
     pload, i,/silent
     
     for jj=0, ntheta-1 do begin
        for ii=0, nrad-1 do begin
           den_2d(ii,jj) = mean(rho(ii,jj,*)) 
        endfor
     endfor    
     

     for jj=0, nz-1 do begin
        z = zaxis(jj)
        for ii=0, nrad-1 do begin
           R = Raxis(ii)          
           r_t = sqrt(R*R + z*z)
           th_t = atan(R, z)           
           if(th_t ge theta(0)) then begin             
              temp = min(abs(rad - r_t),   x0)
              temp = min(abs(theta - th_t),y0)
  
              if x0 eq nx1-1 then begin
              dr   = rad(x0) - rad(x0-1)
              endif else if x0 eq 0 then begin
              dr   = rad(x0+1) - rad(x0)
              endif else begin
              dr   = abs(rad(x0-1) - rad(x0+1))/2d0
              endelse

              ip = x0 + (r_t - rad(x0))/(dr/2.0) + 0.0
              jp = y0 + (th_t - theta(y0))/(dth/2.0) + 0.0 
              data2d(ii,jj) = bilinear(den_2d, ip, jp)
           endif else data2d(ii,jj) = 0.0           
        endfor
     endfor  
     
;     temp = polar_surface1(den_2d, rad, theta1, /grid, spacing=[dbR, dz], bounds=[Rmin,zmin,Rmax,zmax],tol=tol,/quintic)
;     data2d = congrid(temp, nrad, nz,/interp)

     result = data2d     
     for ii=0, nrad-1 do begin
        data1d(ii) = int_tabulated(zaxis, result(ii,*))
     endfor
     
     if not keyword_set(nopert) then begin
        data1d /= data1d_0
        data1d -= 1.0
        print, 'max dSigma', max(abs(data1d))
     endif else begin
        data1d /= data1d(rgrid)
     endelse
     
     name2 = string(i, format='(I03)')
     title = textoidl('t=')+strcompress(string(time(i), format='(f5.1)'),/remove_all)+textoidl('P_0')
     
     set_plot, 'ps'
     device, filename=strcompress('pdisk_sigma1d_'+name2+'.ps',/remove_all) $
             ,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches
     plot, Raxis, data1d,xmargin=[8.5,1.5],ymargin=[3.2,1.8], ystyle=0, xstyle=1 $
           ,charsize=1.5, thick=4, xrange=xrange, title=title, xtitle=textoidl('R/r_0'), yrange=yrange,$
           linestyle = 0, ytitle =ytitle, xtickinterval=xtickinterval, ytickinterval=ytickinterval,charthick=2
     device,/close 
     
     file =strcompress('pdisk_sigma1d_'+name2+'.dat',/remove_all)
     openw,1, file
     for ii=0, nrad-1 do begin
        printf,1,Raxis(ii),data1d(ii),format='(2(e22.15,x))'
     endfor
     close,1
     
  endfor
end
