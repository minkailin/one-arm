pro pdisk_Q1d, start=start, finish=finish, xrange=xrange, yrange=yrange, $
               r0=r0, legend=legend, label=label, smth=smth $
               ,ytickinterval=ytickinterval, smallh=smallh, smallq=smallq, mp=mp, xtickinterval=xtickinterval 
  
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
  
  if not keyword_set(finish) then finish=start 
  if not keyword_set(r0) then r0=1.0
 
  ytitle = textoidl('Q')
  
  pload, 0 
  
  nrad  = nx1
  ntheta= nx2
  nphi  = nx3
  
  rad   = x1
  theta = x2
  phi   = x3
  
;assume uniform spacing in r and theta
  dr  = rad(1) - rad(0)
  dth = theta(1) - theta(0)
  
  time  = t/(2d0*!dpi*r0^1.5)
  
  den_2d = dblarr(nrad, ntheta)
  vphi_2d = dblarr(nrad, ntheta)
  
  zmin = 0.0
  zmax = rad(nx1-1)*tan(!dpi/2.0 - theta(0))
  dz   = (zmax - zmin)/(ntheta-1.0)
  zaxis= zmin + dindgen(ntheta)*dz
  
  temp = min(abs(rad - xrange(0)*r0), r1)
  temp = min(abs(rad - xrange(1)*r0), r2)
  Rmin = rad(r1)
  Rmax = rad(r2)
  nR   = r2 - r1 + 1
  dbR   = (Rmax - Rmin)/(nR - 1.0)
  Raxis= Rmin + dindgen(nR)*dbR
  temp = min(abs(Raxis - r0), rgrid)
  
  
  denRZ = dblarr(nR, ntheta)
  vphiRZ = dblarr(nR, ntheta)

  sigm1d = dblarr(nR)
  vphi1d = dblarr(nR)
  
  for i=start, finish do begin
     pload, i,/silent
     
     for jj=0, ntheta-1 do begin
        for ii=0, nrad-1 do begin
           den_2d(ii,jj) = mean(rho(ii,jj,*)) 
           vphi_2d(ii,jj)= mean(rho(ii,jj,*)*vx3(ii,jj,*));momentum density
        endfor
     endfor    
     
     for jj=0, ntheta-1 do begin
        z = zaxis(jj)
        for ii=0, nR-1 do begin
           R = Raxis(ii)
  
           r_t = sqrt(R*R + z*z)
           th_t = atan(R, z)

           if(th_t ge theta(0)) then begin
              
              temp = min(abs(rad - r_t),   x0)
              temp = min(abs(theta - th_t),y0)
              
              ip = x0 + (r_t - rad(x0))/(dr/2.0) + 0.0
              jp = y0 + (th_t - theta(y0))/(dth/2.0) + 0.0 
              
              denRZ(ii,jj) = bilinear(den_2d, ip, jp)
              vphiRZ(ii,jj)= bilinear(vphi_2d, ip, jp)
           endif else begin
              denRZ(ii,jj) = 0.0
              vphiRZ(ii,jj) = 0.0
           endelse
        endfor
     endfor  
        
     for ii=0, nR-1 do begin
        sigm1d(ii) = int_tabulated(zaxis, denRZ(ii,*))
        vphi1d(ii) = int_tabulated(zaxis, vphiRZ(ii,*))/sigm1d(ii)  
     endfor
     
     if keyword_set(smth) then begin
        sigm1d = smooth(sigm1d, smth)
        vphi1d = smooth(vphi1d, smth)
     endif

     kappa = deriv(Raxis, Raxis^2*vphi1d^2)/Raxis^3
     kappa = sqrt(kappa)
     cs = smallh*r0^(-0.5)*(r0/Raxis)^(smallq/2.0)
     
     toomreQ = kappa*cs
     toomreQ/= !dpi*2.0*sigm1d  ;factor of two to accoutn for lower half disk
     
     name2 = string(i, format='(I03)')
     title = textoidl('t=')+strcompress(string(time(i), format='(f5.1)'),/remove_all)+textoidl('P_0')
     
     set_plot, 'ps'
     device, filename=strcompress('pdisk_Q1d_'+name2+'.ps',/remove_all) $
             ,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches
     plot, Raxis, toomreQ,xmargin=[8.5,1.5],ymargin=[3.5,1.5], ystyle=0, xstyle=1 $
           ,charsize=1.5, thick=4, xrange=xrange, title=title, xtitle=textoidl('R/r_0'), yrange=yrange,$
           linestyle = 0, ytitle =ytitle, xtickinterval=xtickinterval, ytickinterval=ytickinterval
     device,/close 
     
  endfor
end
