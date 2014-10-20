pro pdisk_sigma, loc=loc, start=start, finish=finish, log=log, nopert=nopert, xrange=xrange, yrange=yrange, $
                 float=float, ct=ct, zslice=zslice, plotrange=plotrange, r0=r0,phi0=phi0, $
                 ytickinterval=ytickinterval, smallh=smallh, mp=mp, xtickinterval=xtickinterval, nonaxi=nonaxi, basic=basic, cart=cart, $
                 nz=nz

COMMON PLUTO_GRID,  nx1,nx2,nx3,x1,x2,x3,$
                    dx1,dx2,dx3, geometry, xpos, ypos,$
                    AMRLevel, AMRBoxes;  ** Chombo data structure **
                                      ;  ** loaded when HDF5LOAD is called **

COMMON PLUTO_VAR, NVAR, rho, vx1, vx2, vx3, $
                             bx1, bx2, bx3, $
                             Ax1, Ax2, Ax3, $
                             bx1s, bx2s, bx3s, pot, $
                         ; ----------------------------------------
                             v1, v2, v3, $   ; Kept for backward
                             b1, b2, b3, $   ; compatibility with 
                             A1, A2, A3, $   ; PLUTO 3 data name
                             b1s, b2s, b3s, $ ;
                             pr,            $ ;
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

if not keyword_set(ct) then ct =5   
if not keyword_set(finish) then finish=start 
if not keyword_set(zslice) then zslice= 0.99
if not keyword_set(r0) then r0=1.0
if not keyword_set(basic) then basic = 0
if not keyword_set(nz) then nz   = ntheta

if not keyword_set(float) then begin
   pload, basic
endif else pload, basic, /float
d0 = rho

nrad  = nx1
ntheta= nx2
nphi  = nx3

rad   = x1
theta = x2
phi   = x3

theta1 = !dpi/2d0 - theta
dth    = theta(1) - theta(0)
tol = 1d-6*max([max(rad),max(theta1)])

if not keyword_set(xrange) then xrange=[min(x1),max(x1)]

if keyword_set(basic) then begin
   for j=0, nx2-1 do begin
      for i=0, nx1-1 do begin
         d0(i,j,*) = mean(d0(i,j,*))
      endfor
   endfor
endif

data     = dblarr(nx1, nx3)
array    = dblarr(nx1, nx3)
dataplot = dblarr(nx1, nx3)

;normalize axis, time

azi   = x3
time  = t/(2d0*!dpi*r0^1.5)

;get surface density of the reference density (axisymmetric)
;theta is uniform spacing

nrad = nx1

den_2d = dblarr(nrad, ntheta)
den_rz = dblarr(nrad, nz)

surf0  = dblarr(nrad, nphi)
surf   = dblarr(nrad, nphi)

zmin = 0.0
zmax = rad(nx1-1)*tan(!dpi/2.0 - theta(0))
dz   = (zmax - zmin)/(nz-1.0)
zaxis= zmin + dindgen(nz)*dz

Rmin = rad(0)
Rmax = x1(nx1-1)*sin(theta(0))
dbR   = (Rmax - Rmin)/(nrad - 1.0)
Raxis= Rmin + dindgen(nrad)*dbR
temp = min(abs(Raxis - r0), rgrid)

if keyword_set(mp) then begin
   rh     = (mp/3.0)^(1.0/3.0)*r0
   rplot  = (Raxis-r0)/rh
   xtitle=textoidl('(R-r_p)/r_h')
endif else begin
   rplot   = Raxis/r0
   xtitle  = textoidl('R/r_0')
endelse

den_2d(*,*) = d0(*,*,0) ;this den field is axsymm   
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
         
         ip = x0 + (r_t - rad(x0))/(dr/2d0) + 0.0
         jp = y0 + (th_t - theta(y0))/(dth/2d0) + 0.0
         den_rz(ii,jj) = bilinear(den_2d, ip, jp)
      endif else den_rz(ii,jj) = 0.0
   endfor
endfor
for ii=0, nrad-1 do begin
   surf0(ii,*) = int_tabulated(zaxis, den_rz(ii,*))
endfor

if keyword_set(cart) then begin
   dazi = dx3(0)
   dnx1 = 2*nrad
   dx   = 2*Raxis(nrad-1)/(dnx1-1.0)
   xaxis = -rad(nrad-1) + dx*dindgen(dnx1)
   yaxis = xaxis
   dataxy = dblarr(dnx1, dnx1)
endif

if keyword_set(phi0) then begin
   ytitle = textoidl('(\phi-\phi_0)/\pi')
   azi1 = x3/!dpi - 1.0
   if not keyword_set(yrange) then yrange=[-1,1]
endif else begin
   ytitle = textoidl('\phi/2\pi')
   azi1 = x3/(2d0*!dpi)
   if not keyword_set(yrange) then yrange=[0,1]
endelse

loadct,ct,/silent

for i=start, finish do begin
   if not keyword_set(float) then begin   
      pload, i,/silent
   endif else pload, i, /silent,/float
   d = rho
   
   if keyword_set(nonaxi) then begin
      for jj=0, nx2-1 do begin
         for ii=0, nx1-1 do begin
            avg = mean(d(ii,jj,*))
            d0(ii,jj,*) = avg
         endfor
      endfor
      den_2d(*,*) = d0(*,*,0)
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
               
               ip = x0 + (r_t - rad(x0))/(dr/2d0) + 0.0
               jp = y0 + (th_t - theta(y0))/(dth/2d0) + 0.0
               den_rz(ii,jj) = bilinear(den_2d, ip, jp)
            endif else den_rz(ii,jj) = 0.0
         endfor
      endfor
      for ii=0, nrad-1 do begin
         surf0(ii,*) = int_tabulated(zaxis, den_rz(ii,*))
      endfor
   endif
   
   for kk=0, nphi-1 do begin
      den_2d(*,*) = d(*,*,kk)
      
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
               
               ip = x0 + (r_t - rad(x0))/(dr/2d0) + 0.0
               jp = y0 + (th_t - theta(y0))/(dth/2d0) + 0.0
               den_rz(ii,jj) = bilinear(den_2d, ip, jp)
            endif else den_rz(ii,jj) = 0.0
         endfor
      endfor
      for ii=0, nrad-1 do begin
         surf(ii,kk) = int_tabulated(zaxis, den_rz(ii,*))
      endfor
   endfor
   
   if not keyword_set(nopert) then begin
      surf /= surf0
      if not keyword_set(log) then   surf -= 1.0
   endif
   if keyword_set(log) then surf=alog10(surf)
   data(*,*) = surf
   
   if keyword_set(phi0) then begin
      
      if(phi0 gt 0.0) then begin
         phiplt = phi0*2d0*!dpi
      endif else begin
         array(*,*) = d(*,nx2-1,*)
         temp = max(array, maxloc)
         ind  = array_indices(array, maxloc)
         phiplt = azi(ind(1)) + !dpi
         print, 'density center is', phiplt/(2.0*!dpi)
      endelse
      
      if(phiplt gt !dpi) then begin
         temp2=min(abs(azi-(phiplt-!dpi)),grid2)
         dataplot(0:nx1-1,0:nx3-1-grid2)   = data(0:nx1-1,grid2:nx3-1)
         dataplot(0:nx1-1,nx3-grid2:nx3-1) = data(0:nx1-1, 0:grid2-1)
      endif
      if(phiplt lt !dpi) then begin
         temp2=min(abs(azi-(phiplt+!dpi)),grid2)
         dataplot(0:nx1-1,nx3-grid2:nx3-1) = data(0:nx1-1,0:grid2-1)
         dataplot(0:nx1-1,0:nx3-1-grid2) = data(0:nx1-1, grid2:nx3-1)
      endif
      if(phiplt eq !dpi) then dataplot = data
      
   endif else dataplot = data
   
   if not keyword_set(plotrange) then begin
      temp = min(abs(xrange(0) - rplot),r1)
      temp = min(abs(xrange(1) - rplot),r2)
      
      temp = min(abs(yrange(0) - azi1),y1)
      temp = min(abs(yrange(1) - azi1),y2)
      
      plotrange0 = [min(dataplot(r1:r2,y1:y2)),max(dataplot(r1:r2,y1:y2))]
   endif else begin
      plotrange0 = plotrange
   endelse
   levels = plotrange0(0) + (plotrange0(1)-plotrange0(0))*dindgen(32)/31d0
   
   name2 = string(i, format='(I03)')
   title1 = textoidl('t=')+strcompress(string(time(i), format='(f5.1)'),/remove_all)+textoidl('P_0')
   title  = title1
   
   if not keyword_set(cart) then begin   
      set_plot,'ps'
      device, filename=strcompress('pdisk_sigma_'+name2+'.ps',/remove_all) $
              ,/color, bits_per_pixel=8, xsize=12, ysize=14
      contour,dataplot,rplot,azi1,/fill,levels=levels, xstyle=1, $
              xtitle=xtitle, ytitle=ytitle,charsize=1.5,xmargin=[7.0,6.0], ymargin=[3.5,2], $
              xtickinterval=xtickinterval,xrange=xrange, yrange=yrange, title=title
      colorbar, position=[0.865, 0.15, 0.90, 0.90],/vertical,/right,range=plotrange0,format='(f5.2)'
      device,/close
   endif else begin
      
      for jj=0, dnx1-1 do begin
         y = yaxis(jj)
         for ii=0, dnx1-1 do begin
            x = xaxis(ii)           
            r_t = sqrt(x^2 + y^2)
            azi_t = pltphi(x, y)            
            if( (r_t ge xrange(0)) and (r_t le Raxis(nrad-1)) )then begin
               
               temp = min(abs(Raxis - r_t),   x0)
               temp = min(abs(azi - azi_t),  y0)
               
               ip = x0 + (r_t - Raxis(x0))/(dx/2.0) + 0.0
               jp = y0 + (azi_t - azi(y0))/(dazi/2.0) + 0.0
               
               dataxy(ii,jj) = bilinear(dataplot, ip, jp)
            endif else begin
               dataxy(ii,jj) = 1d10
            endelse
         endfor
      endfor
      
      set_plot, 'ps'
      device, filename=strcompress('pdiskxy_sigma_'+name2+'.ps',/remove_all) $
              ,/color, bits_per_pixel=8,xsize=14, ysize=12
      contour, dataxy, xaxis, yaxis, /isotropic,/fill,levels=levels,title=title $
               ,ymargin=[2,2],xmargin=[4,4], charsize=1., xtickinterval=xtickinterval, ytickinterval=xtickinterval, xrange=[-1,1]*xrange(1), yrange=[-1,1]*xrange(1), xstyle=1, ystyle=1
      colorbar, position=[0.85, 0.09, 0.9, 0.91],/vertical,/right,range=plotrange0,format='(f5.2)'
      device,/close
   endelse
   
endfor

end
