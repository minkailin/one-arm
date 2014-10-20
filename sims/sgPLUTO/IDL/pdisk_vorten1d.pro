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

omegaz = omegar - omegatheta ; geometric factors already incoroporated. 

return, omegaz
end


pro pdisk_vorten1d, loc=loc, start=start, finish=finish, xrange=xrange, yrange=yrange, $
         r0=r0, legend=legend, label=label, nopert=nopert, smth=smth $
        ,ytickinterval=ytickinterval, smallh=smallh, mp=mp, xtickinterval=xtickinterval 

COMMON PLUTO_GRID,  nx1,nx2,nx3,x1,x2,x3,$
                    dx1,dx2,dx3, geometry, xpos, ypos,$
                    AMRLevel, AMRBoxes;  ** Chombo data structure **
                                      ;  ** loaded when HDF5LOAD is called **

COMMON PLUTO_VAR, NVAR, rho, vx1, vx2, vx3, $
                             bx1, bx2, bx3, $
                             Ax1, Ax2, Ax3, $
                             bx1s, bx2s, bx3s,$
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



if not keyword_set(finish) then finish=start 
if not keyword_set(zslice) then zslice= 0.99
if not keyword_set(r0) then r0=1.0
if keyword_set(nopert) then begin
   ytitle=textoidl('<\eta_z>_\phi/<\eta_z>_\phi(r_0)')
endif else begin
   ytitle = textoidl('<\eta_z>_\phi/\eta_z(t=0) - 1')
endelse

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

den_2d  = dblarr(nrad, ntheta)
vort_2d = dblarr(nrad, ntheta)

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

densxy = dblarr(nR, ntheta)
vortxy = dblarr(nR, ntheta)

den_1d   = dblarr(nR)
vort_1d  = dblarr(nR) 
data1d   = dblarr(nR)

vortz       = get_vorticity(vx1, vx2, vx3, x1, x2, x3)
for jj=0, ntheta-1 do begin
   for ii=0, nrad-1 do begin 
      den_2d(ii,jj) = mean(rho(ii,jj,*))      
      vort_2d(ii,jj) = mean(vortz(ii,jj,*))  
   endfor
endfor

for jj=0, ntheta-1 do begin
   z = zaxis(jj)
   for ii=0, nR-1 do begin
      R = Raxis(ii)
      
;   recpol, R, z, r_t, psi_t
      r_t = sqrt(R*R + z*z)
      th_t = atan(R, z)
      
      if(th_t ge theta(0)) then begin
         
         temp = min(abs(rad - r_t),   x0)
         temp = min(abs(theta - th_t),y0)
         
         ip = x0 + (r_t - rad(x0))/(dr/2.0) + 0.0
         jp = y0 + (th_t - theta(y0))/(dth/2.0) + 0.0
   
         densxy(ii,jj) = bilinear(den_2d, ip, jp)
         vortxy(ii,jj) = bilinear(vort_2d, ip, jp)
      endif else begin
         densxy(ii,jj) = -1
         vortxy(ii,jj) = -1
      endelse

   endfor
endfor

for ii=0, nR-1 do begin
   mask = where(densxy(ii,*) gt 0.0)
   den_1d(ii) = int_tabulated(zaxis(mask), densxy(ii,mask))
   vort_1d(ii) = int_tabulated(zaxis(mask), vortxy(ii,mask))
endfor

data1d_0 = vort_1d/den_1d
if keyword_set(smth) then data1d_0 = smooth(data1d_0, smth)

for i=start, finish do begin
   pload, i,/silent
   
   vortz       = get_vorticity(vx1, vx2, vx3, x1, x2, x3)
   for jj=0, ntheta-1 do begin
      for ii=0, nrad-1 do begin
         den_2d(ii,jj) = mean(rho(ii,jj,*)) 
         vort_2d(ii,jj) = mean(vortz(ii,jj,*))  
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
            
            densxy(ii,jj) = bilinear(den_2d, ip, jp)
            vortxy(ii,jj) = bilinear(vort_2d, ip, jp)
         endif else begin
            densxy(ii,jj) = -1
            vortxy(ii,jj) = -1
         endelse
         
      endfor
   endfor  
  
   for ii=0, nR-1 do begin
      mask = where(densxy(ii,*) gt 0.0)
      den_1d(ii) = int_tabulated(zaxis(mask), densxy(ii, mask))
      vort_1d(ii) = int_tabulated(zaxis(mask), vortxy(ii,mask))
   endfor
   
   data1d = vort_1d/den_1d
if keyword_set(smth) then   data1d = smooth(data1d, smth)
   ;; data1d = vortxy(*,15)/den_1d

   if keyword_set(nopert) then begin
      temp = min(abs(Raxis-r0), rgrid)
      data1d /= data1d(rgrid)
   endif else begin
      data1d /= data1d_0
      data1d -= 1.0
   endelse

   print, 'max deta_z', max(abs(data1d))

   name2 = string(i, format='(I03)')
   title = textoidl('t=')+strcompress(string(time(i), format='(f5.1)'),/remove_all)+textoidl('P_0')
   
   set_plot, 'ps'
   device, filename=strcompress('pdisk_vorten1d_'+name2+'.ps',/remove_all) $
          ,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches
   plot, Raxis, data1d,xmargin=[8.5,1.5],ymargin=[3.5,1.5], ystyle=1, xstyle=1  $
        ,charsize=1.5, thick=4, xrange=xrange, title=title, xtitle=textoidl('R/r_0'), yrange=yrange,$
        linestyle = 0, ytitle =ytitle, xtickinterval=xtickinterval, ytickinterval=ytickinterval

;   oplot, [0.98,0.98], [-10,10],linestyle=2
   device,/close 


   file =strcompress('pdisk_vorten1d_'+name2+'.dat',/remove_all)
   openw,1, file
   for ii=0, nR-1 do begin
      printf,1,Raxis(ii),data1d(ii),format='(2(e22.15,x))'
   endfor
   close,1

  endfor
end
