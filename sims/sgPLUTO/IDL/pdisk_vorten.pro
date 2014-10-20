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


pro pdisk_vorten, start=start, finish=finish, xrange=xrange, yrange=yrange, scale=scale, label=label,colxrange=colxrange, $
         r0=r0, log=log, smth=smth,nopert=nopert,ct=ct,ytickinterval=ytickinterval, xtickinterval=xtickinterval,inv=inv, plotrange=plotrange

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
if not keyword_set(r0) then r0=1.0
if not keyword_set(ct) then ct=5
if not keyword_set(yrange) then yrange=[0,1]
if not keyword_set(scale) then scale=1.0

pload, 0 

nrad  = nx1
ntheta= nx2
nphi  = nx3

rad   = x1
theta = x2
phi   = x3

azi1=x3/(2.0*!dpi)

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

vorten2d=dblarr(nR, nx3)

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
if keyword_set(inv) then data1d_0 = 1.0/data1d_0

if keyword_set(smth) then data1d_0 = smooth(data1d_0, smth);axisymmetric vortensity 

for i=start, finish do begin
   pload, i,/silent
   
   vortz       = get_vorticity(vx1, vx2, vx3, x1, x2, x3)

for kk=0, nx3-1 do begin

    den_2d(*,*) = rho(*,*,kk)
    vort_2d(*,*) = vortz(*,*,kk) 
   
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
   if keyword_set(inv) then data1d = 1.0/data1d

if keyword_set(smth) then   data1d = smooth(data1d, smth)


   if not keyword_set(nopert) then begin
      data1d /= data1d_0
   if not keyword_set(log) then   data1d -= 1.0
   endif

   if keyword_set(log) then data1d = alog10(data1d)

   vorten2d(0:nR-1, kk) = data1d(*) 

endfor      

    dataplot = vorten2d*scale


   if not keyword_set(plotrange) then begin
      temp = min(abs(xrange(0) - Raxis),r1)
      temp = min(abs(xrange(1) - Raxis),r2)

      if keyword_set(colxrange) then begin
      temp = min(abs(colxrange(0) - Raxis),cr1)
      temp = min(abs(colxrange(1) - Raxis),cr2)
      endif else begin
      cr1 = r1
      cr2 = r2
      endelse

      temp = min(abs(yrange(0) - azi1),y1)
      temp = min(abs(yrange(1) - azi1),y2)

      plotrange0 = [min(dataplot(r1:r2,y1:y2)),max(dataplot(cr1:cr2,y1:y2))]
   endif else begin
      plotrange0 = plotrange
   endelse
   levels = plotrange0(0) + (plotrange0(1)-plotrange0(0))*dindgen(32)/31d0


   name2 = string(i, format='(I03)')
   title1 = textoidl('t=')+strcompress(string(time(i), format='(f5.1)'),/remove_all)+textoidl('P_0')
   title  = title1 ;+ textoidl(', tan\psi=')+height_string+'h'

   loadct,ct,/silent
   set_plot,'ps'
   device, filename=strcompress('pdisk_vorten_'+name2+'.ps',/remove_all) $
           ,/color, bits_per_pixel=8, xsize=12, ysize=14
   contour,dataplot, Raxis, azi1,/fill,levels=levels, xstyle=1, $
           ytitle=textoidl('\phi/2\pi'), xtitle=textoidl('R/r_0'),charsize=1.5,xmargin=[7.0,6.0], ymargin=[3.5,2], $
           xtickinterval=xtickinterval,ytickinterval=ytickinterval, xrange=xrange, yrange=yrange, title=title
   colorbar, position=[0.852, 0.15, 0.90, 0.90],/vertical,/right,range=plotrange0,format='(f5.2)'

   if keyword_set(label) then begin
   xyouts, xrange(0)+(xrange(1)-xrange(0))/20., yrange(0)+(yrange(1)-yrange(0))/10., textoidl(label), charsize=2,charthick=8, color=255
   endif

   device,/close

endfor
end
