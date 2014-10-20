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

pro pdisk_vort1d, start=start, finish=finish, xrange=xrange, yrange=yrange, $
        zslice=zslice, r0=r0, legend=legend, label=label, azislice=azislice, nopert=nopert $
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

vortz_2d = dblarr(nrad, ntheta)
dvortz_2d = dblarr(nrad, ntheta)
vort_mid = dblarr(nrad, nphi)
two_omega= dblarr(nrad, ntheta)

nzslice = ntheta*zslice

zmin = 0.0
zmax = rad(nx1-1)*tan(!dpi/2.0 - theta(0))
dz   = (zmax - zmin)/(ntheta-1.0)
zaxis= zmin + dindgen(ntheta)*dz
print, 'chosen heights are:'
for i=0, n_elements(nzslice)-1 do begin
print, zaxis(nzslice(i))/(smallh*r0)
endfor
print,'in units of H(r_0)'

temp = min(abs(rad - xrange(0)*r0), r1)
temp = min(abs(rad - xrange(1)*r0), r2)
Rmin = rad(r1)
Rmax = rad(r2)
nR   = r2 - r1 + 1
dbR   = (Rmax - Rmin)/(nR - 1.0)
Raxis= Rmin + dindgen(nR)*dbR

data2d = dblarr(nR, ntheta)
data1d = dblarr(nR)

vortz0 = get_vorticity(vx1, vx2, vx3, x1, x2, x3)


for i=start, finish do begin
   pload, i,/silent
   
   vortz = get_vorticity(vx1, vx2, vx3, x1, x2, x3)
   dvortz = vortz/vortz0 - 1.0

   for jj=0, nx2-1 do begin
    for ii=0, nx1-1 do begin
        two_omega(ii,jj) = 2.0*mean(vx3(ii,jj,*))/(x1(ii)*sin(x2(jj))) 
        vortz_2d(ii,jj) = mean(vortz(ii,jj,*))
        vortz_2d(ii,jj)/= two_omega(ii,jj)
        dvortz_2d(ii,jj) = mean(dvortz(ii,jj,*))
    endfor
   endfor   

   if not keyword_set(nopert) then vortz_2d = dvortz_2d

   if not keyword_set(azislice) then begin
      if not keyword_set(nopert) then begin
         ytitle = textoidl('<\Delta(curl v)_z>_\phi')
      endif else begin
         ytitle = textoidl('<(curl v)_z>_\phi/2<\Omega>_\phi')
      endelse
   endif else begin

      if (azislice gt 0.0) then begin
         nphislice = azislice*nphi 
      endif else begin  
         for ii=0, nrad-1 do begin 
            vort_mid(ii,*) = vortz(ii,ntheta-1,*)/two_omega(ii,ntheta-1) - vortz_2d(ii,ntheta-1)
         endfor
         temp = min(vort_mid, minloc)
         ind  = array_indices(vort_mid, minloc)
         nphislice = ind(1)     
      endelse
      print, 'azislice/2.pi=', x3(nphislice)/(2.0*!dpi)
      
      if not keyword_set(nopert) then begin
         ytitle = textoidl('Ro')
         vortz_2d(*,*) = vortz(*,*,nphislice)/two_omega(*,*) - vortz_2d(*,*) 
      endif else begin
         ytitle =  textoidl('(curl v)_z/2<\Omega>_\phi')
         vortz_2d(*,*) = vortz(*,*,nphislice)/two_omega(*,*)
      endelse

   endelse
   
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

   data2d(ii,jj) = bilinear(vortz_2d, ip, jp)
   endif else data2d(ii,jj) = 0.0

   endfor
   endfor  
   result = data2d


   data1d = result(*,nzslice(0))

   name2 = string(i, format='(I03)')
   title = textoidl('t=')+strcompress(string(time(i), format='(f5.1)'),/remove_all)+textoidl('P_0')
if keyword_set(azislice) then begin
   azistring = string(phi(nphislice)/(2.0*!dpi), format='(f4.2)')
   title = title + textoidl(', \phi/2\pi=')+azistring
endif
   set_plot, 'ps'
   device, filename=strcompress('pdisk_vort1d_'+name2+'.ps',/remove_all) $
          ,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches
   plot, Raxis, data1d,xmargin=[8.5,1.5],ymargin=[3.5,1.5], ystyle=0  $
        ,charsize=1.5, xstyle=1,thick=4, xrange=xrange, title=title, xtitle=textoidl('R/r_0'), yrange=yrange,$
        linestyle = 0, ytitle =ytitle, xtickinterval=xtickinterval, ytickinterval=ytickinterval

  for k=1, n_elements(nzslice)-1 do begin
  data1d(*) = result(*,nzslice(k))
  mask = where(data1d ne 0.0)
  oplot, Raxis(mask), data1d(mask), thick=4, linestyle=k
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

  endfor
end
