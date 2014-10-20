function aziavg, array ;array is 3D in rad, theta, phi spherical 

COMMON PLUTO_GRID,  nx1,nx2,nx3,x1,x2,x3,$
                    dx1,dx2,dx3, geometry, xpos, ypos,$
                    AMRLevel, AMRBoxes

array_2d = dblarr(nx1, nx2)

for i=0, nx1-1 do begin
   for j=0, nx2-1 do begin
      array_2d(i,j) = mean(array(i,j,*))
   endfor
endfor

return, array_2d ; this is 2D in rad, theta spherical 
end

function convert_to_Rz, array, Raxis, zaxis ; array is output from above
  
COMMON PLUTO_GRID,  nx1,nx2,nx3,x1,x2,x3,$
                    dx1,dx2,dx3, geometry, xpos, ypos,$
                    AMRLevel, AMRBoxes

ntheta    = n_elements(zaxis)
nR        = n_elements(Raxis)
array_out = dblarr(nR, ntheta)

rad   = x1
theta = x2

dr    = dx1(0)
dth   = dx2(0)

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
         
         array_out(ii,jj) = bilinear(array, ip, jp)
      endif else array_out(ii,jj) = 0.0
      
   endfor
endfor

return, array_out
end

function zint, array2d, zaxis;array is output from above (2D array in R, z coordinates)

COMMON PLUTO_GRID,  nx1,nx2,nx3,x1,x2,x3,$
                    dx1,dx2,dx3, geometry, xpos, ypos,$
                    AMRLevel, AMRBoxes

nR        = n_elements(array2d(*,0))
ntheta    = n_elements(array2d(0,*))

array_1d = dblarr(nR)

for i=0, nR-1 do begin
   array_1d(i) = array2d(i, 0);; int_tabulated(zaxis(0:5), array2d(i,0:5))
endfor

return, array_1d
end


pro pdisk_alpha, range=range, start=start, finish=finish, smallq=smallq, smallh=smallh, r0=r0, dvrad=dvrad

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
if not keyword_set(smallq) then smallq = 0.0
if not keyword_set(range) then range=[0.8,1.2]

pload, 0

rad   = x1
theta = x2
phi   = x3 

nrad  = nx1
ntheta= nx2
nphi  = nx3
time  = t

;assume uniform spacing in r and theta
dr  = rad(1) - rad(0)
dth = theta(1) - theta(0)

;construct R, z axis
zmin = 0.0
zmax = r0*tan(!dpi/2.0 - theta(0))
dz   = (zmax - zmin)/(ntheta-1.0)
zaxis= zmin + dindgen(ntheta)*dz

temp = min(abs(rad - range(0)*r0), r1)
temp = min(abs(rad - range(1)*r0), r2)
Rmin = rad(r1)
Rmax = rad(r2)
nR   = r2 - r1 + 1
dbR   = (Rmax - Rmin)/(nR - 1.0)
Raxis= Rmin + dindgen(nR)*dbR


data       = dblarr(3,finish-start+1)

vR  = dblarr(nx1, nx2, nx3)
csq = dblarr(nx1, nx2)

rho1d         = dblarr(nR)
rho_csq1d     = dblarr(nR)
rho_vR_vphi1d = dblarr(nR)
rho_vR1d      = dblarr(nR)
rho_vphi1d    = dblarr(nR)
alpha1d       = dblarr(nR)

if keyword_set(smallh) then begin
   for j=0, ntheta-1 do begin
      for i=0, nrad-1 do begin
         bigR = rad(i)*sin(theta(j))
         csq(i,j) = smallh*smallh*(bigR/r0)^(-smallq)
      endfor
   endfor
endif


;get initial cylindrical radial velocity
alpha_init = dblarr(nx1, nx2)
for k=0, nphi-1 do begin
   for i=0, nrad-1 do begin
      vR(i,*,k) = vx1(i,*,k)*sin(theta(*)) + vx2(i,*,k)*cos(theta(*))
   endfor
endfor

omega1d = vx3(*,nx2-1,0)/x1
domega1d = deriv(alog(x1), alog(omega1d))

for j=0, nx2-1 do begin
   for i=0, nx1-1 do begin
      bigR = rad(i)*sin(theta(j))
      alpha_init(i,j) = vR(i,j,0)*bigR^(-0.5) ;second factor is keplerian speed
      alpha_init(i,j)/= csq(i,j)
   endfor
endfor

alpha1d_0   = alpha_init(*,nx2-1,0)/domega1d(*)

print, 'background viscous alpha is=', mean(alpha1d_0(r1:r2))

stop

for n=start, finish do begin
   data(0, n-start)     = t(n)

;read data
   pload, n,/silent
   
;get cylindrical radial velocity
   for k=0, nphi-1 do begin
      for i=0, nrad-1 do begin
         vR(i,*,k) = vx1(i,*,k)*sin(theta(*)) + vx2(i,*,k)*cos(theta(*))
      endfor
   endfor

;azi average of density, then vertically integrate
   result       = aziavg(rho)
   result_Rz    = convert_to_Rz(result, Raxis, zaxis)
   rho1d        = zint(result_Rz, zaxis)


;azi average of density*csq, then vertically integrate
   result       = aziavg(rho)*csq
   result_Rz    = convert_to_Rz(result, Raxis, zaxis)
   rho_csq1d    = zint(result_Rz, zaxis)

;azi average of rho*vR, then vertically integrate
   result    = aziavg(rho*vR)
   result_Rz = convert_to_Rz(result, Raxis, zaxis)
   rho_vR1d  = zint(result_Rz, zaxis)

;azi average of rho*vphi, then vertically integrate
   result    = aziavg(rho*vx3)
   result_Rz = convert_to_Rz(result, Raxis, zaxis)
   rho_vphi1d  = zint(result_Rz, zaxis)

;azi average of rho*vR*vphi, then vertically integrate
   result         = aziavg(rho*vR*vx3)
   result_Rz      = convert_to_Rz(result, Raxis, zaxis)
   rho_vR_vphi1d  = zint(result_Rz, zaxis)

;compute reynolds stress and normalize it
   alpha1d = rho_vR_vphi1d - (rho_vphi1d/rho1d)*rho_vR1d
   alpha1d/= rho_csq1d

;alpha relative to background
   alpha1d -= alpha1d_0

;average over radius to get single value
   data(1, n-start) = mean(alpha1d)
   print, 'max alpha =', max(alpha1d)
endfor

;computing running time average
for k=1, finish-start do begin
   data(2, k) = int_tabulated(data(0,0:k), data(1, 0:k))
   data(2, k) /= data(0,k) - data(0,0)
endfor

openw,1, 'pdisk_alpha.dat'
for i=0, finish-start do begin
   printf, 1, data(0:2,i), format='(3(e22.15,x))'
endfor
close,1

end
