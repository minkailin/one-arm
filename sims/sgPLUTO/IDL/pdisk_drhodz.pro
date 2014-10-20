pro pdisk_drhodz, start=start, finish=finish, smallq=smallq, smallh=smallh, r0=r0, omit=omit

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

pload, 0

rad   = x1
theta = x2
phi   = x3 

nrad  = nx1
ntheta= nx2
nphi  = nx3
time  = t

temp = min(abs(rad - r0), rzero)
temp = min(abs(rad - 0.8*r0), r1)
temp = min(abs(rad - 1.2*r0), r2)

if not keyword_set (omit) then omit = 0
t1   = omit
t2   = ntheta-2 - omit

rho_2d     = dblarr(nrad, ntheta)
drho_dtheta= dblarr(nrad, ntheta, nphi)
drho_dr    = dblarr(nrad, ntheta, nphi)

bigH        = dblarr(nrad, ntheta)

data       = dblarr(5,finish-start+1)

if keyword_set(smallh) then begin
for j=0, ntheta-1 do begin
for i=r1, r2 do begin
bigR   = rad(i)*sin(theta(j))
omegak = bigR^(-1.5)
csq    = sqrt(smallh*smallh*(bigR/r0)^(-smallq))
bigH(i,j) = csq/omegak
endfor
endfor
endif

buff = 5

for n=start, finish do begin
data(0, n-start)     = t(n)

;read data
   pload, n,/silent

;azimuthal average of density field 
for j=0, ntheta-1 do begin
for i=r1-buff, r2+buff do begin
rho_2d(i,j) = mean(rho(i,j,*))
endfor
endfor

;get  density/<density>_phi 
for k=0, nphi-1 do begin
rho(r1-buff:r2+buff,*,k) = rho(r1-buff:r2+buff,*,k)/rho_2d(r1-buff:r2+buff,*)
endfor

;radial and theta gradient
for k=0, nphi-1 do begin

for j=0, ntheta-1 do begin
drho_dr(r1-buff:r2+buff,j,k) = cos(theta(j))*deriv(rad(r1-buff:r2+buff), rho(r1-buff:r2+buff,j,k))
endfor

for i=r1-buff, r2+buff do begin
drho_dtheta(i,*,k) =sin(theta(*))*deriv(theta, rho(i,*,k))/rad(i)
endfor

endfor

;get absolute value of d(Delta rho)/dz
rho = abs( drho_dr - drho_dtheta )

;multiply by bigH to non-dimensionalize
for k=0, nphi-1 do begin
rho(r1-buff:r2+buff,*,k) *= bigH(r1-buff:r2+buff,*)
endfor

;average over theta and rad to get single value

data(1, n-start) = mean(rho(r1:r2, t1:t2))

data(2, n-start) = mean(data(1, 0:n-start)) ;running-time average 


;max values 
data(3, n-start) = max(rho(r1:r2, t1:t2))

data(4, n-start) = mean(data(3, 0:n-start)) ;running-time average


print, 'max[drhodz] at ', string(n,format='(I03)'),' is', max(rho(r1:r2,t1:t2))

endfor

openw,1, 'pdisk_drhodz.dat'
for i=0, finish-start do begin
printf, 1, data(0:4,i), format='(5(e22.15,x))'
endfor
close,1

end
