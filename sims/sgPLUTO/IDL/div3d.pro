;+
;
; NAME:       div3d
;
; AUTHOR:       P. Tzeferacos, A. Mignone
;
; PURPOSE:      to compute the divergence of a vector field
;
; SYNTAX:       array = div(u,v,w,x1,x2,x3,dx1,dx2,dx3,/geo)
;
; DESCRIPTION:  The input parameters are:
;
;
;               u = first  vector component (3-D array)
;               v = second vector component (3-D array)
;               w = third  vector component (3-D array)
;               x1 = 1-D array of n1 points containing the x1 zone centers
;               x2 = 1-D array of n2 points containing the x2 zone centers
;               x3 = 1-D array of n3 points containing the x3 zone centers
;               dx1 = 1-D array of n1 points containing mesh widths
;               dx2 = 1-D array of n2 points containing mesh widths
;               dx3 = 1-D array of n3 points containing mesh widths
;              /geo = specify in which geometry the divergence
;                     should be computed. Options are:
;
;                     /cartesian, /spherical, /polar
;
; HISTORY:      July 12, 2007
;
;-

function div3d, u1, u2, u3, x1, x2, x3, dx1, dx2, dx3, geo = geo


if NOT KEYWORD_SET(geo) THEN geo='cartesian'

sz = size(u1)

ndim  = sz[0]
n1    = sz[1]
n2    = sz[2]
n3    = sz[3]

du1dx1   = fltarr(n1,n2,n3)
du2dx2   = fltarr(n1,n2,n3)
du3dx3   = fltarr(n1,n2,n3)


A1       = fltarr(n1,n2)
B1       = fltarr(n1,n2,n3)
A2       = fltarr(n1,n2,n3)
A3       = fltarr(n1,n2,n3)

div      = fltarr(n1,n2,n3)






;
;calculate derivatives
;

;
;second order
;

;dx1

dq = 1.0/dx1

for i = 1,n1 - 2 do begin
  du1dx1(i,*,*) = 0.5*(u1(i+1,*,*) - u1(i-1,*,*))*dq(i)
endfor


;dx2

dq  = 1.0/dx2

for j = 1,n2 - 2 do begin
  du2dx2(*,j,*) = 0.5*(u2(*,j+1,*) - u2(*,j-1,*))*dq(j)
endfor


;dx3

dq  = 1.0/dx3

for k = 1,n3 - 2 do begin
  du3dx3(*,*,k) = 0.5*(u3(*,*,k+1) - u3(*,*,k-1))*dq(k)
endfor


;
; Boundary derivatives (first order)
;



;dx1
dq = 1.0/dx1

ibeg = 0
  du1dx1(ibeg,*,*) = (u1(ibeg+1,*,*) - u1(ibeg,*,*))*dq(ibeg)

iend = n1-1
  du1dx1(iend,*,*) = (u1(iend,*,*) - u1(iend-1,*,*))*dq(iend)



;dx2
dq  = 1.0/dx2

jbeg = 0
  du2dx2(*,jbeg,*) = (u2(*,jbeg+1,*) - u2(*,jbeg,*))*dq(jbeg)

jend = n2-1
  du2dx2(*,jend,*) = (u2(*,jend,*) - u2(*,jend-1,*))*dq(jend)


;dx3
dq  = 1.0/dx3

kbeg = 0
  du3dx3(*,*,kbeg) = (u3(*,*,kbeg+1) - u3(*,*,kbeg))*dq(kbeg)

kend = n3-1
  du3dx3(*,*,kend) = (u3(*,*,kend) - u3(*,*,kend-1))*dq(kend)


IF (STRCMP(geo,'cartesian',3,/FOLD_CASE)) THEN BEGIN

div = du1dx1 + du2dx2 + du3dx3

ENDIF

IF (STRCMP(geo,'polar',3,/FOLD_CASE)) THEN BEGIN

;------------------------------
;create      1 / x1 == A2
;------------------------------

for i = ibeg,iend do begin
A2(i,*,*) = 1.0/x1(i)
endfor

div = A2*u1 + du1dx1 + A2*du2dx2 + du3dx3


ENDIF

IF (STRCMP(geo,'spherical',3,/FOLD_CASE)) THEN BEGIN

;------------------------------
;Create 1 /( x1 sin(x2) ) == A1
;       1 / x1            == A2
;       cos(x2)           == A3
;------------------------------

for i = ibeg,iend do begin
A1(i,*)   = 1.0/(x1(i)*sin(x2(*)))
A2(i,*,*) = 1.0/x1(i)
endfor

for j = jbeg,jend do begin
A3(*,j,*) = cos(x2(j))
endfor

;
;make A1 a 3D array
;

for k = kbeg,kend do begin
B1(*,*,k)=A1(*,*)
endfor

div = 2.0*A2*u1 + du1dx1 + B1*du2dx2*A3 + A2*du2dx2 + B1*du3dx3


ENDIF

return, div

end









