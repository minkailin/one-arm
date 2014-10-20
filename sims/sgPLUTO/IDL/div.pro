;+
;
; NAME:       div
;
; AUTHOR:     A. Mignone
;
; PURPOSE:     to compute the divergence of a vector field
;
; SYNTAX:       array = div(u,v,x1,x2,dx1,dx2,/geo)
;
; DESCRIPTION: div returns an array of the same size
;              of u (or v) containing the divergence
;              of the vector field (u,v) computed as:
;
;              d(Au)/dV + d(Av)/dV
;
;              The input parameters are:
;
;
;               u = first  vector component (2-D array)
;               v = second vector component (2-D array)
;               x1 = 1-D array of n1 points containing the x1 zone centers
;               x2 = 1-D array of n2 points containing the x2 zone centers
;               dx1 = 1-D array of n1 points containing mesh widths
;               dx2 = 1-D array of n2 points containing mesh widths
;              /geo = specify in which geometry the divergence
;                     should be computed. Options are:
;
;                     /cartesian, /spherical, /cylindrical, /polar
;
;  EXAMPLES
;
;    #1 Compute and display the divergence of the magnetic field in cartesian
;       coordinates of the 1st output file:
;    IDL> display, div(b1(1), b2(1), x1, x2, dx1, dx2, /cartesian), /vbar
;
;    #2 Compute and display the initial divergence of the velocity in
;       cylindrical coordinates:
;    IDL> display, div(v1(0), v2(0), x1, x2, dx1, dx2, /cylindrical), /vbar;
;
;
; HISTORY:      Sept 22, 2003
;
;-

function div, u1, u2, x1, x2, dx1, dx2, $
              cartesian   = cartesian,$
              cylindrical = cylindrical,$
              spherical   = spherical, $
              polar       = polar


sz = size(u1)


ndim = sz[0]
n1   = sz[1]
n2   = sz[2]

du1 = fltarr(n1,n2)
du2 = fltarr(n1,n2)

A1 = fltarr(n1)
A2 = fltarr(n2)

dV1 = fltarr(n1,n2)
dV2 = fltarr(n1,n2)


; ------------------------------------------------
;  define area and volume elements for the
;  different coordinate systems
; ------------------------------------------------

IF (KEYWORD_SET(cartesian)) THEN BEGIN
  A1(*) = 1.0
  A2(*) = 1.0
  dV1   = dx1 # A2
  dV2   = A1  # dx2
ENDIF

IF (KEYWORD_SET(cylindrical)) THEN BEGIN
  A1    = x1
  A2(*) = 1
  dV1   = (x1*dx1) # A2
  for i = 0,n1-1 do dV2(i,*) = dx2(*)
ENDIF

IF (KEYWORD_SET(polar)) THEN BEGIN
  A1    = x1
  A2(*) = 1
  dV1   = x1(*) # A2
  dV2   = x1 # dx2
ENDIF

IF (KEYWORD_SET(spherical)) THEN BEGIN
  A1  = x1*x1
  A2  = sin(x2)
  for j = 0,n2-1 do dV1(*,j) = A1*dx1
  dV2 = x1 # (sin(x2)*dx2)
ENDIF

; ------------------------------------------------
;              Make divergence
; ------------------------------------------------

for i = 1,n1 - 2 do begin
  du1(i,*) = 0.5*(A1(i+1)*u1(i+1,*) - A1(i-1)*u1(i-1,*)) $
                          /dV1(i,*)
endfor
for j = 1,n2 - 2 do begin
  du2(*,j) = 0.5*(A2(j+1)*u2(*,j+1) - A2(j-1)*u2(*,j-1)) $
                          /dV2(*,j)
endfor

;
; Interpolate at boundaries
;

i = 0
for j = 0,n2-1 do begin
  du1(i,j) = interpol(du1(i+1:i+4,j),x1(i+1:i+4),x1(i),/quadratic)
endfor

i = n1-1
for j = 0,n2-1 do begin
  du1(i,j) = interpol(du1(i-4:i-1,j),x1(i-4:i-1),x1(i),/quadratic)
endfor


j = 0
for i = 0,n1-1 do begin
  du2(i,j) = interpol(du2(i,j+1:j+4),x2(j+1:j+4),x2(j),/quadratic)
endfor

j = n2-1
for i = 0,n1-1 do begin
  du2(i,j) = interpol(du2(i,j-4:j-1),x2(j-4:j-1),x2(j),/quadratic)
endfor

return, (du1+du2)
end
