;+
;
; NAME:       curl
;
; AUTHOR:    Andrea Mignone (mignone@to.astro.it)
;
; PURPOSE:      compute the curl of a vector field
;
; SYNTAX:      array = curl(u,v,x1,x2,dx1,dx2,geo=geo)
;
; DESCRIPTION: curl returns an array of the same size
;              of u (or v) containing the curl
;              of the vector field (u,v) computed as:
;
;              d(u)/dy - d(v)/dx
;
;              The input parameters are:
;
;
;               u   = first  vector component (2-D array)
;               v   = second vector component (2-D array)
;               x1  = 1-D array of n1 points containing the x1 zone centers
;               x2  = 1-D array of n2 points containing the x2 zone centers
;               dx1 = 1-D array of n1 points containing mesh widths
;               dx2 = 1-D array of n2 points containing mesh widths
;              geo = specify in which geometry the divergence
;                     should be computed. Options are:
;                     'cartesian', 'spherical', 'cylindrical', 'polar'
;
;  EXAMPLES
;
;    #1 Compute and display the curl of the initial magnetic field (toroidal
;       current) in cartesian coordinates:
;    IDL> display, curl(b1(0), b2(0), x1, x2, dx1, dx2, geo="cartesian"), /vbar
;
;
; HISTORY:      Sept 22, 2003
;
;-

function curl, u1, u2, x1, x2, dx1, dx2, $
              geo   = geo

if NOT KEYWORD_SET(geo) THEN geo='cartesian'

sz = size(u1)

ndim = sz[0]
n1   = sz[1]
n2   = sz[2]

du1 = fltarr(n1,n2)
du2 = fltarr(n1,n2)
A1  = fltarr(n1)
A2  = fltarr(n2)
dV1 = fltarr(n1,n2)
dV2 = fltarr(n1,n2)

; ------------------------------------------------
;  define area and volume elements for the
;  different coordinate systems
; ------------------------------------------------

IF (STRCMP(geo,'cartesian',3,/FOLD_CASE)) THEN BEGIN

  dq  = 1.0/dx2
  for j = 1,n2 - 2 do begin
    du1(*,j) = 0.5*(u1(*,j+1) - u1(*,j-1))*dq(j)
  endfor
  dq = 1.0/dx1
  for i = 1,n1 - 2 do begin
    du2(i,*) = 0.5*(u2(i+1,*) - u2(i-1,*))*dq(i)
  endfor

;
; Interpolate at boundaries
;

  j = 0
  for i = 0,n1-1 do begin
    du1(i,j) = interpol(du1(i,j+1:j+4),x2(j+1:j+4),x2(j))
  endfor

  j = n2-1
  for i = 0,n1-1 do begin
    du1(i,j) = interpol(du1(i,j-4:j-1),x2(j-4:j-1),x2(j))
  endfor

  i = 0 
  for j = 0,n2-1 do begin
    du2(i,j) = interpol(du2(i+1:i+4,j),x1(i+1:i+4),x1(i))
  endfor

  i = n1-1
  for j = 0,n2-1 do begin
    du2(i,j) = interpol(du2(i-4:i-1,j),x1(i-4:i-1),x1(i))
  endfor
ENDIF

IF (STRCMP(geo,'cylindrical',3,/FOLD_CASE)) THEN BEGIN
  print," ! curl does not work in cylindrical     coordinates"
  stop
ENDIF

IF (STRCMP(geo,'polar',3,/FOLD_CASE)) THEN BEGIN

  dq  = 1.0/dx2
  FOR j = 1,n2 - 2 DO BEGIN
    du1(*,j) = 0.5*(u1(*,j+1) - u1(*,j-1))*dq(j)
  endfor
  dq = 1.0/dx1
  for i = 1,n1 - 2 do begin
    du2(i,*) = 0.5*(u2(i+1,*)*x1(i+1) - u2(i-1,*)*x1(i-1))*dq(i)
  endfor

  for i=1,n1-2 do begin
     du1(i,*) = du1(i,*)/x1(i)
     du2(i,*) = du2(i,*)/x1(i)
  endfor

  ; ---------------------------
  ;  Interpolate at boundaries  
  ; ---------------------------

  j = 0
  for i = 0,n1-1 do begin
    du1(i,j) = interpol(du1(i,j+1:j+4),x2(j+1:j+4),x2(j))
  endfor

  j = n2-1
  for i = 0,n1-1 do begin
    du1(i,j) = interpol(du1(i,j-4:j-1),x2(j-4:j-1),x2(j))
  endfor

  i = 0
  for j = 0,n2-1 do begin
    du2(i,j) = interpol(du2(i+1:i+4,j),x1(i+1:i+4),x1(i))
  endfor

  i = n1-1
  for j = 0,n2-1 do begin
    du2(i,j) = interpol(du2(i-4:i-1,j),x1(i-4:i-1),x1(i))
  endfor
ENDIF; geo == POLAR

IF (STRCMP(geo,'spherical',3,/FOLD_CASE)) THEN BEGIN
  print," ! curl does not work in spherical coordinates"
  stop
ENDIF

; ------------------------------------------------
;              Make Curl
; ------------------------------------------------

return, (du2-du1)
end
