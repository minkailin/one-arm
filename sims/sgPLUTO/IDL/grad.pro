;+
;
; NAME:        grad
;
; AUTHOR:      A. Mignone
;
; PURPOSE:     to compute the gradient of a scalar
;
; SYNTAX:      array = grad(u,v,x1,x2,dx1,dx2,/geo)
;
; DESCRIPTION:
;
;               phi = scalar function (2-D array)
;               x1  = 1-D array of n1 points containing the x1 zone centers
;               x2  = 1-D array of n2 points containing the x2 zone centers
;               dx1 = 1-D array of n1 points containing mesh widths
;               dx2 = 1-D array of n2 points containing mesh widths
;              /geo = specify in which geometry the divergence
;                     should be computed. Options are:
;
;                     /cartesian, /spherical, /cylindrical, /polar
;
; HISTORY:      Sept 22, 2003
;
;-

function grad, phi, x1, x2, dx1, dx2, $
               polar = polar


sz = size(phi)

ndim     = sz[0]
n1       = sz[1]
n2       = sz[2]
grad_phi = fltarr(n1,n2,2)


; ------------------------------------------------
;  define area and volume elements for the
;  different coordinate systems
; ------------------------------------------------

h2      = fltarr(n1,n2)
h2(*,*) = 1.0

IF (KEYWORD_SET(polar)) THEN BEGIN
  FOR j = 0, n2-1 DO BEGIN
    h2(*,j) = x1(*)
  ENDFOR
ENDIF


; ------------------------------------------------
;              Make gradient
; ------------------------------------------------

FOR j = 0, n2 - 1 DO BEGIN
 scrh = phi(*, j)
 grad_phi(*, j, 0) = deriv(x1, scrh)
ENDFOR
FOR i = 0, n1 - 1 DO BEGIN
 scrh = phi(i, *)
 grad_phi(i, *, 1) = deriv(x2, scrh)/h2(i,*)
ENDFOR

return, grad_phi
end
