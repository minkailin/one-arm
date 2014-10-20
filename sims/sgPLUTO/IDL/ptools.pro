; TODO: A.GRAD(B), A.B, A x B

;+
;
; NAME:        PAVERAGE
;
; AUTHOR:      A. Mignone (mignone@ph.unito.it)
;
; SYNTAX:      q_av = PAVERAGE(q, /geometry, dir=dir)
;
; PURPOSE:     compute the volume-average of a 2D or 3D array "q" over the 
;              computational grid defined during the most recent call to the 
;              PLOAD function.
;              On output, "q_av" is either a scalar - if the average is taken
;              over all directions - or an array of reduced size if the 
;              average considers only some of the coordinates.
;              The dimensionality of the array decreases by one each time 
;              a coordinate axis is averaged down.
;
; ARGUMENTS:   
;
;   q          2D or 3D array with (NX1,NX2) or (NX1,NX2,NX3) points, 
;              respectively.
;                      
; KEYWORDS:
;
;   /geometry   Specify the geometry to be adopted. 
;               Options are: /polar, /cylindrical (2D) and /spherical.
;               If none is given, a Cartesian mesh will be assumed.
;
;   dir         an array of integers giving the directions over which 
;               the average has to be taken. 
;               For instance, dir=[0,1,1] will average along x2 and x3 and
;               the output will be a 1D array function of x1.
;              
; EXAMPLES:
;
;  * Example #1: compute the average of density in 3D spherical coordinates,
;
;    IDL> rho_av = PAVERAGE(rho,/spherical)
;
;  * Example #2: compute the average of vorticity along the phi and z 
;                directions in 3D polar coordinates:
;
;    IDL> w = PCURL(v1,v2,/polar)
;    IDL> wr = PAVERAGE(w,dir=[0,1,1],/polar)
;
; LAST MODIFIED:   Jan 29, 2012 by A. Mignone (mignone@ph.unito.it)
;
;-
FUNCTION PAVERAGE, q, polar=polar, spherical=spherical,dir=dir

 COMMON PLUTO_GRID

; --------------------------------------------------------
;   check if the array size is correct
; --------------------------------------------------------

 ndim = SIZE(q,/N_DIMENSIONS)
 sz   = SIZE(q,/DIMENSIONS)
 n123 = [NX1,NX2,NX3]
 FOR idim = 0, ndim-1 DO BEGIN
   IF (sz[idim] NE n123[idim]) THEN BEGIN
     PRINT,"! Scalar dimensions does not match with PLUTO grid"
     RETURN,-1
   ENDIF
 ENDFOR

; --------------------------------------------------------------------------
;  Polar coordinates:    dV = (rp^2-rm^2)/2 * dphi * dz
; --------------------------------------------------------------------------

 IF (KEYWORD_SET(polar)) THEN BEGIN
   rp = x1 + 0.5*dx1 
   rm = x1 - 0.5*dx1 
  
   dVr   = 0.5*(rp^2 - rm^2)
   dVphi = dx2
   dVz   = dx3
   dV = FLTARR(NX1,NX2,NX3)
   FOR k=0,NX3-1 DO dV(*,*,k) = (dVr#dVphi)*dVz(k)

 ENDIF ELSE IF (KEYWORD_SET(cylindrical)) THEN BEGIN

; --------------------------------------------------------------------------
;  Cylindrical (2D) coordinates:    dV = (rp^2-rm^2)/2 * dz
; --------------------------------------------------------------------------

   rp = x1 + 0.5*dx1 
   rm = x1 - 0.5*dx1 
  
   dVr = 0.5*(rp^2 - rm^2)
   dVz = dx2
   dV  = FLTARR(NX1,NX2)

 ENDIF ELSE IF (KEYWORD_SET(spherical)) THEN BEGIN

; --------------------------------------------------------------------------
;  Spherical coordinates:  dV = (rp^3-rm^3)/3 * (cos(thm)-cos(thp)) * dphi
; --------------------------------------------------------------------------

   rp = x1 + 0.5*dx1 & thp = x2 + 0.5*dx2 
   rm = x1 - 0.5*dx1 & thm = x2 - 0.5*dx2
  
   dVr   = (rp^3 - rm^3)/3.0
   dVth  = cos(thm) - cos(thp)
   dVphi = dx3
   dV = FLTARR(NX1,NX2,NX3)
   FOR k=0,NX3-1 DO dV(*,*,k) = (dVr#dVth)*dVphi(k)
 ENDIF ELSE BEGIN

; --------------------------------------------------------------------------
;  Cartesian coordinates:  dV = dx*dy*dz
; --------------------------------------------------------------------------

   dV = FLTARR(NX1,NX2,NX3)
   FOR k=0,NX3-1 DO dV(*,*,k) = (dx1#dx2)*dx3(k)
 ENDELSE
 
; -----------------------------------------------------------------
;   Compute average with respect to all or some of the coordinates.
; -----------------------------------------------------------------

 dV = REFORM(dV)
 IF (NOT KEYWORD_SET(dir)) THEN BEGIN
   RETURN, TOTAL(q*dV,/double)/TOTAL(dV,/double)
 ENDIF ELSE BEGIN
   q = q*dV
   FOR d=ndim,1,-1 DO BEGIN
     IF (dir[d-1] NE 0) THEN BEGIN
       q  = TOTAL(q,d)
       dV = TOTAL(dV,d)
     ENDIF
   ENDFOR
   RETURN, q/dV
 ENDELSE
END

;+
;
; NAME:        PDIFF
;
; AUTHOR:      A. Mignone (mignone@ph.unito.it)
;
; SYNTAX:      dq = PDIFF(q, /X1_DIR, /X2_DIR, /X3_DIR, /PERIODIC)
;
; PURPOSE:     compute the derivative of a 2D or 3D array "q" in the 
;              direction given by either one of /X1_DIR, /X2_DIR or /X3_DIR.
;              The computational grid is defined from the most recent call 
;              to the PLOAD function.
;              On output, "dq" is an array of the same size containing
;              the derivative of the function computed using standard
;              central finite differences. At boundaries, one sided 
;              backward of forward differences are used, unless the /PERIODIC
;              keyword is used.
;
; ARGUMENTS:   
;
;   q          2D or 3D array with (NX1,NX2) or (NX1,NX2,NX3) points, 
;              respectively.
;                      
;   Xn_DIR     either one of X1_DIR, X2_DIR, X3_DIR setting the direction 
;              (x1, x2 or x3) along which the derivative should be
;              taken.
;
; KEYWORDS:
;
;   PERIODIC   assume the domain is periodic in the specified direction. 
;              Central difference is used throughout, including boundaries.
;              
; LAST MODIFIED:   Jan 29, 2012 by A. Mignone (mignone@ph.unito.it)
;
;-
FUNCTION PDIFF, q, X1_DIR=X1_DIR, X2_DIR=X2_DIR, X3_DIR=X3_DIR,$
                periodic=periodic

 COMMON PLUTO_GRID

 dq = FLTARR(NX1,NX2,NX3)
 IF (KEYWORD_SET(X1_DIR)) THEN BEGIN

   IF (KEYWORD_SET(periodic)) THEN BEGIN
     inv_dh = 1.0/(2.0*dx1)
     dq     = SHIFT(q,-1) - SHIFT(q,1)
     FOR i=0,NX1-1 DO dq(i,*,*) = dq(i,*,*)*inv_dh(i)
     RETURN,dq
   ENDIF

   inv_dh = FLTARR(NX1)

   inv_dh(0) = 1.0/(-x1(2) + 4.0*x1(1) - 3*x1(0))
   dq(0,*,*) = (-q(2,*,*) + 4.0*q(1,*,*) - 3.0*q(0,*,*))*inv_dh(0)

   FOR i = 1, NX1-2 DO inv_dh(i) = 1.0/(x1(i+1)-x1(i-1))
   FOR i = 1, NX1-2 DO dq(i,*,*) = (q(i+1,*,*) - q(i-1,*,*))*inv_dh(i)

   inv_dh(i) = 1.0/(x1(i-2) - 4.0*x1(i-1) + 3.0*x1(i))
   dq(i,*,*) = (q(i-2,*,*) - 4.0*q(i-1,*,*) + 3.0*q(i,*,*))*inv_dh(i)
   RETURN,dq
 ENDIF

 IF (KEYWORD_SET(X2_DIR)) THEN BEGIN

   IF (KEYWORD_SET(periodic)) THEN BEGIN
     inv_dh = 1.0/(2.0*dx2)
     dq     = SHIFT(q,0,-1) - SHIFT(q,0,1)
     FOR j=0,NX2-1 DO dq(*,j,*) = dq(*,j,*)*inv_dh(j)
     RETURN,dq
   ENDIF

   inv_dh = FLTARR(NX2)

   inv_dh(0) = 1.0/(-x2(2) + 4.0*x2(1) - 3*x2(0))
   dq(*,0,*) = (-q(*,2,*) + 4.0*q(*,1,*) - 3.0*q(*,0,*))*inv_dh(0)

   FOR j = 1, NX2-2 DO inv_dh(j) = 1.0/(x2(j+1)-x2(j-1))
   FOR j = 1, NX2-2 DO dq(*,j,*) = (q(*,j+1,*) - q(*,j-1,*))*inv_dh(j)

   inv_dh(j) = 1.0/(x2(j-2) - 4.0*x2(j-1) + 3.0*x2(j))
   dq(*,j,*) = (q(*,j-2,*) - 4.0*q(*,j-1,*) + 3.0*q(*,j,*))*inv_dh(j)
   RETURN,dq
 ENDIF

 IF (KEYWORD_SET(X3_DIR)) THEN BEGIN

   IF (KEYWORD_SET(periodic)) THEN BEGIN
     inv_dh = 1.0/(2.0*dx3)
     dq     = SHIFT(q,0,0,-1) - SHIFT(q,0,0,1)
     FOR k=0,NX3-1 DO dq(*,*,k) = dq(*,*,k)*inv_dh(k)
     RETURN,dq
   ENDIF

   inv_dh = FLTARR(NX3)

   inv_dh(0) = 1.0/(-x3(2) + 4.0*x3(1) - 3*x3(0))
   dq(*,*,0) = (-q(*,*,2) + 4.0*q(*,*,1) - 3.0*q(*,*,0))*inv_dh(0)

   FOR k = 1, NX3-2 DO inv_dh(k) = 1.0/(x3(k+1)-x3(k-1))
   FOR k = 1, NX3-2 DO dq(*,*,k) = (q(*,*,k+1) - q(*,*,k-1))*inv_dh(k)

   inv_dh(k) = 1.0/(x3(k-2) - 4.0*x3(k-1) + 3.0*x3(k))
   dq(*,*,k) = (q(*,*,k-2) - 4.0*q(*,*,k-1) + 3.0*q(*,*,k))*inv_dh(k)
   RETURN,dq
 ENDIF

END

;+
;
; NAME:        PCURL
;
; AUTHOR:      A. Mignone (mignone@ph.unito.it)
;
; SYNTAX:      curl_v = PCURL(v1,v1[,v3],/geometry)
;
; PURPOSE:     compute the curl of a vector field {v1,v2} (in 2D) or 
;              {v1,v2,v3} (in 2.5D or 3D) on the computational grid defined
;              during the most recent call to the PLOAD function.
;              On output, "curl_v" contains the three components of the curl
;              stored as 
;
;                 dq(*,*,*,0) -> (curl_v)_x1
;                 dq(*,*,*,1) -> (curl_v)_x2
;                 dq(*,*,*,2) -> (curl_v)_x3
;              
;              Note: if only 2 vector components are supplied, the output
;                    consists in a single 2D array w which, in the case
;                    Cartesian coordinates, is w = dv2/dx1 - dv1/dx2.
; ARGUMENTS:   
;
;   v1,v2[,v3]   2D or 3D arrays with (NX1,NX2) or (NX1,NX2,NX3) points, 
;                respectively, giving the vector components along the 
;                three coordinate directions.
;                       
; KEYWORDS:
;
;   geometry   specify the geometry. 
;              Possible options are /polar or /spherical. 
;              If none is given, Cartesian coordinates are assumed.
;                          
; EXAMPLES:
;
;  * Example #1: compute the vorticity in 2D polar coordinates:
;
;    IDL> drho = PCURL(v1, v2, /polar)
;
;  * Example #2: compute the three components of the electric current in 3D
;                Cartesian coordinates and take the modulus:
;
;    IDL> J = PCURL(B1, B2, B3)
;    IDL> J2 = J(*,*,*,0)^2 + J(*,*,*,1)^2 + J(*,*,*,2)^2
;
; LAST MODIFIED:   Jan 29, 2012 by A. Mignone (mignone@ph.unito.it)
;
;-
FUNCTION PCURL, v1, v2, v3, polar=polar, spherical=spherical 

 COMMON PLUTO_GRID

; --------------------------------------------------------
;   check if the array size is correct
; --------------------------------------------------------

 ndim = SIZE(v1,/N_DIMENSIONS)
 sz   = SIZE(v1,/DIMENSIONS)
 n123 = [NX1,NX2,NX3]
 FOR idim = 0, ndim-1 DO BEGIN
   IF (sz[idim] NE n123[idim]) THEN BEGIN
     PRINT,"! Scalar dimensions does not match with PLUTO grid"
     RETURN,-1
   ENDIF
 ENDFOR

 IF (N_PARAMS() EQ 2) THEN BEGIN

   w = FLTARR(NX1,NX2)

   ; --------------------------------------------------------------
   ;  2D Polar/Spherical coord, w = 1/r*d(r*v2)/dx1 - 1/r*dv1/dx2
   ; --------------------------------------------------------------

   IF (KEYWORD_SET(polar)) THEN BEGIN
     Ar = x1#replicate(1.0,NX2)
     w  = (PDIFF(Ar*v2, /X1_DIR) - PDIFF(v1, /X2_DIR, /periodic))/Ar
     RETURN,w
   ENDIF 

   IF (KEYWORD_SET(spherical)) THEN BEGIN
     Ar = x1#replicate(1.0,NX2)
     w  = (PDIFF(Ar*v2, /X1_DIR) - PDIFF(v1, /X2_DIR))/Ar
     RETURN,w
   ENDIF 

   ; ----------------------------------------------------------------
   ;  2D Cartesian, w = dv2/dx1 - dv1/dx2
   ; ----------------------------------------------------------------

   w  = PDIFF(v2, /X1_DIR) - PDIFF(v1, /X2_DIR)
   RETURN,w

 ENDIF; N_PARAMS() EQ 2

 w = FLTARR(NX1,NX2,NX3,3)

; ----------------------------------------------------------------
;       3D Polar Coordinates,
;
;   w = curl(v1,v2,v3) = [1/r*dv3/dx2 - dv2/dx3, 
;                             dv1/dx3 - dv3/dx1,
;                         1/r*d(r*v2)/dx1 - 1/r*dv1/dx2]
; ----------------------------------------------------------------

 IF (KEYWORD_SET(polar)) THEN BEGIN
   A  = x1#replicate(1.0,NX2)
   Ar = FLTARR(NX1,NX2,NX3)
   FOR k=0,NX3-1 DO Ar(*,*,k) = A
   
   w(*,*,*,0) =   PDIFF(v3, /X2_DIR, /periodic)/Ar 
   w(*,*,*,1) = - PDIFF(v3, /X1_DIR) 
   w(*,*,*,2) =  (PDIFF(Ar*v2, /X1_DIR) - PDIFF(v1, /X2_DIR, /periodic))/Ar 

   IF (ndim EQ 2) THEN RETURN, REFORM(w)

   w(*,*,*,0) -=  PDIFF(v2, /X3_DIR) 
   w(*,*,*,1) +=  PDIFF(v1, /X3_DIR) 

   RETURN,w
 ENDIF

; ----------------------------------------------------------------
;        3D   Spherical Coordinates, s = sin(theta)
;
;   w = curl(v1,v2,v3) = [1/(r*s)*d(s*v3)/dx2 - 1/(r*s)*dv2/dx3, 
;                         1/(r*s)*dv1/dx3     - 1/r*d(r*v3)/dx1,
;                         1/r*d(r*v2)/dx1     - 1/r*dv1/dx2]
; ----------------------------------------------------------------

 IF (KEYWORD_SET(spherical)) THEN BEGIN
   Ar   = FLTARR(NX1,NX2,NX3)
   Ath  = FLTARR(NX1,NX2,NX3)
   Arth = FLTARR(NX1,NX2,NX3)

   A  = x1#replicate(1.0,NX2)      & FOR k=0,NX3-1 DO Ar(*,*,k)   = A 
   A  = replicate(1.0,NX1)#sin(x2) & FOR k=0,NX3-1 DO Ath(*,*,k)  = A 
   A  = x1#sin(x2)                & FOR k=0,NX3-1 DO Arth(*,*,k) = A 
   
   w(*,*,*,0) =   PDIFF(Ath*v3, /X2_DIR)/Arth
   w(*,*,*,1) = - PDIFF( Ar*v3, /X1_DIR)/Ar 
   w(*,*,*,2) =  (PDIFF( Ar*v2, /X1_DIR) - PDIFF(v1, /X2_DIR))/Ar 

   IF (ndim EQ 2) THEN RETURN, REFORM(w)

   w(*,*,*,0) -=  PDIFF(v2, /X3_DIR, /periodic)/Arth
   w(*,*,*,1) +=  PDIFF(v1, /X3_DIR, /periodic)/Arth 

   RETURN,w
 ENDIF

; ----------------------------------------------------------------
;        3D   Cartesian Coordinates (default)
;
;       w = curl(v1,v2,v3) = [dv3/dx2 - dv2/dx3, 
;                             dv1/dx3 - dv3/dx1,
;                             dv2/dx1 - dv1/dx2]
; ----------------------------------------------------------------

 w(*,*,*,0) =   PDIFF(v3, /X2_DIR)
 w(*,*,*,1) = - PDIFF(v3, /X1_DIR) 
 w(*,*,*,2) =   PDIFF(v2, /X1_DIR) - PDIFF(v1, /X2_DIR) 

 IF (ndim EQ 2) THEN RETURN, REFORM(w)

 w(*,*,*,0) -=  PDIFF(v2, /X3_DIR)
 w(*,*,*,1) +=  PDIFF(v1, /X3_DIR)

 RETURN,w
END

;+
;
; NAME:        PGRAD
;
; AUTHOR:      A. Mignone (mignone@ph.unito.it)
;
; SYNTAX:      dq = PGRAD(q,/GEOMETRY)
;
; PURPOSE:     compute the gradient of a scalar function "q" on the 
;              computational grid defined during the most recent call to
;              the PLOAD function.
;              On output "dq" contains the three components of the gradient
;              stored as 
;
;                 dq(*,*,*,0) -> dq/dx1
;                 dq(*,*,*,1) -> dq/dx2
;                 dq(*,*,*,2) -> dq/dx3
;              
; ARGUMENTS:   
;
;   q          a 2D or 3D array with (NX1,NX2) or (NX1,NX2,NX3) points, 
;              respectively. Grid information must already have been stored
;              during the most recent call to PLOAD function.
;
; KEYWORDS:
;
;   GEOMETR   specify the geometry. 
;             Possible options are /POLAR or /SPHERICAL. 
;             If none is given, Cartesian coordinates are assumed.
;                          
; EXAMPLES:
;
;  * Example #1: compute the gradient of rho in spherical coordinates:
;
;    IDL> drho = PGRAD(rho, /SPHERICAL)
;
; LAST MODIFIED:   Sep 25, 2012 by A. Mignone (mignone@ph.unito.it)
;
;-
FUNCTION PGRAD, q, polar=polar, spherical=spherical 

 COMMON PLUTO_GRID

; --------------------------------------------------------
;   check if the array size is correct
; --------------------------------------------------------

 ndim = SIZE(q,/N_DIMENSIONS)
 sz   = SIZE(q,/DIMENSIONS)
 n123 = [NX1,NX2,NX3]
 FOR idim = 0, ndim-1 DO BEGIN
   IF (sz[idim] NE n123[idim]) THEN BEGIN
     PRINT,"! Scalar dimensions does not match with PLUTO grid"
     RETURN,-1
   ENDIF
 ENDFOR

 dq = FLTARR(NX1,NX2,NX3,ndim)

; ----------------------------------------------------------------
;                    Polar coordinates:  
;
;           grad(q) = [dq/dr, 1/r*dq/dphi, dq/dz]
; ----------------------------------------------------------------

 IF (KEYWORD_SET(polar)) THEN BEGIN
   Ar = FLTARR(NX1,NX2,NX3)
   A  = x1#replicate(1.0,NX2)    
   FOR k=0,NX3-1 DO Ar(*,*,k) = A 

   dq = FLTARR(NX1,NX2,NX3,ndim)
   dq(*,*,*,0) = PDIFF(q, /X1_DIR)
   dq(*,*,*,1) = PDIFF(q, /X2_DIR)/Ar
   IF (ndim EQ 2) THEN RETURN, REFORM(dq)
   
   dq(*,*,*,2) = PDIFF(q, /X3_DIR)
   RETURN,dq
 ENDIF

; ----------------------------------------------------------------
;                      Spherical coordinates:  
;
;     grad(q) = [dq/dr, 1/r*dq/dtheta, 1/(r*sin(theta))*dq/dphi]
; ----------------------------------------------------------------

 IF (KEYWORD_SET(spherical)) THEN BEGIN
   Ar = FLTARR(NX1,NX2,NX3)
   A  = x1#replicate(1.0,NX2)      
   FOR k=0,NX3-1 DO Ar(*,*,k) = A 

   dq(*,*,*,0) = PDIFF(q, /X1_DIR)
   dq(*,*,*,1) = PDIFF(q, /X2_DIR)/Ar
   IF (ndim EQ 2) THEN RETURN, REFORM(dq)

   Arth = FLTARR(NX1,NX2,NX3)
   A    = x1#sin(x2)  
   FOR k=0,NX3-1 DO Arth(*,*,k) = A 
   dq(*,*,*,2) = PDIFF(q, /X3_DIR)/Arth
   RETURN, dq
 ENDIF

; ----------------------------------------------------------------
;           Cartesian Coordinates (default)
;
;           grad(q) = [dq/dx, dq/dy, dq/dz]
; ----------------------------------------------------------------

 dq = FLTARR(NX1,NX2,NX3,ndim)
 dq(*,*,*,0) = PDIFF(q, /X1_DIR)
 dq(*,*,*,1) = PDIFF(q, /X2_DIR)
 IF (ndim EQ 2) THEN RETURN, REFORM(dq)

 dq(*,*,*,2) = PDIFF(q, /X3_DIR)
 RETURN, dq

END
