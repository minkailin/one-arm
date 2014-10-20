#include"param.h"

/* ###################################################################
 
   FILE: limiters.c

   PURPOSE: Compute slope limiters needed in the interpolation 
            routines; slopes are returned as

                           dV|
            dV_LIM[i]  =   --|  Delta x_i 
                           dx|i 
  
            where dV/dx depends on the choice of the limiter and the 
            grid type (uniform or not).


   ARGUMENTS:

     dwp (in)    : array w/forward  undivided differences 
     dwm (in)    : array w/backward undivided differences
     dw_lim(out) : array of limited slopes
     ibeg (in)     : starting point
     iend (in)     : ending point
     GG (in)     : grid structure in the current sweep direction 


   NOTE:   In developing new slopes limiters, keep in mind that 
           in order to satisfy the TVD condition you must have
           (see TORO, pg 477):

        min(U[i], U[i+1]) <= UR = U[i] + dV_LIM[i] <= max (U[i], U[i+1])

        min(U[i], U[i-1]) <= UL = U[i] - dV_LIM[i] <= max (U[i], U[i-1])
 
   ################################################################### */


/* ###################################################################### */
void flat_lim (double *dwp, double *dwm, double *dw_lim,
               int ibeg, int iend, struct GRID *GG)
/* 
 #
 #   Zero-slope; first order interpolation
 #
 ######################################################################## */
{
  int i;

  for (i = ibeg; i <= iend; i++){
    dw_lim[i] = 0.0;
  }
}

/* ###################################################################### */
void minmod_lim (double *dwp, double *dwm, double *dw_lim,
                 int ibeg, int iend, struct GRID *GG)
/* 
 #
 # Evaluate slopes using the minmod limiter:
 #
 #
 #                   xm  when |xm| < |xp|  AND xp*xm > 0
 #  minmod(xp, xm) = xp  when |xp| < |xm|  AND xp*xm > 0
 #                   0   otherwise
 #
 #  Here xp and xm are approximations to forward and backward first 
 #  derivative.
 # 
 ######################################################################## */
{
  int i;
  double xp, xm;

  for (i = ibeg; i <= iend; i++){
    xp = dwp[i]*GG->c1[i];  /* forward derivative  * dx/2  */
    xm = dwm[i]*GG->c2[i];  /* backward derivative * dx/2 */
    if (xp*xm > 0.0){
      dw_lim[i] = fabs(xp) < fabs(xm) ? xp : xm;
    }else{
      dw_lim[i] = 0.0;
    }
  }
}

/* ###################################################################### */
void vanleer_lim (double *dwp, double *dwm, double *dw_lim,
                  int ibeg, int iend, struct GRID *GG)
/* 
 #
 #    vl = 2*x*y/(x + y)   when x*y > 0
 #    vl = 0               otherwise
 #
 ######################################################################## */
{
  int i;
  double scrh, s;

  if (GG->uniform) {     /*  uniform grids  */

    for (i = ibeg; i <= iend; i++){
      s         = dwp[i]*dwm[i];
      dw_lim[i] = s > 0.0 ? 2.0*s/(dwp[i] + dwm[i]):0.0;        
    }

  } else {    /*    non-uniform grids  */

    for (i = ibeg; i <= iend; i++){
      s  = dwp[i]*dwm[i];
      if (s > 0.0) {
        if (GG->c1[i]*dwp[i] < GG->c2[i]*dwm[i]){
          scrh      = (2.0 - GG->c1[i])/GG->c2[i];
          dw_lim[i] = 2.0*s/(scrh*dwp[i] + dwm[i]);
        }else{
          scrh      = (2.0 - GG->c2[i])/GG->c1[i];
          dw_lim[i] = 2.0*s/(dwp[i] + dwm[i]*scrh);
        }
      }else{
        dw_lim[i] = 0.0;
      }
    }

  }
}

/* ###################################################################### */
void woodward_lim (double *dwp, double *dwm, double *dw_lim,
                   int ibeg, int iend, struct GRID *GG)
/* 
 #
 #
 # Evaluate slopes using the woodward limiter:
 #
 #      woodward_lim (xp, xm, xc) = MinMod [2*xp , 2*xm , xc]
 #
 #  Where xp, xm and xc are, respectively, approximations to 
 #  the forward, backward and central derivatives.
 #
 ######################################################################## */
{
  int i;
  double xc;

  if (GG->uniform){

    for (i = ibeg; i <= iend; i++){
      if ( dwp[i]*dwm[i] > 0.0 ){
        xc        = 0.5*(dwp[i] + dwm[i]);
        dw_lim[i] = 2.0*(fabs(dwp[i]) < fabs(dwm[i]) ? dwp[i]:dwm[i]);
        dw_lim[i] = fabs(dw_lim[i]) < fabs(xc) ? dw_lim[i]: xc;
      }else{
        dw_lim[i] = 0.0;
      }
    }
  
  } else {     /*    non-uniform grids    */

    for (i = ibeg; i <= iend; i++){
      if (dwp[i]*dwm[i] > 0.0 ){
        xc        = dwp[i]*GG->c3[i] + dwm[i]*GG->c4[i];
        dw_lim[i] = 2.0*(fabs(dwp[i]) < fabs(dwm[i]) ? dwp[i]:dwm[i]);
        dw_lim[i] = fabs(dw_lim[i]) < fabs(xc) ? dw_lim[i]:xc;
      }else{
        dw_lim[i] = 0.0;
      }
    }
  }
}

/* ###################################################################### */
void gminmod_lim (double *dwp, double *dwm, double *dw_lim,
                  int ibeg, int iend, struct GRID *GG)
/* 
 #
 #
 # Evaluate slopes using the general mindmod limiter:
 #
 #      gminmod (xp, xm, xc) = MinMod [a*xp , a*xm , xc]
 #
 #  Where xp, xm and xc are, respectively, approximations to 
 #  the forward, backward and central derivatives, and 
 #  1 < a < 2.
 #
 # setting a = 1  yields  the minmod limiter
 # setting a = 2  yields  the woodward limiter
 #
 #
 ######################################################################## */
{
  int i;
  double xc, a = 1.5;

  if (GG->uniform){

    for (i = ibeg; i <= iend; i++){
      if ( dwp[i]*dwm[i] > 0.0 ){
        xc        = 0.5*(dwp[i] + dwm[i]);
        dw_lim[i] = a*(fabs(dwp[i]) < fabs(dwm[i]) ? dwp[i]:dwm[i]);
        dw_lim[i] = fabs(dw_lim[i]) < fabs(xc) ? dw_lim[i]: xc;
      }else{
        dw_lim[i] = 0.0;
      }
    }
  
  } else {     /*    non-uniform grids    */

    for (i = ibeg; i <= iend; i++){
      if (dwp[i]*dwm[i] > 0.0 ){
        xc        = dwp[i]*GG->c3[i] + dwm[i]*GG->c4[i];
        dw_lim[i] = a*(fabs(dwp[i]) < fabs(dwm[i]) ? dwp[i]:dwm[i]);
        dw_lim[i] = fabs(dw_lim[i]) < fabs(xc) ? dw_lim[i]:xc;
      }else{
        dw_lim[i] = 0.0;
      }
    }
  }
}

#define EPS_VA  1.e-18
/* ###################################################################### */
void vanalbada_lim (double *dwp, double *dwm, double *dw_lim,
                    int ibeg, int iend, struct GRID *GG)
/* 
 #
 #
 #   va(x,y) = x*y*(x + y)/(x*x + y*y)
 #
 #
 # Notice although this limiter does not require the switch x*y > 0,
 # we nevertheless employ it.
 # 
 ######################################################################## */
{
  int  i;
  double s, xp, xm;
  double dwp2, dwm2, xp2, xm2;

  if (GG->uniform) {     /*  uniform grids  */
  
    for (i = ibeg; i <= iend; i++){
      if ( (s = dwp[i]*dwm[i]) > 0.0) {    
        dwp2 = dwp[i]*dwp[i];
        dwm2 = dwm[i]*dwm[i];
        dw_lim[i] =   (dwp[i]*(dwm2 + EPS_VA) + dwm[i]*(dwp2 + EPS_VA))
                    / (dwp2 + dwm2 + EPS_VA);
      } else {
        dw_lim[i] = 0.0;
      }
    }
      
  }else{   /*   non-uniform grid    */

    for (i = ibeg; i <= iend; i++){
      xp = dwp[i]*GG->c1[i];
      xm = dwm[i]*GG->c2[i];
      if ( (s = xp*xm) > 0.0) {
        xp2 = xp*xp;
        xm2 = xm*xm;
        dw_lim[i] =   (xp*(xm2 + EPS_VA) + xm*(xp2 + EPS_VA)) 
                    / (xp2 + xm2 + EPS_VA);
      } else {
        dw_lim[i] = 0.0;
      }
    }                            
  }
}
#undef EPS_VA
/* ###################################################################### */
void umist_lim (double *dwp, double *dwm, double *dw_lim,
                int ibeg, int iend, struct GRID *GG)
/* 
 #
 #
 # Evaluate slopes using the Umist Limiter :
 #
 #   UMIST (x,y) = MinMod [2*xp , 2*xm , 1/4 (xp + 3xm) , 1/4 (xm + 3xp)]
 #
 #
 ######################################################################## */
{
  int  i;
  double x1,x2, xp, xm;

  if (GG->uniform) {
  
    for (i = ibeg; i <= iend; i++){
      if ( dwp[i]*dwm[i] > 0.0) {    
        x1 = 0.25*(dwp[i] + 3.0*dwm[i]);
        x2 = 0.25*(dwm[i] + 3.0*dwp[i]);                                   
        dw_lim[i] = 2.0*(fabs(dwp[i]) < fabs(dwm[i]) ? dwp[i]:dwm[i]);
        dw_lim[i] = fabs(dw_lim[i]) < fabs(x1) ? dw_lim[i]:x1;
        dw_lim[i] = fabs(dw_lim[i]) < fabs(x2) ? dw_lim[i]:x2;
      } else {
        dw_lim[i] = 0.0;
      }
    }

  }else{
          
    for (i = ibeg; i <= iend; i++){
      if ( dwp[i]*dwm[i] > 0.0) {
        xp = dwp[i]*GG->c1[i];
        xm = dwm[i]*GG->c2[i];
        x1 = 0.25*(xp + 3.0*xm);
        x2 = 0.25*(xm + 3.0*xp);
        dw_lim[i] = 2.0*(fabs(dwp[i]) < fabs(dwm[i]) ? dwp[i]:dwm[i]);
        dw_lim[i] = fabs(dw_lim[i]) < fabs(x1) ? dw_lim[i]: x1;
        dw_lim[i] = fabs(dw_lim[i]) < fabs(x2) ? dw_lim[i]: x2;
      } else {
        dw_lim[i] = 0.0;
      }
    }

  }
}

/* ================================================ */
/* double weightmm_lim (x,y,dh)

    Weighted Minmod limiter :

    WGHT_MM(x,y) = MinMod [ x , omega * y]

    where omega must be

    1 <= omega <= (eta-3)/(eta-1)

    Weighted minmod can also be a TVB minmod limiter by
    setting MTVB to a positive integer value; if done this
    limiter becomes :

    WGHT_MM(x,y) = MinMod [ x , omega * y + M*dh^2]

    where dh is the actual GRID point distance.


{
      int nvar ;
      parameter (nvar = 2);
      double   x,y,u(nvar),sx,dh;
      double   xlim;

      sx = 0.d0;
      if (mtvb.ne.0) {
         sx = dh*dh;
         sx = dble(mtvb)*sx;
         sx = DSIGN (sx,x);
      }

      u(1) = x;
      u(2) = y*omega + sx ;
      MINMOD (u,nvar,xlim);

      weightmm_lim = xlim;

      return ??;
} */


/* ###################################################################### */
void colella_lim (double *dwp, double *dwm, double *dw_lim,
                  int ibeg, int iend, struct GRID *GG)
/* 
 #
 #
 #   Evaluate slopes using Colella's fourth-order slope
 #
 #
 ######################################################################## */
{
  int i;
  static double *s;
  static double *dqf, *dqc, *dqlim;
  double scrh;

  if (s == NULL){
    s     = vector(NMAX_POINT);
    dqf   = vector(NMAX_POINT);
    dqc   = vector(NMAX_POINT);
    dqlim = vector(NMAX_POINT);
  }

  if (GG->uniform){

    for (i = ibeg-1; i <= iend+1; i++){
      if ( dwp[i]*dwm[i] > 0.0 ){
        dqc[i]   = 0.5*(dwp[i] + dwm[i]);
        s[i]     = DSIGN(dqc[i]);
        dqlim[i] = 2.0*MIN(fabs(dwp[i]), fabs(dwm[i]));        
      }else{
        dqlim[i] = s[i] = 0.0;
      }
      dqf[i]   = MIN(fabs(dqc[i]), dqlim[i]) * s[i];
    }

    for (i = ibeg; i <= iend; i++){
      scrh = 4./3.*dqc[i] - 1./6.*(dqf[i+1] + dqf[i-1]);
      dw_lim[i] = MIN(fabs(scrh), dqlim[i])*s[i];        
    }
  
  } else {     /*    non-uniform grids    */

    printf ("Sorry, fourth order slopes are not yet implemented for\n");
    printf ("non uniform grids \n");
    exit(1);

  }

}
/* ###################################################################### */
void colella_lim_2 (double *dwp, double *dwm, double *dw_lim,
                  int ibeg, int iend, struct GRID *GG)
/* 
 #
 #
 #   Evaluate slopes using Colella's fourth-order slope
 #
 #  Ref:  Miller, G.H and P. COlella, 
 #        "A high-order Eulerian Godunov Method for 
 #         Elastic-Plastic Flow in Solids", JCP 167,131-176 (2001)
 #    
 #                             +
 #
 #        Saltzman, J, " An Unsplit 3D Upwind Method for 
 #                       Hyperbolic Conservation Laws", 
 #                       JCP 115, 153-168 (1994)
 #
 ######################################################################## */
{
  int i;
  static double *s;
  static double *dqf, *dqc, *dqlim;
  double scrh;

  if (s == NULL){
    s     = vector(NMAX_POINT);
    dqf   = vector(NMAX_POINT);
    dqc   = vector(NMAX_POINT);
    dqlim = vector(NMAX_POINT);
  }

  for (i = ibeg-1; i <= iend+1; i++){
    dqc[i] = dwp[i] + dwm[i];
    if ( dwp[i]*dwm[i] > 0.0 ){
      s[i]     = DSIGN(dqc[i]);
      dqlim[i] = MIN(fabs(dwp[i]), fabs(dwm[i]));        
    }else{
      dqlim[i] = s[i] = 0.0;
    }
    scrh   = GG->dxv[i]/(GG->dxv[i-1] + 2.0*GG->dxv[i] + GG->dxv[i+1]);
    dqf[i] = 2.0*MIN(scrh*fabs(dqc[i]), dqlim[i])*s[i];
  }
  for (i = ibeg; i <= iend; i++){
    scrh = GG->dxv[i]/(0.25*GG->dxv[i-1] + GG->dxv[i] + 0.25*GG->dxv[i+1]);     
    scrh = scrh*(dqc[i] - 0.25*(dqf[i+1] + dqf[i-1]));

    dw_lim[i] = MIN(fabs(scrh), 2.0*dqlim[i])*s[i];        
  }  

}
