#include "pluto.h"

/* *********************************************************************** */
void CheckPrimStates(double **vM, double **vP, double **v0,  int beg, int end)
/*
 *
 *  PURPOSE
 *
 *   check if primitive states vL and vR are physically
 *   admissible. 
 *   Replace their value with v0 otherwise.
 *
 *
 ************************************************************************** */
{
  int    i, nv, switch_to_1st;
  double scrhm, scrhp;
  double *ac, *ap, *am;

  for (i = beg; i <= end; i++){
  
    switch_to_1st = 0;

    ac = v0[i];
    am = vM[i];
    ap = vP[i];

  /*  ----  Prevent unphysical states by revertin to first
            order in time and space,  i.e.  set dw = 0      ----  */

    #if EOS != ISOTHERMAL && EOS != BAROTROPIC
     switch_to_1st = (ap[PRS] < 0.0) || (am[PRS] < 0.0) ;
    #endif
    switch_to_1st = switch_to_1st || 
                    (ap[RHO] < 0.0) || (am[RHO] < 0.0) ;

    /*  ----  Check for superluminal velocities  ---- */

    #if (PHYSICS == RHD && USE_FOUR_VELOCITY == NO) || PHYSICS == RMHD 
     scrhm = EXPAND(am[VX1]*am[VX1], + am[VX2]*am[VX2], + am[VX3]*am[VX3]);
     scrhp = EXPAND(ap[VX1]*ap[VX1], + ap[VX2]*ap[VX2], + ap[VX3]*ap[VX3]);
     switch_to_1st = switch_to_1st || (scrhm >= 1.0);
     switch_to_1st = switch_to_1st || (scrhp >= 1.0);
    #endif

    if (switch_to_1st){
/*
      WARNING (
        print (" ! CheckPrimStates: Unphysical state, ");
        Where (i,NULL);
      )
*/              
      #ifdef STAGGERED_MHD 
       scrhp = ap[BXn];
       scrhm = am[BXn];
      #endif

      for (nv = 0; nv < NVAR; nv++){
        am[nv] = ap[nv] = ac[nv];
      }
 
      #ifdef STAGGERED_MHD
       ap[BXn] = scrhp;
       am[BXn] = scrhm;
      #endif
      
    }
  }
}


/*
#if GEOMETRY == CYLINDRICAL
if (NSWEEP == 1) {
 for (nv = 0; nv < NVAR; nv++){
   if (nv == VX) scrh = fabs(v1[3][nv] + v1[4][nv]);
   else          scrh = fabs(v1[3][nv] - v1[4][nv]);
   
   #if PHYSICS == RMHD 
    if (nv == BZ) scrh = fabs(v1[3][nv] + v1[4][nv]);
   #endif
   
   if (scrh > 1.e-8){
     printf ("symmetry violated, z : %d,  var: %d\n", *nyp, nv);
     SHOW(rhs,4);
     SHOW(rhs,3);
     printf (" --- centered:\n");
     SHOW(vv,0); SHOW(vv,1); SHOW(vv,2); SHOW(vv,3); 
     SHOW(vv,4); SHOW(vv,5); SHOW(vv,6); SHOW(vv,7);
     printf (" --- edge\n");
     SHOW(vl,3); SHOW(vr,3);
     SHOW(vr,2); SHOW(vl,4);

     printf ("Source: \n");
     SHOW(src,3);SHOW(src,4);
     scrhp =  fl[4][MXn]*GG->A[4]/GG->dV[4] - src[4][MXn];
     scrhm = -fr[2][MXn]*GG->A[2]/GG->dV[3] - src[3][MXn];
     printf ("%12.6e  %12.6e\n",scrhp, scrhm);
     exit(1);
   }
 }
}
#endif
*/


/* ************************************************************************* */
void CheckConsStates(double **uM, double **uP, double **u0, int beg, int end)
/*
 *
 *  PURPOSE
 *
 *   check if conservative states uR and uL are
 *   physically admissible. 
 *   Replace their left and right interpolated values 
 *   inside the cell with uu0 otherwise.
 *
 *
 *************************************************************************** */
{
  int   i, nv;
  real  scrhp = 1.0, scrhm = 1.0;
  real  *q0, *qm, *qp;

  for (i = beg; i <= end; i++){
  
    q0 = u0[i];
    qm = uM[i];
    qp = uP[i];

    #if PHYSICS == HD && EOS != ISOTHERMAL
     
     scrhp = EXPAND(qp[MX1]*qp[MX1], + qp[MX2]*qp[MX2], + qp[MX3]*qp[MX3]);
     scrhp = qp[ENG] - 0.5*scrhp/qp[RHO];

     scrhm = EXPAND(qm[MX1]*qm[MX1], + qm[MX2]*qm[MX2], + qm[MX3]*qm[MX3]);
     scrhm = qm[ENG] - 0.5*scrhm/qm[RHO];

    #elif PHYSICS == MHD && EOS != ISOTHERMAL && EOS != BAROTROPIC

     scrhp  = EXPAND(qp[MX1]*qp[MX1], + qp[MX2]*qp[MX2], + qp[MX3]*qp[MX3]);
     scrhp /= qp[RHO];
     scrhp += EXPAND(qp[BX1]*qp[BX1], + qp[BX2]*qp[BX2], + qp[BX3]*qp[BX3]);
     scrhp = qp[ENG] - 0.5*scrhp;

     scrhm  = EXPAND(qm[MX1]*qm[MX1], + qm[MX2]*qm[MX2], + qm[MX3]*qm[MX3]);
     scrhm /= qm[RHO];
     scrhm += EXPAND(qm[BX1]*qm[BX1], + qm[BX2]*qm[BX2], + qm[BX3]*qm[BX3]);
     scrhm  = qm[ENG] - 0.5*scrhm;

    #elif PHYSICS == RHD

     scrhp = EXPAND(qp[MX1]*qp[MX1], + qp[MX2]*qp[MX2], + qp[MX3]*qp[MX3]);
     scrhp = qp[ENG]*qp[ENG] - scrhp - qp[RHO]*qp[RHO]; 

     scrhm = EXPAND(qm[MX1]*qm[MX1], + qm[MX2]*qm[MX2], + qm[MX3]*qm[MX3]);
     scrhm = qm[ENG]*qm[ENG] - scrhm - qm[RHO]*qm[RHO]; 

    #endif

    if (scrhp  < 0.0 || scrhm  < 0.0 ||
        #if EOS != ISOTHERMAL && EOS != BAROTROPIC
         qm[ENG] < 0.0 || qp[ENG] < 0.0 ||
        #endif
        qm[RHO] < 0.0 || qp[RHO] < 0.0) {
      WARNING (
        print (" ! CheckConstStates: Unphysical state, ");
        Where (i,NULL);
      )
      for (nv = 0; nv < NVAR; nv++) {
         qm[nv] = qp[nv] = q0[nv];
      }
    }
  }
}
