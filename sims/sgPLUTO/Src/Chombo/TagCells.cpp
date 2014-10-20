#include <cstdio>
#include <string>
using std::string;

#include "PatchPluto.H"
#include "LoHiSide.H"

static void computeRefVar(double ***UU[], double ***q, double, RBox *Ubox);

#if (EOS != ISOTHERMAL) && (AMR_EN_SWITCH == NO)
 #define REF_VAR  ENG
#else
 #define REF_VAR  RHO
#endif

// #define REF_VAR -1  /* means user-defined */

#define REF_CRIT 2   /* 1 == first derivative, 2 == second derivative */

/* ************************************************************************* */
void PatchPluto::computeRefGradient(FArrayBox& gFab, FArrayBox& UFab, const Box& b)
/*!
 * Compute numerical gradient of the solution in order to tag 
 * zones for refinement.
 * The gradient is computed by standard finite differences using
 *
 * - REF_CRIT equal to 1 --> compute (normalized) gradient using 1st 
 *                           derivative of the solution;
 * - REF_CRIT equal to 2 --> compute (normalized) gradient using 2nd 
 *                           derivative of the solution (default);
 *                           This approach is based on Lohner (1987).
 *
 * Zones will be flagged for refinement whenever grad[k][j][i] exceeds 
 * the threshold value specified by the 'Refine_thresh' parameter read in
 * pluto.ini.
 *
 * Derivatives are computed using the conserved variable U[REF_VAR] 
 * where REF_VAR is taken to be energy density (default).
 * However, by setting REF_VAR = -1, you can provide your own 
 * physical variable through the function computeRefVar().
 * 
 * \authors C. Zanni   (zanni@oato.inaf.it)\n
 *          A. Mignone (mignone@ph.unito.it)
 * \date    Oct 11, 2012
 *************************************************************************** */
{
  CH_assert(m_isDefined);

  int nv, i, j, k;

  double rp, rm, r, tp, tm;
  double x, dqx_p, dqx_m, dqx, d2qx, den_x;
  double y, dqy_p, dqy_m, dqy, d2qy, den_y;
  double z, dqz_p, dqz_m, dqz, d2qz, den_z;

  double gr1, gr2, eps = 0.01;
  double ***UU[NVAR], ***q, ***grad;
  RBox  Ubox, Gbox;
  
  rp = rm = r = 1.0;
  tp = tm = 1.0;

/* -- check ref criterion -- */

  #if REF_CRIT != 1 && REF_CRIT != 2
   print ("! TagCells.cpp: Refinement criterion not valid\n");
   QUIT_PLUTO(1);
  #endif

/* -----------------------------------------------
    The solution array U is defined on the box 
    [Uib, Uie] x [Ujb, Uje] x [Ukb, Uke], which 
    differs from that of gFab ([Gib,...Gke]), 
    typically one point larger in each direction. 
   ----------------------------------------------- */
    
  Ubox.jb = Ubox.je = Ubox.kb = Ubox.ke = 0;
  Gbox.jb = Gbox.je = Gbox.kb = Gbox.ke = 0;

  D_EXPAND(Ubox.ib = UFab.loVect()[IDIR]; Ubox.ie = UFab.hiVect()[IDIR]; ,
           Ubox.jb = UFab.loVect()[JDIR]; Ubox.je = UFab.hiVect()[JDIR]; ,
           Ubox.kb = UFab.loVect()[KDIR]; Ubox.ke = UFab.hiVect()[KDIR]; );

  D_EXPAND(Gbox.ib = gFab.loVect()[IDIR]; Gbox.ie = gFab.hiVect()[IDIR]; ,
           Gbox.jb = gFab.loVect()[JDIR]; Gbox.je = gFab.hiVect()[JDIR]; ,
           Gbox.kb = gFab.loVect()[KDIR]; Gbox.ke = gFab.hiVect()[KDIR]; );

  for (nv = 0; nv < NVAR; nv++){
    UU[nv] = ArrayBoxMap(Ubox.kb, Ubox.ke, 
                         Ubox.jb, Ubox.je, 
                         Ubox.ib, Ubox.ie, UFab.dataPtr(nv));
  }
  grad = ArrayBoxMap(Gbox.kb, Gbox.ke, 
                     Gbox.jb, Gbox.je, 
                     Gbox.ib, Gbox.ie, gFab.dataPtr(0));

/* ---------------------------------------------
     set refinement variable
   --------------------------------------------- */

  #if REF_VAR >= 0
   q = UU[REF_VAR];
  #else
   q = ArrayBox(Ubox.kb, Ubox.ke, Ubox.jb, Ubox.je, Ubox.ib, Ubox.ie);
   computeRefVar(UU, q, m_dx, &Ubox);
  #endif

/* ----------------------------------------------------------------
    Main spatial loop for zone tagging based on 1st (REF_CRIT = 1) 
    or 2nd (REF_CRIT = 2) derivative error norm. 
   ---------------------------------------------------------------- */

  BOX_LOOP(&Gbox, k, j, i){
    z = (k + 0.5)*m_dx + g_domBeg[KDIR];
    y = (j + 0.5)*m_dx + g_domBeg[JDIR];
    x = (i + 0.5)*m_dx + g_domBeg[IDIR];

    #if GEOMETRY == CYLINDRICAL
     Real xll = g_domBeg[0] + double(i-1)*m_dx;
     Real xl  = xll + m_dx;
     Real xr  = xl  + m_dx;
     Real xrr = xr  + m_dx;
     rm = xr*xr-xl*xl;
     rp = rm/(xrr*xrr-xr*xr);
     rm = rm/(xl*xl-xll*xll);
/*
     rp = (i+0.5)/(i+1.5);
     rm = (i+0.5)/(i-0.5);
*/
    #endif

    D_EXPAND(dqx_p =    q[k][j][i+1]*rp - q[k][j][i];
             dqx_m = - (q[k][j][i-1]*rm - q[k][j][i]);  ,
             dqy_p =    q[k][j+1][i] - q[k][j][i];
             dqy_m = - (q[k][j-1][i] - q[k][j][i]);     ,
             dqz_p =    q[k+1][j][i] - q[k][j][i];
             dqz_m = - (q[k-1][j][i] - q[k][j][i]);)

  /* --------------------------------------------------------------
      Physical boundary values are not up to date and should be 
      excluded from gradient computation. 
      In this case, left and right derivatives are set equal to 
      each other. This will not trigger refinement in the leftmost 
      and rightmost internal zones (using 2nd derivative) but we 
      really don't care since buffer size will do the job.
     -------------------------------------------------------------- */
      
    D_EXPAND(if (i == 0) dqx_m = dqx_p;  ,
             if (j == 0) dqy_m = dqy_p;  ,
             if (k == 0) dqz_m = dqz_p;)

    D_EXPAND(if (i == m_domain.size(IDIR)-1) dqx_p = dqx_m;  ,
             if (j == m_domain.size(JDIR)-1) dqy_p = dqy_m;  ,
             if (k == m_domain.size(KDIR)-1) dqz_p = dqz_m;)

  /* -----------------------------------------------
         Compute gradient using 1st derivative 
      ---------------------------------------------- */

    #if REF_CRIT == 1
     D_EXPAND(dqx = dqx_p + dqx_m;  ,
              dqy = dqy_p + dqy_m;  ,
              dqz = dqz_p + dqz_m;)

     D_EXPAND(den_x = fabs(q[k][j][i+1]*rp) + fabs(q[k][j][i-1]*rm);  ,
              den_y = fabs(q[k][j+1][i])    + fabs(q[k][j-1][i]);     ,
              den_z = fabs(q[k+1][j][i])    + fabs(q[k-1][j][i]);)

     gr1  = D_EXPAND(dqx*dqx, + dqy*dqy, + dqz*dqz);
     gr1 /= D_EXPAND(den_x*den_x, + den_y*den_y, + den_z*den_z);

     grad[k][j][i] = sqrt(gr1);
    #endif

  /* -----------------------------------------------
         Compute gradient using 2nd derivative 
      ---------------------------------------------- */

    #if REF_CRIT == 2
     D_EXPAND(d2qx = dqx_p - dqx_m;  ,
              d2qy = dqy_p - dqy_m;  ,
              d2qz = dqz_p - dqz_m;)

     D_EXPAND(
       den_x = 2.0*fabs(q[k][j][i]) + fabs(q[k][j][i+1]*rp) + fabs(q[k][j][i-1]*rm);
       den_x = fabs(dqx_p) + fabs(dqx_m) + eps*den_x;    ,

       den_y = 2.0*fabs(q[k][j][i]) + fabs(q[k][j+1][i]) + fabs(q[k][j-1][i]);
       den_y = fabs(dqy_p) + fabs(dqy_m) + eps*den_y;    ,

       den_z = 2.0*fabs(q[k][j][i]) + fabs(q[k+1][j][i]) + fabs(q[k-1][j][i]);
       den_z = fabs(dqz_p) + fabs(dqz_m) + eps*den_z;
     )

     gr2  = D_EXPAND(d2qx*d2qx, + d2qy*d2qy, + d2qz*d2qz);
     gr2 /= D_EXPAND(den_x*den_x, + den_y*den_y, + den_z*den_z);

     grad[k][j][i] = sqrt(gr2);
    #endif
  }

/* --------------------------------------------------------------
    Ok, grad[] has been computed. Now Free memory.
   -------------------------------------------------------------- */
   
  FreeArrayBoxMap(grad, Gbox.kb, Gbox.ke, Gbox.jb, Gbox.je, Gbox.ib, Gbox.ie);
  for (nv = 0; nv < NVAR; nv++){
    FreeArrayBoxMap(UU[nv], Ubox.kb, Ubox.ke, Ubox.jb, Ubox.je, Ubox.ib, Ubox.ie);
  }
  #if REF_VAR == -1
   FreeArrayBox(q, Ubox.kb, Ubox.jb, Ubox.ib);
  #endif
}

/* ********************************************************************* */
void computeRefVar(double ***UU[], double ***q, double dx, RBox *Ubox)
/*!
 * Compute a user-defined array q(U) function of the conserved
 * variables.
 *
 *
 *********************************************************************** */
{
  int nv, i, j, k;
  double x1, x2, x3, r;
  double us[NVAR], vs[NVAR];
  double pm = 0.0, Kin; 

  r = 1.0;

  BOX_LOOP(Ubox, k, j, i) {
    #if GEOMETRY == CYLINDRICAL 
     r = g_domBeg[IDIR] + dx*(i+0.5); /* we could use m_dx if we were inside
                                         a class member */
                                          
    #endif
    for (nv = 0; nv < NVAR; nv++) us[nv] = UU[nv][k][j][i]/r; 

    /* -- recover pressure from conservative vars (MHD only) -- */

    #if AMR_EN_SWITCH == NO
     Kin  = EXPAND(us[MX1]*us[MX1], + us[MX2]*us[MX2], + us[MX3]*us[MX3]); 
     Kin *= 0.5/us[RHO];

     #if PHYSICS == MHD       
      pm = EXPAND(us[BX1]*us[BX1], + us[BX2]*us[BX2], + us[BX3]*us[BX3]); 
      pm = 0.5*pm;
     #endif
     q[k][j][i] = (us[ENG] - Kin - pm)*(g_gamma - 1.0);
    #else
     pm = us[RHO];
     q[k][j][i] = us[ENG]*pow(pm, g_gamma - 1.0);
    #endif

     q[k][j][i] *= r;
 
    #if PHYSICS == RMHD || PHYSICS == RHD
     print ("! TagCells.cpp: refinement valid only for classical hydro\n");
     QUIT_PLUTO(1);
    #endif
  }
}
#undef REF_VAR
#undef REF_CRIT
