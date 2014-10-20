#include "pluto.h"

/* ***************************************************************** */
void Flux (double **u, double **v, double *a2, double **fx, 
           double *p, int beg, int end)
/*
 *
 *
 *
 ******************************************************************* */
{
  int    nv, i;
  double vn;

  for (i = beg ; i <= end; i++) {

    #if USE_FOUR_VELOCITY == YES
     vn = EXPAND( v[i][VXn]*v[i][VXn], 
                + v[i][VXt]*v[i][VXt],
                + v[i][VXb]*v[i][VXb]);

     vn = v[i][VXn]/sqrt(1.0 + vn);
    #else  
     vn = v[i][VXn];
    #endif

    fx[i][RHO]  = u[i][RHO]*vn;
    EXPAND(fx[i][MX1] = u[i][MX1]*vn;  ,
           fx[i][MX2] = u[i][MX2]*vn;  ,
           fx[i][MX3] = u[i][MX3]*vn;)
    fx[i][ENG] = u[i][MXn];
    p[i] = v[i][PRS];
  }
}



