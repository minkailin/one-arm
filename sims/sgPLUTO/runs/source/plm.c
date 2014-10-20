/*
 *    Example taken from the GNU Scientific Library Reference Manual
 *        Edition 1.1, for GSL Version 1.1
 *            9 January 2002
 *               URL: gsl/ref/gsl-ref_23.html#SEC364
 *               */

/* 
 *   Compile and link with:
 *       gcc -c J0_test.c
 *           gcc -o J0_test J0_test.o -lgsl -lgslcblas -lm
 *           */    
    
/* The answer should be J0(5) = -1.775967713143382920e-01 */


#include <stdio.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_legendre.h>
int
main (void)
{
  int l, m;
  double x = 1.523175;
  double cosx;

  cosx = cos(x);

  l = 1;
  m = 1;

  double y = gsl_sf_legendre_Plm (l, m, cosx);

  printf("P11(%g) = %.18e\n", cosx, y);

  return 0;
}
