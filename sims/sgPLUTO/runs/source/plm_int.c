#include <stdio.h>
     #include <math.h>
     #include <gsl/gsl_integration.h>
     #include <gsl/gsl_sf_legendre.h>
     
     double f (double x, void * params) {
       double alpha = *(double *) params;
       int l=8, m=8;
       double cosx;
       cosx = cos(x);
       double f = gsl_sf_legendre_sphPlm(l, m, cosx);
       return f;
     }
     
     int
     main (void)
     {
       gsl_integration_workspace * w 
         = gsl_integration_workspace_alloc (4000);
  
       size_t    neval;     
       double result, error;
       double a=0.0, b=cos(1.27);    

       gsl_function F;
       F.function = &f;

/*     
       gsl_integration_qags (&F, a, b, 1.0e-9, 1e-9, 4000,
                             w, &result, &error); 
*/
/*      
       gsl_integration_qag(&F, a, b, 1.0e-7, 1.0e-7, 4000, 1, w, &result, &error); 
*/

       gsl_integration_qng(&F, a, b, 1.0e-6, 1.0e-6, &result, &error, &neval);

    
       printf ("result          = % .18f\n", result);
       printf ("estimated error = % .18f\n", error);
     
       gsl_integration_workspace_free (w);
     
       return 0;
     }

