/* ############################################################
      
     FILE:     cooling.h

     PURPOSE:  contains common definitions for the 
               whole CODE

     Notice: the order is absolutely important and MUST NOT
             be changed !!!
   ############################################################ */


#define C_IONS  3   /* in [1,5] */
#define N_IONS  3   /* in [1,5] */
#define O_IONS  3   /* in [1,5] */
#define Ne_IONS 3   /* in [1,5] */
#define S_IONS  3   /* in [1,5] */
#define Fe_IONS 0   /* in [0,3] */

#if C_IONS == 0
 #define C_EXPAND(a,b,c,d,e)  
#elif C_IONS == 1       
 #define C_EXPAND(a,b,c,d,e)  ,a
#elif C_IONS == 2     
 #define C_EXPAND(a,b,c,d,e)  ,a,b
#elif C_IONS == 3       
 #define C_EXPAND(a,b,c,d,e)  ,a,b,c 
#elif C_IONS == 4
 #define C_EXPAND(a,b,c,d,e)  ,a,b,c,d
#elif C_IONS == 5      
 #define C_EXPAND(a,b,c,d,e)  ,a,b,c,d,e
#endif

#if N_IONS == 0
 #define N_EXPAND(a,b,c,d,e)  
#elif N_IONS == 1      
 #define N_EXPAND(a,b,c,d,e)  ,a
#elif N_IONS == 2     
 #define N_EXPAND(a,b,c,d,e)  ,a,b
#elif N_IONS == 3       
 #define N_EXPAND(a,b,c,d,e)  ,a,b,c 
#elif N_IONS == 4
 #define N_EXPAND(a,b,c,d,e)  ,a,b,c,d
#elif N_IONS == 5      
 #define N_EXPAND(a,b,c,d,e)  ,a,b,c,d,e
#endif

#if O_IONS == 0
 #define O_EXPAND(a,b,c,d,e)  
#elif O_IONS == 1      
 #define O_EXPAND(a,b,c,d,e)  ,a
#elif O_IONS == 2     
 #define O_EXPAND(a,b,c,d,e)  ,a,b
#elif O_IONS == 3       
 #define O_EXPAND(a,b,c,d,e)  ,a,b,c 
#elif O_IONS == 4
 #define O_EXPAND(a,b,c,d,e)  ,a,b,c,d
#elif O_IONS == 5      
 #define O_EXPAND(a,b,c,d,e)  ,a,b,c,d,e
#endif

#if Ne_IONS == 0
 #define Ne_EXPAND(a,b,c,d,e)  
#elif Ne_IONS == 1      
 #define Ne_EXPAND(a,b,c,d,e)  ,a
#elif Ne_IONS == 2     
 #define Ne_EXPAND(a,b,c,d,e)  ,a,b
#elif Ne_IONS == 3       
 #define Ne_EXPAND(a,b,c,d,e)  ,a,b,c 
#elif Ne_IONS == 4
 #define Ne_EXPAND(a,b,c,d,e)  ,a,b,c,d
#elif Ne_IONS == 5      
 #define Ne_EXPAND(a,b,c,d,e)  ,a,b,c,d,e
#endif

#if S_IONS == 0
 #define S_EXPAND(a,b,c,d,e)  
#elif S_IONS == 1      
 #define S_EXPAND(a,b,c,d,e)  ,a
#elif S_IONS == 2     
 #define S_EXPAND(a,b,c,d,e)  ,a,b
#elif S_IONS == 3       
 #define S_EXPAND(a,b,c,d,e)  ,a,b,c 
#elif S_IONS == 4
 #define S_EXPAND(a,b,c,d,e)  ,a,b,c,d
#elif S_IONS == 5      
 #define S_EXPAND(a,b,c,d,e)  ,a,b,c,d,e
#endif

#if Fe_IONS == 0
 #define Fe_EXPAND(a,b,c)  
#elif Fe_IONS == 1      
 #define Fe_EXPAND(a,b,c)  ,a
#elif Fe_IONS == 2     
 #define Fe_EXPAND(a,b,c)  ,a,b
#elif Fe_IONS == 3       
 #define Fe_EXPAND(a,b,c)  ,a,b,c 
#endif

/* **********************************************************************
     Ions are labeled progressively, depending on how many ionization 
     stages are effectively included in the network through the previous 
     X_EXPAND macros. 
     Elements are ordered as {H, He, C, N, O, Ne, S, Fe} and must be 
     carefully respected everywhere in the code. 
     Hydrogen and Helium are always included.
   ********************************************************************** */

enum {
  HI = NFLX, HeI, HeII
  C_EXPAND(CI, CII, CIII, CIV, CV)
  N_EXPAND(NI, NII, NIII, NIV, NV)
  O_EXPAND(OI, OII, OIII, OIV, OV)
  Ne_EXPAND(NeI, NeII, NeIII, NeIV, NeV)
  S_EXPAND(SI, SII, SIII, SIV, SV)
  Fe_EXPAND(FeI, FeII, FeIII)
};

#define NIONS  (3+C_IONS+N_IONS+O_IONS+Ne_IONS+S_IONS+Fe_IONS)

/*
#define HI    (NFLX)
#define HeI   (NFLX + 1)
#define HeII  (NFLX + 2)
#define CI    (NFLX + 3)
#define CII   (NFLX + 4)
#define CIII  (NFLX + 5)
#define CIV   (NFLX + 6)
#define CV    (NFLX + 7)
#define NI    (NFLX + 8)
#define NII   (NFLX + 9)
#define NIII  (NFLX + 10)
#define NIV   (NFLX + 11)
#define NV    (NFLX + 12)
#define OI    (NFLX + 13)
#define OII   (NFLX + 14)
#define OIII  (NFLX + 15)
#define OIV   (NFLX + 16)
#define OV    (NFLX + 17)
#define NeI   (NFLX + 18)
#define NeII  (NFLX + 19)
#define NeIII (NFLX + 20)
#define NeIV  (NFLX + 21)
#define NeV   (NFLX + 22)
#define SI    (NFLX + 23)
#define SII   (NFLX + 24)
#define SIII  (NFLX + 25)
#define SIV   (NFLX + 26)
#define SV    (NFLX + 27)
#if INCLUDE_Fe == YES
 #define FeI    (NFLX + 28)
 #define FeII   (NFLX + 29)
 #define FeIII  (NFLX + 30)
#endif
*/

real GetMaxRate          (double *, double *, double);
real MeanMolecularWeight (double *);
double H_MassFrac (void);
real CompEquil            (double, double, double *);
real find_N_rho ();
void Radiat (double *, double *);

void CHECK_NORMALIZATION (double *, char *);
void NORMALIZE_IONS (double *);
/*
int Ros4_expl, Ros4_impl, Ros4_sup_dt;
*/








