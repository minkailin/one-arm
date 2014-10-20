/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief Disk-Planet problem.

  Disk-Planet setup.
 
  Reference paper:
   "A Conservative orbital advection scheme for simulations
    of magnetized shear flows with the PLUTO Code"
    Mignone et al, A&A (2012)
 
  -------------------------------------------------------------
   Independent of the coordinate system (polar/spherical), we
   adopt the following conventions:
 
    r = spherical radius
    R = cylindrical radius
    z = cylindrical height
   th = meridional angle

  \author A. Mignone (mignone@ph.unito.it)
  \date   Aug 16, 2012
*/
/* ///////////////////////////////////////////////////////////////////// */

#include "pluto.h"
#include "float.h"
#define MIN_DENSITY 1.0e-9
#define g_OmegaZ  0.0
#define bigG 1.0
#define mstar 1.0


static double ran2(long int *idum);

/*MKL: accretion disk functions*/
static double grav_pot3_cylin(const double bigR, const double z);
static double azivel(const double bigR, const double z);
static double domega0(const double bigR, const double z);
static double dlog_rho0(const double bigR);
static double inner_hole(const double bigR);
static double dloghole(const double bigR);
static double d2loghole(const double bigR);
static double radvel(const double x1, const double x2);
static double polvel(const double x1, const double x2);

int pert;
double sig0, bwidth, dampin, dampout, planeton, switchon, planeton_switchon, min_rho;

/* ********************************************************************* */
void Init (double *us, double x1, double x2, double x3)
/*
 *
 *
 *
 *********************************************************************** */
{
  double r, th, R, z, H, OmegaK, cs;
  double scrh;
  long int iseed = -1;

  iseed = -1 - (int) x1*x2*x3*1.0e4;

  #if EOS == IDEAL
   g_gamma = g_inputParam[gmma];
  #endif

  g_unitLength   = g_inputParam[r0];
  g_unitDensity  = density3D(g_unitLength, 0.0);
  g_unitVelocity = g_unitLength*omega_k(g_unitLength);

  
  #if GEOMETRY == POLAR
   R  = x1;
   #if DIMENSIONS == 2
    z  = 0.0;
    r  = R;
    th = 0.5*CONST_PI;
   #else
    z  = x3;
    r  = sqrt(R*R + z*z);
    th = atan2(R,z);
   #endif
  #elif GEOMETRY == SPHERICAL
   r  = x1;
   th = x2;
   R  = r*sin(th);
   z  = r*cos(th);
  #endif
  
  if(g_inputParam[visc_steady] < 0){//qout is toomre at outer boundary
  sig0      = sqrt(csq(g_inputParam[rout]))*omega_k(g_inputParam[rout]);
  sig0     /=
    g_inputParam[qout]*CONST_PI*bigG*pow(g_inputParam[rout]/g_unitLength,
					 -g_inputParam[smallp])*inner_hole(g_inputParam[rout]); 

  double nmax;
  nmax  = g_inputParam[max_H]*g_inputParam[smallh]*g_inputParam[rout]/bigH(g_inputParam[rout]);
  if(sg_on > 0.0) sig0/=erf(nmax/sqrt(2.0));
  }

  min_rho = MIN_DENSITY*sig0/bigH(g_unitLength);

  if(g_inputParam[visc_steady] > 0){//qout is outer boundary of domain or dead zone
  double rmax;

  if(g_inputParam[densbeta] < 0.0){
  rmax = g_inputParam[rout]; 
  if( g_inputParam[visc_rdead2] < g_inputParam[rout]) rmax = g_inputParam[visc_rdead2]; //if outer dead edge within domain, then Q is Q at outer edge of dead zone
  sig0      = sqrt(csq(rmax))*omega_k(rmax);
  sig0     /=
    g_inputParam[qout]*CONST_PI*bigG*pow(rmax/g_unitLength,
                                         -g_inputParam[smallp])*inner_hole(rmax);
  } else {
  double rs;
  rmax = g_unitLength;
  rs = pow(3.0*g_inputParam[densbeta]/g_inputParam[smallp], -1.0/3.0)*g_unitLength; 
  sig0      = sqrt(csq(rmax))*omega_k(rmax);
  sig0     /=
    g_inputParam[qout]*CONST_PI*bigG*pow(rs/rmax,
                                         g_inputParam[smallp])*inner_hole(rmax);   
  }

  double nmax;
  nmax  = g_inputParam[max_H]*g_inputParam[smallh]*rmax/bigH(rmax);
  if(sg_on > 0.0) sig0/=erf(nmax/sqrt(2.0));
  }

  dampin = g_inputParam[damp_in]*g_inputParam[rmin];
  dampout= g_inputParam[damp_out]*g_inputParam[rout];
  planeton= g_inputParam[planet_on]*2.0*CONST_PI;
  switchon= g_inputParam[switch_on]*2.0*CONST_PI;
  planeton_switchon = planeton + switchon;  
    
  us[RHO] = density3D(R, z);
  us[VX1] = us[VX2] = 0.0;
  
  pert = 0;
 
  //background vel field (non-SG)
  us[iVPHI] = azivel(R, z);

  #if EOS == IDEAL
   us[PRS] = us[RHO]*csq(R);
  #elif EOS == ISOTHERMAL
   g_isoSoundSpeed = sqrt(csq(g_unitLength));
  #endif
}

static double grav_pot3_cylin(const double bigR, const double z)
{
  double d;
  double star_pot;
  
  d =sqrt(bigR*bigR + z*z);

  star_pot  = -bigG*mstar/d;
  return star_pot;
}

double omega_k(const double bigR)
{
  double omega_sq;

  omega_sq = bigG*mstar/pow(bigR, 3.0);
  return sqrt(omega_sq);
}

double csq(const double bigR)
{
  double soundspeed_sq;

  soundspeed_sq = pow(g_inputParam[smallh], 2.0);
  soundspeed_sq*= bigG*mstar/g_unitLength; 
  soundspeed_sq*= pow(bigR/g_unitLength, -g_inputParam[smallq]); 

  return soundspeed_sq;
}

double bigH(const double bigR)
{
  return sqrt(csq(bigR))/omega_k(bigR);
}

static double inner_hole(const double bigR)
{  
  double res; 

  res = 1.0; 

  if(g_inputParam[visc_steady] < 0){ 
  double Hin, Rmin, sintheta;
  sintheta = 1.0/sqrt(1.0 +
		      pow(g_inputParam[max_H]*g_inputParam[smallh],2.0));
  Rmin = g_inputParam[rmin]*sintheta;			
  Hin  = g_inputParam[HIN]*bigH(Rmin);
  res = 1.0 - sqrt(Rmin/(bigR + Hin));
  }

  if(g_inputParam[visc_steady] > 0){
  if(g_inputParam[densbeta] < 0.0){
  double dr1, fr1, dr2, fr2, eps;

  eps = 1.0/g_inputParam[visc_amp]; 

  dr1 = g_inputParam[visc_width]*bigH(g_inputParam[visc_rdead1]);
  fr1 = (1.0 - eps)*(1.0 + tanh((bigR-g_inputParam[visc_rdead1])/dr1))/2.0 + eps;

  dr2 = g_inputParam[visc_width]*bigH(g_inputParam[visc_rdead2]);
  fr2 = (1.0 - eps)*(1.0 - tanh((bigR-g_inputParam[visc_rdead2])/dr2))/2.0 + eps;

  res = fr1*fr2;
  } else {
  double rs;
  rs = pow(3.0*g_inputParam[densbeta]/g_inputParam[smallp], -1.0/3.0)*g_unitLength;

  res = exp(-g_inputParam[densbeta]*pow(rs/bigR, 3.0));
  }
  }
  return res;
}

static double dloghole(const double bigR)
{  
  double res;
  
  res = 0.0;   

  if(g_inputParam[visc_steady] < 0){
  double Hin, Rmin, sintheta, dhole;
  sintheta = 1.0/sqrt(1.0 +
		      pow(g_inputParam[max_H]*g_inputParam[smallh],2.0));
  Rmin = g_inputParam[rmin]*sintheta;			
  Hin  = g_inputParam[HIN]*bigH(Rmin);
  dhole=  0.5*sqrt(Rmin)*pow(bigR + Hin, -1.5);
  res  =  dhole/inner_hole(bigR);
  }

 if(g_inputParam[visc_steady] > 0){
 if(g_inputParam[densbeta] < 0.0){
 double eps, dr1, dfr1, dr2, dfr2, sech1, sech2;

 eps = 1.0/g_inputParam[visc_amp];

  dr1 = g_inputParam[visc_width]*bigH(g_inputParam[visc_rdead1]);
  sech1 = 1.0/cosh( (bigR-g_inputParam[visc_rdead1])/dr1);
  dfr1 = (1.0 - eps)*pow(sech1, 2.0);
  dfr1/= 2.0*dr1;

  dr2 = g_inputParam[visc_width]*bigH(g_inputParam[visc_rdead2]);
  sech2 =  1.0/cosh( (bigR-g_inputParam[visc_rdead2])/dr2);
  dfr2 = -(1.0 - eps)*pow(sech2, 2.0);
  dfr2/= 2.0*dr2;

  res = (dfr1 + dfr2)/inner_hole(bigR);
 } else {
  double rs;
  rs = pow(3.0*g_inputParam[densbeta]/g_inputParam[smallp], -1.0/3.0)*g_unitLength;
  res =3.0*g_inputParam[densbeta]*pow(rs/bigR, 3.0)/bigR;
  }
  }
  return res;
}

static double d2loghole(const double bigR)
{  
  double res;

  res = 0.0; 

  if(g_inputParam[visc_steady] < 0){
  double Hin, Rmin, sintheta, d2hole;
  sintheta = 1.0/sqrt(1.0 +
		      pow(g_inputParam[max_H]*g_inputParam[smallh],2.0));
  Rmin = g_inputParam[rmin]*sintheta;			
  Hin  = g_inputParam[HIN]*bigH(Rmin);
  d2hole =  -(3.0/4.0)*sqrt(Rmin)*pow(bigR+Hin,-2.5);
  res =  d2hole/inner_hole(bigR) - pow(dloghole(bigR), 2.0);
  }

  if(g_inputParam[visc_steady] > 0){
  if(g_inputParam[densbeta] < 0.0){
  double eps, dr1, dr2, d2fr1, d2fr2, sech1, sech2, arg1, arg2;

  eps = 1.0/g_inputParam[visc_amp];

  dr1 = g_inputParam[visc_width]*bigH(g_inputParam[visc_rdead1]);
  arg1 = (bigR-g_inputParam[visc_rdead1])/dr1;
  sech1 = 1.0/cosh(arg1);
  d2fr1 = (1.0 - eps)*pow(sech1, 2.0)*tanh(arg1);
  d2fr1/= -dr1*dr1;

  dr2 = g_inputParam[visc_width]*bigH(g_inputParam[visc_rdead2]);
  arg2 = (bigR-g_inputParam[visc_rdead2])/dr2;
  sech2 = 1.0/cosh(arg2);
  d2fr2 = -(1.0 - eps)*pow(sech2, 2.0)*tanh(arg2);
  d2fr2/= -dr2*dr2;

  res = (d2fr1+d2fr2)/inner_hole(bigR) - pow(dloghole(bigR), 2.0);
  } else {
  double rs;
  rs = pow(3.0*g_inputParam[densbeta]/g_inputParam[smallp], -1.0/3.0)*g_unitLength;  
  res = -12.0*g_inputParam[densbeta]*pow(rs/bigR, 3.0)/bigR/bigR;
  }
  }
  return res; 
}


double surface_density(const double bigR)
{
  double sig_s;
  
  if(g_inputParam[densbeta] < 0.0){  
  sig_s = sig0*pow(bigR/g_unitLength, -g_inputParam[smallp]);
  } else {
  double rs;
  rs = pow(3.0*g_inputParam[densbeta]/g_inputParam[smallp], -1.0/3.0)*g_unitLength;
  sig_s = sig0*pow(rs/bigR, g_inputParam[smallp]); 
  }
  return sig_s*inner_hole(bigR);
}

double density3D(const double bigR, const double z)
{
  double vertical, sigma;

  vertical = -(grav_pot3_cylin(bigR, z) - grav_pot3_cylin(bigR,0.0));
  vertical/= csq(bigR); 
 
  sigma = surface_density(bigR);
  sigma/= sqrt(2.0*CONST_PI)*bigH(bigR);
  
  return sigma*exp(vertical);
}

static double dlog_rho0(const double bigR)
{
  double dsigma_s, dbump, dbigH, dtaper, dtaper_in;
  double dcs, domega_k, Hin, Rmin, sintheta;

  //power law contribution
  dsigma_s = -g_inputParam[smallp]/bigR;

  //inner hole contribution
  dsigma_s+= dloghole(bigR);
    
  //sound-speed gradients
  dcs     =  -0.5*g_inputParam[smallq]/bigR;
  domega_k = -1.5/bigR;

  dbigH   = dcs - domega_k;
  
  return dsigma_s - dbigH;
}

static double azivel(const double bigR, const double z)
{
  double dcs2, cs2, dphi0, vsq; 
  
  cs2     = csq(bigR);
  dcs2    =-g_inputParam[smallq]/bigR;
  dphi0   = bigR*pow(omega_k(bigR),2.0);
  
  vsq     = cs2*dlog_rho0(bigR) +
    dcs2*(grav_pot3_cylin(bigR,z)-grav_pot3_cylin(bigR,0.0)+cs2) +
    dphi0; 
  vsq    *= bigR;
  return sqrt(vsq); 
}

static double domega0(const double bigR, const double z)
{
  double cs2, omega, dlogH, d2logH, dPhi, dPhi0, d2Phi0; 
  double result, dlogcs2, d2logcs2, dlogsigma, d2logsigma;  

  cs2     =  csq(bigR);  
  dlogcs2 = -g_inputParam[smallq]/bigR;
  d2logcs2=  g_inputParam[smallq]/(bigR*bigR);
 
  dlogsigma = -g_inputParam[smallp]/bigR + dloghole(bigR); 
  d2logsigma=  g_inputParam[smallp]/(bigR*bigR) + d2loghole(bigR);

  dlogH = (3.0 - g_inputParam[smallq])/(2.0*bigR); 
  d2logH=-(3.0 - g_inputParam[smallq])/(2.0*bigR*bigR); 

  dPhi = bigG*mstar*bigR/pow(bigR*bigR + z*z, 1.5);  
  dPhi0= bigG*mstar/pow(bigR, 2.0); 
  d2Phi0=-2.0*bigG*mstar/pow(bigR, 3.0); 

  omega = azivel(bigR, z)/bigR;

  result = cs2*dlogcs2*(dlogcs2 + dlogsigma - dlogH);
  result+= cs2*(d2logcs2 + d2logsigma - d2logH);
  result+= d2logcs2*(grav_pot3_cylin(bigR, z) - grav_pot3_cylin(bigR, 0.0) );
  result+= dlogcs2*(dPhi - dPhi0);
  result+= d2Phi0;   

  result-= omega*omega;
  result/= 2.0*omega*bigR;

  return result;
}


double kinematic_visc(const double x1, const double x2)
{
  /*assume spherical polars*/
  double z, fz, R, fr ;
  double visc, Rmin, sintheta, Rmax;

  R = x1*sin(x2);
  z = x1*cos(x2);
 
  sintheta = 1.0/sqrt(1.0 +
                      pow(g_inputParam[max_H]*g_inputParam[smallh],2.0));
  Rmin = g_inputParam[rmin]*sintheta;

  Rmax = g_inputParam[rout]; 

  //reference viscosity value  
   if(g_inputParam[nu_alpha] > 0.0) {
    visc  = g_inputParam[nu_alpha]*sqrt(csq(Rmin))*bigH(Rmin);
   } 
   if(g_inputParam[nu] > 0.0){
    visc  = g_inputParam[nu]*g_unitLength*g_unitLength*omega_k(g_unitLength); 
   }
  

  fz = 1.0;

  if(g_inputParam[visc_steady] < 0){ //we're not starting with steady state, impose fixed structured viscosity profile
  double fr, drnu1, drnu2, fz1, fz2; 

  // nu.sigma.r^3.omega' = constant. reference viscosity is dead zone viscosity 
  
  visc *= surface_density(Rmin)*pow(Rmin,3.0)*domega0(Rmin,0.0);
  visc /= surface_density(R)*pow(R,3.0)*domega0(R,0.0);

  drnu1 = g_inputParam[visc_width]*bigH(g_inputParam[visc_rdead1]);
  drnu2 = g_inputParam[visc_width]*bigH(g_inputParam[visc_rdead2]);

  fz1 = 1.0;
  fz2 = 1.0;

  //viscosity jump
  if(g_inputParam[visc_rdead1] > Rmin){//impose an inner viscosity edge 
    fr   = 0.5*(1.0 + tanh( (R - g_inputParam[visc_rdead1])/drnu1 ) );
    fz1 =  fr +  (1.0 - fr)*g_inputParam[visc_amp];
  }
  
  if(g_inputParam[visc_rdead2] < g_inputParam[rout]){//impose an outer viscosity edge 
    fr  = 0.5*(1.0 + tanh((R - g_inputParam[visc_rdead2])/drnu2) );
    fz2 = fr*g_inputParam[visc_amp] + (1.0 - fr);
  }
    fz = fz1*fz2;  
  }

  if(g_inputParam[visc_steady] > 0){//we're starting with steady state. 
  
  // nu.sigma.r^3.omega' = constant. reference viscosit is outer edge visc
  visc *= surface_density(Rmax)*pow(Rmax,3.0)*domega0(Rmax,0.0);
  visc /= surface_density(R)*pow(R,3.0)*domega0(R,0.0);
  }

  return visc*fz;
}

double pep_eos(const double r, const double theta, const double phi)
{
 double const hp=0.5, eosn=3.5;
 double Hp, dp, x, y, z, bigR, t, omega_p, rsm, xp, yp;
 double csmod;

 if(g_inputParam[mplanet] > 0.0) {
 bigR = r*sin(theta);
 x = bigR*cos(phi);
 y = bigR*sin(phi);
 z = r*cos(theta);
 
 t = g_time;

  if(g_intStage == 2){
#if TIME_STEPPING == RK2 
    t += 0.5*g_dt;
#elif TIME_STEPPING == RK3
    t += 0.25*g_dt;
#endif
  }
  if(g_intStage == 3) t += 1.0/3.0*g_dt;

   omega_p = omega_k(g_unitLength);
   xp =  g_unitLength*cos(omega_p*t + CONST_PI);
   yp =  g_unitLength*sin(omega_p*t + CONST_PI);

   rsm = g_unitLength*pow(g_inputParam[mplanet]/mstar/3.0, 1.0/3.0);
   rsm *= g_inputParam[soft];
   dp = sqrt((x-xp)*(x-xp) + (y-yp)*(y-yp) + z*z + rsm*rsm);

   Hp = hp*dp;

   csmod = Hp*sqrt( 1.0 + (g_inputParam[mplanet]/mstar)*pow(bigR/dp, 3.0) );
   csmod/=pow( pow(bigH(bigR), eosn) + pow(Hp, eosn), 1.0/eosn);
   csmod-= 1.0;
   csmod*= turn_on(t);
   csmod+= 1.0; 
  } else {
   csmod = 1.0;
  }
  
   return csmod;
}

double turn_on(const double time)
{
  double temp;
  
  if(time < planeton){

     return 0.0;
   } else if( (time >= planeton) && (time < planeton_switchon)){
     temp = sin(0.5*CONST_PI*(time - planeton)/switchon);
     temp*= temp;

     return temp;
   } else {

     return 1.0;
   }
} 


/* ********************************************************************* */
void Analysis (const Data *d, Grid *grid)
/* 
 *
 *
 *********************************************************************** */
{

}
/* ********************************************************************* */
void UserDefBoundary (const Data *d, RBox *box, int side, Grid *grid) 
/*! 
 *  Assign user-defined boundary conditions.
 *
 * \param [in,out] d  pointer to the PLUTO data structure containing
 *                    cell-centered primitive quantities (d->Vc) and 
 *                    staggered magnetic fields (d->Vs, when used) to 
 *                    be filled.
 * \param [in] box    pointer to a RBox structure containing the lower
 *                    and upper indices of the ghost zone-centers/nodes
 *                    or edges at which data values should be assigned.
 * \param [in] side   specifies the boundary side where ghost zones need
 *                    to be filled. It can assume the following 
 *                    pre-definite values: X1_BEG, X1_END,
 *                                         X2_BEG, X2_END, 
 *                                         X3_BEG, X3_END.
 *                    The special value side == 0 is used to control
 *                    a region inside the computational domain.
 * \param [in] grid  pointer to an array of Grid structures.
 *
 *********************************************************************** */
{
  int   i, j, k, nv, counter;
  double *x1, *x2, *x3, *x1r, *x1l, *dx1, *x2r, *dx2, R,vOmegaK, v[256];
  static int do_once = 1;
  double r, z, damp_time, vphi_zero, d_zero, p_zero;
  double vphi_bound, vphi_sg1, vphi_sg2, vphi_sg; 
  double taper, dt;  
 
  x1 = grid[IDIR].x;
  x2 = grid[JDIR].x;
  x3 = grid[KDIR].x;

  x1r=grid[IDIR].xr;
  x1l=grid[IDIR].xl;
  dx1=grid[IDIR].dx;

  x2r=grid[JDIR].xr;
  dx2=grid[JDIR].dx;

  if(g_intStage == 1) dt = g_dt;
  if(g_intStage == 2){
#if TIME_STEPPING == RK2 
    dt = 0.5*g_dt;
#elif TIME_STEPPING == RK3
    dt = 0.25*g_dt;
#endif
  }
  if(g_intStage == 3) dt = 2.0/3.0*g_dt;

  
   if (side == X1_BEG){
    
    //modified open boundary 
    X1_BEG_LOOP(k,j,i){ 
      r = x1[i]; 
      R = r*sin(x2[j]); 
      z = r*cos(x2[j]);

      //modified open boundary      
/*      
      for (nv = 0; nv < NVAR; nv++){ 
	d->Vc[nv][k][j][i] = d->Vc[nv][k][j][IBEG]; 
      } 
      if( (d->Vc[VX1][k][j][IBEG] > 0.0) ) d->Vc[VX1][k][j][i] = 0.0;  
*/      

      //unperturbed
      d->Vc[RHO][k][j][i] = density3D(R, z)*d->beta2d[j][i]; 
      d->Vc[VX1][k][j][i] = 0.0;
      d->Vc[VX2][k][j][i] = 0.0;

      //star and temp. grad. 
      vphi_bound = -g_inputParam[smallq]*csq(R); 
      vphi_bound+=  mstar*bigG/r;
      
      //self-g (on fly)
      vphi_sg1 = (d->Phi[k][j][IBEG] - d->Phi[k][j][IBEG-1])/(x1[IBEG] - x1[IBEG-1]); 
      vphi_sg2 = (d->Phi[k][j][IBEG+1] - d->Phi[k][j][IBEG])/(x1[IBEG+1] - x1[IBEG]); 
      vphi_sg  = vphi_sg1 + (r - x1l[IBEG])*(vphi_sg2 - vphi_sg1)/dx1[IBEG];       
      vphi_bound += r*vphi_sg;
 
      //den. grad. 
      vphi_sg1 = log( d->Vc[RHO][k][j][IBEG]/d->Vc[RHO][k][j][IBEG-1] )/(x1[IBEG] - x1[IBEG-1]); 
      vphi_sg2 = log( d->Vc[RHO][k][j][IBEG+1]/d->Vc[RHO][k][j][IBEG] )/(x1[IBEG+1] - x1[IBEG]);
      vphi_sg  = vphi_sg1 + (r - x1l[IBEG])*(vphi_sg2 - vphi_sg1)/dx1[IBEG];     
      vphi_bound += r*vphi_sg*csq(R);       
  
      d->Vc[VX3][k][j][i] = sqrt(vphi_bound);     
  
    } 
    
   }

   if (side == X1_END){
    
    //modified open boundary 
    X1_END_LOOP(k,j,i){ 
      r = x1[i]; 
      R = r*sin(x2[j]); 
      z = r*cos(x2[j]);

      //modified open boundary      
/*         
      for (nv = 0; nv < NVAR; nv++){ 
	d->Vc[nv][k][j][i] = d->Vc[nv][k][j][IEND]; 
      } 
      if( (d->Vc[VX1][k][j][IEND] < 0.0) ) d->Vc[VX1][k][j][i] = 0.0;  
*/      

      //unperturbed
      d->Vc[RHO][k][j][i] = density3D(R, z)*d->beta2d[j][i]; 
      d->Vc[VX1][k][j][i] = 0.0;
      d->Vc[VX2][k][j][i] = 0.0;

      //star and temp. grad. 
      vphi_bound = -g_inputParam[smallq]*csq(R); 
      vphi_bound+= mstar*bigG/r;
      
      //self-g (on fly)
      vphi_sg1 = (d->Phi[k][j][IEND] - d->Phi[k][j][IEND-1])/(x1[IEND] - x1[IEND-1]); 
      vphi_sg2 = (d->Phi[k][j][IEND+1] - d->Phi[k][j][IEND])/(x1[IEND+1] - x1[IEND]); 
      vphi_sg  = vphi_sg1 + (r - x1l[IEND])*(vphi_sg2 - vphi_sg1)/dx1[IEND];       
      vphi_bound += r*vphi_sg;
 
      //den. grad. 
      vphi_sg1 = log( d->Vc[RHO][k][j][IEND]/d->Vc[RHO][k][j][IEND-1] )/(x1[IEND] - x1[IEND-1]); 
      vphi_sg2 = log( d->Vc[RHO][k][j][IEND+1]/d->Vc[RHO][k][j][IEND] )/(x1[IEND+1] - x1[IEND]);
      vphi_sg  = vphi_sg1 + (r - x1l[IEND])*(vphi_sg2 - vphi_sg1)/dx1[IEND];     
      vphi_bound += r*vphi_sg*csq(R);       
  
      d->Vc[VX3][k][j][i] = sqrt(vphi_bound);        
      
    } 
    
   }
  

/*
  if (side == X1_BEG){
    
    //modified open boundary 
    X1_BEG_LOOP(k,j,i){ 
      r = x1[i]; 
      R = r*sin(x2[j]); 
      z = r*cos(x2[j]);
     
      for (nv = 0; nv < NVAR; nv++){ 
	d->Vc[nv][k][j][i] = d->Vc[nv][k][j][IBEG]; 
      } 
          
      vphi_bound = -g_inputParam[smallq]*csq(R); 
      vphi_bound+= mstar*bigG/r; 
      
      //linear extrapolation for dphi/dr  
      vphi_sg1 = (d->Phi[k][j][IBEG] - d->Phi[k][j][IBEG-1])/(x1[IBEG] - x1[IBEG-1]); 
      vphi_sg2 = (d->Phi[k][j][IBEG+1] - d->Phi[k][j][IBEG])/(x1[IBEG+1] - x1[IBEG]); 
      vphi_sg  = vphi_sg1 - (x1l[IBEG] - r)*(vphi_sg2 - vphi_sg1)/dx1[IBEG];       
      vphi_bound += r*vphi_sg;       
      d->Vc[VX3][k][j][i] = sqrt(vphi_bound); 
      
      if( (d->Vc[VX1][k][j][IBEG] > 0.0) ) d->Vc[VX1][k][j][i] = 0.0;  
      
    } 
    
  }
  

  if (side == X1_END){
    double cs, cs2, dW, dU, dV;
    double z1, z2, z3;
   
     //godon nrbc 
     int ip1, im1; 
     double rho_a, rho_e, vphi_a, vphi_e, vrad_a, vtheta_a; 
     double cs_g, cs_e, cs_a, Pplus, Pminus, grad;  
     
    X1_END_LOOP(k,j,i){//if outgoing wave, set to active characteristic, if incoming, set variables to extrapolated (above)

       //first set ghosts to vrad=0, vtheta=0, rho and vphi linear extrapolate 
       grad = (x1[i] - x1[IEND-1])/(x1[IEND] - x1[IEND-1]);

       d->Vc[RHO][k][j][i] = d->Vc[RHO][k][j][IEND];//d->Vc[RHO][k][j][IEND-1] + (d->Vc[RHO][k][j][IEND] - d->Vc[RHO][k][j][IEND-1])*grad;
       d->Vc[VX1][k][j][i] = 0.0;
       d->Vc[VX2][k][j][i] = 0.0;
       d->Vc[VX3][k][j][i] = d->Vc[VX3][k][j][IEND];//d->Vc[VX3][k][j][IEND-1] + (d->Vc[VX3][k][j][IEND] - d->Vc[VX3][k][j][IEND-1])*grad;

       //now apply NRBC 

       im1 = i - 1 ;    //last active cell 
       ip1 = i     ;    //when setting ith cell, should use i-1 cell for outing (above, the prev. cel) and i+1 for incoming. but here use i instead of i+1

     	//last active cell quantities 
     	r        = x1[im1]; 
     	R        = r*sin(x2[j]); 
     	cs_a     = sqrt(csq(R)); 
     	rho_a    = d->Vc[RHO][k][j][im1]; 
     	vrad_a   = d->Vc[VX1][k][j][im1]; 
     	vtheta_a = d->Vc[VX2][k][j][im1]; 
     	vphi_a   = d->Vc[VX3][k][j][im1]; 
	
     	//ghost quantities 
     	r        = x1[i]; 
     	R        = r*sin(x2[j]); 
     	cs_g     = sqrt(csq(R)); 
	
     	//second ghost eqm quantities 
     	r        = x1[ip1]; 
     	R        = r*sin(x2[j]); 
     	cs_e     = sqrt(csq(R)); 
     	rho_e    = d->Vc[RHO][k][j][ip1]; 
     	vphi_e   = d->Vc[VX3][k][j][ip1]; 
	
     	//rho and vrad 
     	if( fabs(vrad_a) < cs_a ){//subsonic radial flow 
     	  d->Vc[RHO][k][j][i] = pow(rho_a, 0.5*cs_a/cs_g)*pow(rho_e,0.5*cs_e/cs_g)*exp(0.5*vrad_a/cs_g); 
     	  d->Vc[VX1][k][j][i] =0.5*( vrad_a + cs_a*log(rho_a) - cs_e*log(rho_e) ); 
     	} else { //supersonic 
     	  if( vrad_a > 0.0){//copy from active 
     	    d->Vc[RHO][k][j][i] = rho_a; 
     	    d->Vc[VX1][k][j][i] = vrad_a; 
     	  } 
     	  if( vrad_a < 0.0){//set to eqm 
     	    d->Vc[RHO][k][j][i] = rho_e; 
     	    d->Vc[VX1][k][j][i] = 0.0; 
     	  } 
     	} 
	

        //vphi according to vrad 
     	 if(vrad_a > 0.0){//copy from active  
           d->Vc[VX2][k][j][i] = vtheta_a; 
     	   d->Vc[VX3][k][j][i] = vphi_a;  
     	 }  
     	 if(vrad_a < 0.0){//set to eqm  
           d->Vc[VX2][k][j][i] = 0.0;
     	   d->Vc[VX3][k][j][i] = vphi_e;  
     	 }  
		
       } 

  }
*/  

  if (side == X2_BEG){
    double rhs1, rhs2, rhs3, rhs;
    double cs2, Phip1, Phip2, Phip, dvphi;
     
    X2_BEG_LOOP(k,j,i){

      dvphi = (d->Vc[VX3][k][JBEG+1][i] - d->Vc[VX3][k][JBEG][i])/(x2[JBEG+1]-x2[JBEG]);

//extrapolate vphi, set zero meridional flow   
      d->Vc[VX3][k][j][i] = d->Vc[VX3][k][JBEG][i]
                           - (x2[JBEG] - x2[j])*dvphi;
      d->Vc[VX2][k][j][i] = 0.0;
      d->Vc[VX1][k][j][i] = 0.0;
      
//approximate vertically static atmosphere
  
      R = x1[i]*sin(x2r[j]);
      cs2 = csq(R); 
        
      rhs1 = g_inputParam[smallq]*cot(x2r[j]);
      
      Phip1 = (d->Phi[k][JBEG][i] - d->Phi[k][JBEG-1][i])/(x2[JBEG] - x2[JBEG-1]);
      Phip2 = (d->Phi[k][JBEG+1][i] - d->Phi[k][JBEG][i])/(x2[JBEG+1] - x2[JBEG]);
      Phip  = Phip1 - (x2r[JBEG-1] - x2r[j])*(Phip2 - Phip1)/dx2[JBEG]; 
      rhs2  =-Phip/cs2;
    
      rhs3  = d->Vc[VX3][k][JBEG][i]
                           - (x2[JBEG] - x2r[j])*dvphi;
      rhs3 = rhs3*rhs3;
      rhs3*= cot(x2r[j])/cs2; 

      rhs  = rhs1 + rhs2 + rhs3; 
      rhs = exp( -(x2[j+1] - x2[j])*rhs );
      d->Vc[RHO][k][j][i] = d->Vc[RHO][k][j+1][i]*rhs;

    }
  }

  if (side == 0){
    double dvrad, dr1, dr2, fr, rho0, drho0, rho1, vphi0, vphi1;
    long int iseed, iseed_m;
  
    //time for random pert?
    //  print1 ("time is %f \n",g_time);
    if((planeton_switchon >= g_time) && (planeton_switchon <= g_time + g_dt) && (pert == 0)){
      pert = 1;
      iseed  =  x1[IEND/2]*x2[JEND/2]*x3[KEND/2]*1.0e4;
      iseed_m= -iseed;
      ran2(&iseed_m);
      
      //random pert to radial velocity
      print1 ("do random pert... \n");

      dr1 = (g_inputParam[visc_rdead2]-g_inputParam[visc_rdead1])/2.0;

      JTOT_LOOP(j){
      ITOT_LOOP(i){ 

      KTOT_LOOP(k){
       d->Vc[VX1][k][j][i] = 0.0;
       d->Vc[VX2][k][j][i] = 0.0; 
      } 

	R = x1[i]*sin(x2[j]);
              	 
        fr = g_inputParam[pert_amp]*(2.0*ran2(&iseed) - 1.0); 
        fr*= exp(-0.5*pow((R - 0.5*(g_inputParam[visc_rdead1]+g_inputParam[visc_rdead2]))/dr1, 2.0));
 
        rho0 = 0.0;
        vphi0= 0.0;

        KDOM_LOOP(k){
        rho0 += d->Vc[RHO][k][j][i];   
        vphi0+= d->Vc[VX3][k][j][i];
        }
        rho0 /= (double) KEND-KBEG+1;
        vphi0/= (double) KEND-KBEG+1; 

        KDOM_LOOP(k){
        rho1 = fr*rho0;
        vphi1=-fr*vphi0;

        drho0= -2.0*rho1*vphi1/vphi0; 

        d->Vc[RHO][k][j][i] = rho0 + drho0 + 2.0*cos(x3[k])*rho1; 
        d->Vc[VX3][k][j][i] = vphi0 + 2.0*cos(x3[k])*vphi1; 
        }

      }
    }
   }
 
    TOT_LOOP(k,j,i){
      r = x1[i];
      R = r*sin(x2[j]);
      z = r*cos(x2[j]);
      
      d_zero = density3D(R, z)*d->beta2d[j][i];
      p_zero = d_zero*csq(R);
//      vphi_zero = azivel(R,z);
      
      damp_time = 2.0*CONST_PI*R/vphi_zero;
      damp_time*= g_inputParam[tdamp]; 
      
      if(r <= dampin){
	
	taper = pow( (dampin - r)/(dampin - g_inputParam[rmin]), 2.0);
	taper/= damp_time;

        if(g_inputParam[den_smooth]>0.0){
          d->Vc[RHO][k][j][i] -= (d->Vc[RHO][k][j][i] - d_zero)*taper*dt;
//          d->Vc[VX3][k][j][i] -= (d->Vc[VX3][k][j][i] - vphi_zero)*taper*dt;

#if EOS == IDEAL
          d->Vc[PRS][k][j][i] -= (d->Vc[PRS][k][j][i] - p_zero)*taper*dt;
#endif
        }

      } else if(r >= dampout){
	
	taper = pow( (r - dampout)/(g_inputParam[rout] - dampout), 2.0);
	taper/= damp_time;
	
      } else taper = 0.0;
      
      d->Vc[VX1][k][j][i] -= (d->Vc[VX1][k][j][i])*taper*dt;
      d->Vc[VX2][k][j][i] -= (d->Vc[VX2][k][j][i])*taper*dt;

//apply density floor
//     if(d->Vc[RHO][k][j][i] < min_rho) d->Vc[RHO][k][j][i] = min_rho;

    }
  }
}

/* ************************************************************** */
void TotalMass (const Data *d, Grid *grid)
{
  int   i, j, k;
  double *dV1, *dV2, *dV3;
  double dV, mass, gmass, mc;
      
  dV1 = grid[IDIR].dV; 
  dV2 = grid[JDIR].dV; 
  dV3 = grid[KDIR].dV;

  mass = 0.0;
  DOM_LOOP(k,j,i){
    dV    = dV1[i]*dV2[j]*dV3[k];
    mass += dV*d->Vc[RHO][k][j][i];
  }
                        
#ifdef PARALLEL
  gmass = 0.;
  MPI_Allreduce (&mass, &gmass, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
  mass = 2.0*gmass; //account for lower half disk not explicitly simulated
#endif
  print1 ("> Disk mass (t=0) = %.16lf \n", mass);
}

#if (BODY_FORCE & VECTOR)
/* ************************************************************************ */
void BodyForceVector(double *v, double *g, double x1, double x2, double x3)
/*
 *
 *
 *
 *************************************************************************** */
{
  g[IDIR] = 0.0;
  g[JDIR] = 0.0;
  g[KDIR] = 0.0;
}
#endif

#if (BODY_FORCE & POTENTIAL)
/* ************************************************************************ */
double BodyForcePotential(double x1, double x2, double x3)
/*
 *
 *
 *
 *************************************************************************** */
{
  double d, R, r, z, th, x, y, phiplanet, rsm;
  double xp, yp, t, phi, mp, temp;

  #if GEOMETRY == POLAR
   R  = x1;
   #if DIMENSIONS == 2
    z  = 0.0;
    r  = R;
    th = 0.5*CONST_PI;
   #else
    z  = x3;
    r  = sqrt(R*R + z*z);
    th = atan2(R,z);
   #endif
   x  = R*cos(x2);
   y  = R*sin(x2);
  #elif (GEOMETRY == SPHERICAL)
   r  = x1;
   th = x2;
   R = r*sin(th);
   z = r*cos(th);
   x = r*sin(th)*cos(x3);
   y = r*sin(th)*sin(x3);
  #endif

/* ---------------------------------------------
             planet position
   --------------------------------------------- */

   double omega_p;
   t = g_time; 
   /* if (t > 0.0) t += g_dt;  */

  if(g_intStage == 2){
#if TIME_STEPPING == RK2 
    t += 0.5*g_dt;
#elif TIME_STEPPING == RK3
    t += 0.25*g_dt;
#endif
  }
  if(g_intStage == 3) t += 1.0/3.0*g_dt;

   omega_p = omega_k(g_unitLength);
   xp =  g_unitLength*cos(omega_p*t + CONST_PI);
   yp =  g_unitLength*sin(omega_p*t + CONST_PI);

   rsm = g_unitLength*pow(g_inputParam[mplanet]/mstar/3.0, 1.0/3.0);
   rsm *= g_inputParam[soft]; 
   d = sqrt((x-xp)*(x-xp) + (y-yp)*(y-yp) + z*z + rsm*rsm);

   if(g_time < planeton){
      mp = 0.0;
   } else if( (g_time >= planeton) && (g_time < planeton_switchon)){
     temp = sin(0.5*CONST_PI*(g_time - planeton)/switchon);
     temp*= temp;
     mp   = g_inputParam[mplanet]*temp;
   } else {
     mp   = g_inputParam[mplanet];
   }

   phiplanet = -bigG*mp/d;

//planet indirect potential. rp=g_unitLength. no z-contrib (planet at midplane)
   phiplanet += bigG*mp*(x*xp + y*yp)/pow(g_unitLength, 3.0);

   phi = grav_pot3_cylin(R, z) + phiplanet;
  return phi;
}
#endif

#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define RNMX (1.0-DBL_EPSILON)

/*! \fn double ran2(long int *idum)
 *  *  \brief  Extracted from the Numerical Recipes in C (version 2) code.  
 *   *   Modified to use doubles instead of floats. - T. A. Gardiner - Aug. 12, 2003
 *    *   
 *     *
 *      * Long period (> 2 x 10^{18}) random number generator of L'Ecuyer
 *       * with Bays-Durham shuffle and added safeguards.  Returns a uniform
 *        * random deviate between 0.0 and 1.0 (exclusive of the endpoint
 *         * values).  Call with idum = a negative integer to initialize;
 *          * thereafter, do not alter idum between successive deviates in a
 *           * sequence.  RNMX should appriximate the largest floating point value
 *            * that is less than 1. 
 *             */

double ran2(long int *idum)
{
  int j;
  long int k;
  static long int idum2=123456789;
  static long int iy=0;
  static long int iv[NTAB];
  double temp;

  if (*idum <= 0) { /* Initialize */
    if (-(*idum) < 1) *idum=1; /* Be sure to prevent idum = 0 */
    else *idum = -(*idum);
    idum2=(*idum);
    for (j=NTAB+7;j>=0;j--) { /* Load the shuffle table (after 8 warm-ups) */
      k=(*idum)/IQ1;
      *idum=IA1*(*idum-k*IQ1)-k*IR1;
      if (*idum < 0) *idum += IM1;
      if (j < NTAB) iv[j] = *idum;
    }
    iy=iv[0];
  }
  k=(*idum)/IQ1;                 /* Start here when not initializing */
  *idum=IA1*(*idum-k*IQ1)-k*IR1; /* Compute idum=(IA1*idum) % IM1 without */
  if (*idum < 0) *idum += IM1;   /* overflows by Schrage's method */
  k=idum2/IQ2;
  idum2=IA2*(idum2-k*IQ2)-k*IR2; /* Compute idum2=(IA2*idum) % IM2 likewise */
  if (idum2 < 0) idum2 += IM2;
  j=(int)(iy/NDIV);              /* Will be in the range 0...NTAB-1 */
  iy=iv[j]-idum2;                /* Here idum is shuffled, idum and idum2 */
  iv[j] = *idum;                 /* are combined to generate output */
  if (iy < 1) iy += IMM1;
  if ((temp=AM*iy) > RNMX) return RNMX; /* No endpoint values */
  else return temp;
}
#undef IM1
#undef IM2
#undef AM
#undef IMM1
#undef IA1
#undef IA2
#undef IQ1
#undef IQ2
#undef IR1
#undef IR2
#undef NTAB
#undef NDIV
#undef RNMX


