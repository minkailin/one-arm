#include "mp.h"

extern real ScalingFactor;

/* Surface density */
real Sigma(r)
     real r;
{
  real density, rmax; 
 
  rmax = RMAX;
  if(RDEAD2 < RMAX) rmax = RDEAD2; //if outer dead edge within domain, then Q is Q at outer edge of dead zone

  SIGMA0 = sqrt(csq(rmax))*pow(rmax, -1.5);//cs*omega
  SIGMA0 /=QOUT*PI*pow(rmax, -SIGMASLOPE)*hole(rmax);//divive by pi*G*...

  density = SIGMA0*pow(r,-SIGMASLOPE);

  return density*hole(r);
}

double RadialTaper(r)
       double r;
{
       double res, DRMIN, DRMAX, R1, R2;
        
R1    = RMIN*DAMPRADIN;//0.5*(RMIN+RDEAD1);
R2    = RMAX*DAMPRADOUT;//0.5*(RMAX+RDEAD2);
DRMIN = R1*DAMPRADIN;
DRMAX = R2*DAMPRADOUT;

if( (r>=DRMIN) && (r<=DRMAX) ){
res = 1.0;
} else if( (r < DRMIN) && (r >= R1) ){
res = (R1 - r)/(R1 - DRMIN);
res = res*res;
}else if( (r > DRMAX) && (r<=R2) ){
res = (R2 - r)/(R2 - DRMAX);
res = res*res;
} else if ( (r < R1) || (r>R2) ){
res = 0.0;
}
return res;
}

double hole(r)
real r;
{
  real dr1, fr1, dr2, fr2, eps, Hin; 

  Hin = bigH(RMIN);
  eps = FIN;

  dr1 = DNU*bigH(RDEAD1);
  fr1 = (1.0 - eps)*(1.0 + tanh((r-RDEAD1)/dr1))/2.0 + eps;

  dr2 = DNU*bigH(RDEAD2);
  fr2 = (1.0 - eps)*(1.0 - tanh((r-RDEAD2)/dr2))/2.0 + eps;


  return fr1*fr2;//*(1.0 - sqrt(RMIN/(r+Hin))); 
}

double dloghole(r)
real r;
{
  real eps; 
  real dr1, dfr1, dr2, dfr2, sech1, sech2;

  eps = FIN;

  dr1 = DNU*bigH(RDEAD1);
  sech1 = 1.0/cosh( (r-RDEAD1)/dr1);
  dfr1 = (1.0 - eps)*pow(sech1, 2.0);
  dfr1/= 2.0*dr1; 

  dr2 = DNU*bigH(RDEAD2);
  sech2 =  1.0/cosh( (r-RDEAD2)/dr2);
  dfr2 = -(1.0 - eps)*pow(sech2, 2.0);
  dfr2/= 2.0*dr2;

  return (dfr1 + dfr2)/hole(r); 
}

double d2loghole(r)
real r;
{
  real eps;
  real dr1, dr2, d2fr1, d2fr2, sech1, sech2, arg1, arg2;

  eps = FIN;

  dr1 = DNU*bigH(RDEAD1);
  arg1 = (r-RDEAD1)/dr1; 
  sech1 = 1.0/cosh(arg1);
  d2fr1 = (1.0 - eps)*pow(sech1, 2.0)*tanh(arg1);
  d2fr1/= -dr1*dr1;

  dr2 = DNU*bigH(RDEAD2);
  arg2 = (r-RDEAD2)/dr2;
  sech2 = 1.0/cosh(arg2);
  d2fr2 = -(1.0 - eps)*pow(sech2, 2.0)*tanh(arg2);
  d2fr2/= -dr2*dr2;
 
  return (d2fr1+d2fr2)/hole(r) - pow(dloghole(r), 2.0);
}

double bigH(r)
real r;
{
 real omegak;

 omegak = pow(r, -1.5); 
 return sqrt(csq(r))/omegak;
}

double csq(r)
real r;
{
    return ASPECTRATIO*ASPECTRATIO*pow(r, 2.0*FLARINGINDEX - 1.0); 
}

double dlogcsq(r)
real r;
{
    return (2.0*FLARINGINDEX - 1.0)/r;
}

double d2logcsq(r)
real r;
{
    return -dlogcsq(r)/r;
}

double vphi0(r)
real r;
{
   real v; 
   v = 1.0/r; //contribution from central star, G=M=1
   v+= r*csq(r)*( dlogcsq(r) - SIGMASLOPE/r + dloghole(r) );
   return sqrt(v); 
}

double omega0(r)
real r;
{
 return vphi0(r)/r; 
}

double domega0(r)
real r;
{
  real om, dom, omksq; 
  
  om = omega0(r); 
  omksq = 1.0/pow(r,3.0); 

  dom = -2.0*omksq;
  dom+=csq(r)*( dlogcsq(r)*(dlogcsq(r) - SIGMASLOPE/r + dloghole(r) ) + d2logcsq(r) + SIGMASLOPE/(r*r) + d2loghole(r) );
  dom-= om*om;
  dom/= 2.0*r*om;
  return dom;
}

double vr0(r)
real r;
{
 real numax, sigmax, v;
/*
 numax = VISCOSITY;
 sigmax= Sigma(RMAX); 
 v = -numax*sigmax*pow(RMAX, 3.0)*domega0(RMAX); 
 v/= RMIN*RMIN*omega0(RMIN) - RMAX*RMAX*omega0(RMAX);
 v/= r*Sigma(r);
 return v; 
*/
 return 0.0;
}


void FillSigma() {
  int i;
  for (i = 0; i < NRAD; i++) {
    SigmaMed[i] = Sigma(Rmed[i]);
    SigmaInf[i] = Sigma(Rinf[i]);
  }
}

void RefillSigma (Surfdens)
     PolarGrid *Surfdens;
{
  int i, j, nr, ns, l;
  real *field;
  real moy;
  nr = Surfdens->Nrad;
  ns = Surfdens->Nsec;
  field = Surfdens->Field;
  for (i = 0; i < nr; i++) {
    moy = 0.0;
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      moy += field[l];
    }
    moy /= (real)ns;
    SigmaMed[i] = moy;
  }
  SigmaInf[0] = SigmaMed[0];
  for (i = 1; i < nr; i++) {
    SigmaInf[i] = (SigmaMed[i-1]*(Rmed[i]-Rinf[i])+\
		   SigmaMed[i]*(Rinf[i]-Rmed[i-1]))/\
      (Rmed[i]-Rmed[i-1]);
  }
}

/* Thermal energy */
real Energy(r)
     real r;
{
  real energy0;
  if (ADIABATICINDEX == 1.0) {
    fprintf (stderr, "The adiabatic index must differ from unity to initialize the gas internal energy. I must exit.\n");
    prs_exit (1);
  }
  else
    /*energy0 = R/MU/(ADIABATICINDEX-1.0)*SIGMA0*pow(ASPECTRATIO,2.0)*pow(r,-SIGMASLOPE-1.0+2.0*FLARINGINDEX);*/
    energy0 = R/MU/(ADIABATICINDEX-1.0)*Sigma(r)*pow(ASPECTRATIO,2.0)*pow(r,-1.0+2.0*FLARINGINDEX);
  return energy0;
}

void FillEnergy() {
  int i;
  for (i = 0; i < NRAD; i++)
    EnergyMed[i] = Energy(Rmed[i]);
}


void RefillEnergy (energy)
     PolarGrid *energy;
{
  int i, j, nr, ns, l;
  real *field;
  real moy;
  nr = energy->Nrad;
  ns = energy->Nsec;
  field = energy->Field;
  for (i = 0; i < nr; i++) {
    moy = 0.0;
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      moy += field[l];
    }
    moy /= (real)ns;
    EnergyMed[i] = moy;
  }
}

/* Cooling time */
real CoolingTime(r)
     real r;
{
  real eps = 0.01, rp=10.0; 
  real omega, omegap;
  real ct0;
 
  ct0 = COOLINGTIME0*pow(r,1.5);
  return ct0;
}

void FillCoolingTime() {
  int i;
  for (i = 0; i < NRAD; i++)
    CoolingTimeMed[i] = CoolingTime(Rmed[i]);
}

/* Heating source term */
real Qplusinit(r)
     real r;
{
  real qp0, viscosity;
  real romegak_prime;
  
  viscosity = FViscosity(r);///Sigma(r);
  romegak_prime = r*domega0(r); 
  
  qp0 = (1.0/2.0)*viscosity*Sigma(r)*pow(romegak_prime, 2.0);
  return qp0;
}

void FillQplus() {
  int i;
  for (i = 0; i < NRAD; i++)
    QplusMed[i] = Qplusinit(Rmed[i]);
}
