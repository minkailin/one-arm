#include "mp.h"

#define CFLSECURITY 0.5		/* Maximum fraction of zone size */
				/* swept in one timestep */

#define CVNR 1.41       	/* Shocks are spread over CVNR zones:       */
                                /* von Neumann-Richtmyer viscosity constant */
				/* Beware of misprint in Stone and Norman's */
				/* paper : use C2^2 instead of C2           */

static PolarGrid *TemperInt;
static PolarGrid *VradNew, *VradInt;
static PolarGrid *VthetaNew, *VthetaInt;
static PolarGrid *EnergyNew, *EnergyInt;
static real timeCRASH;  
extern boolean Corotating;
extern boolean Adiabatic;
extern boolean SelfGravity, SGZeroMode;
extern boolean ZMPlus;
real PhysicalTime=0.0, OmegaFrame, PhysicalTimeInitial;
int FirstGasStepFLAG=1;
static int AlreadyCrashed = 0, GasTimeStepsCFL;
extern boolean FastTransport, IsDisk;
extern boolean PERTURBED;
Pair DiskOnPrimaryAcceleration;

boolean DetectCrash (array)
     PolarGrid *array;
{
  int i, j, l, nr, ns;
  real *ptr;
  boolean bool = NO;
  nr = array->Nrad;
  ns = array->Nsec;
  ptr= array->Field;
#pragma omp parallel for private(j,l) shared(bool)
  for (i = 0; i < nr; i++) {
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      if (ptr[l] < 0.0) {
	bool = YES;
      }
    }
  }
  return bool;
}
 
void DensityFloor(Rho)
	PolarGrid *Rho;
{
	int i, j, l, nr, ns;
	real *rho;
	nr = Rho -> Nrad;
	ns = Rho -> Nsec;
	rho = Rho -> Field;
#pragma omp parallel for private(j,l) 
  for (i = 0; i < nr; i++) {
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      if (rho[l] < DENFLOOR*SIGMA0)
        rho[l] = DENFLOOR*SIGMA0;
    }
  }
}




void AzimuthalAverage (Vrad, Vtheta, Rho, Energy)
     PolarGrid *Vrad, *Vtheta, *Rho, *Energy;
{
  int i, j, l, nr, ns;
  real *vrad, *vtheta, *dens, *energ;
  real vrad0, vtheta0, dens0, energ0;
  vrad = Vrad->Field;
  vtheta = Vtheta->Field;
  dens = Rho->Field;
  energ = Energy->Field;
  nr = Rho->Nrad;
  ns = Rho->Nsec;

  for (i = 0; i < nr; i++) {

       vrad0   =  0.0;
       vtheta0 =  0.0;
       dens0   =  0.0;
       energ0  =  0.0;

         for(j = 0; j < ns; j++){
            l = i*ns + j;
            vrad0  += vrad[l];
            vtheta0+= vtheta[l];
            dens0  += dens[l];
            if(Adiabatic) energ0+=energ[l];
         }
             vrad0 /= (real)ns;
             vtheta0/= (real)ns;
             dens0  /= (real)ns;
           if(Adiabatic) energ0/=(real)ns;

    for (j = 0; j < ns; j++) {
        l = i*ns + j;
        vrad[l]   = vrad0;
        vtheta[l] = vtheta0;
        dens[l]   = dens0;
        if (Adiabatic)
          energ[l]  = energ0;
      }
    }
}

void Filter (Vrad, Vtheta, Rho, Energy, dt)
     PolarGrid *Vrad, *Vtheta, *Rho, *Energy;
     real dt;
{
  int i, j, l, nr, ns;
  real *vrad, *vtheta, *dens, *energ;
  real dens0, dens1_real, dens1_imag;
  real vr_0, vr1_real, vr1_imag;
  real vt0, vt1_real, vt1_imag;
  real energ0, energ1_real, energ1_imag;
  real phi, cosphi, sinphi;
  vrad = Vrad->Field;
  vtheta = Vtheta->Field;
  dens = Rho->Field;
  energ = Energy->Field;
  nr = Rho->Nrad;
  ns = Rho->Nsec;

  for (i = 0; i < nr; i++) {
      
    dens0      = 0.0;
    dens1_real = 0.0;
    dens1_imag = 0.0;

    vr_0      = 0.0;
    vr1_real = 0.0;
    vr1_imag = 0.0;

    vt0      = 0.0;
    vt1_real = 0.0;
    vt1_imag = 0.0;
    
    for(j = 0; j < ns; j++){
      l = i*ns + j;
	   phi = 2.0*PI*(real)j/ns;
	   cosphi = cos(phi);
	   sinphi = sin(phi);
	   
	   dens0      += dens[l];
	   dens1_real += dens[l]*cosphi;
           dens1_imag +=-dens[l]*sinphi;
	   
	   vr_0      += vrad[l];
	   vr1_real += vrad[l]*cosphi;
           vr1_imag +=-vrad[l]*sinphi;
	   
	   vt0      += vtheta[l];
	   vt1_real += vtheta[l]*cosphi;
           vt1_imag +=-vtheta[l]*sinphi;
	   
	   
	   if(Adiabatic) {
	     energ0+=energ[l];
	     energ1_real += energ[l]*cosphi;
	     energ1_imag +=-energ[l]*sinphi;
	   }
    }
    dens0      /= (real)ns;
    dens1_real /= (real)ns;
    dens1_imag /= (real)ns;
    
    vr_0      /= (real)ns;
    vr1_real /= (real)ns;
    vr1_imag /= (real)ns;
    
    vt0      /= (real)ns;
    vt1_real /= (real)ns;
    vt1_imag /= (real)ns;
    
    if(Adiabatic) {
      energ0/=(real)ns;
      energ1_real /= (real)ns;
      energ1_imag /= (real)ns;
    }
    
    for (j = 0; j < ns; j++) {
      l = i*ns + j;
      phi = 2.0*PI*(real)j/ns;
      cosphi = cos(phi);
      sinphi = sin(phi);
      
      dens[l] = SigmaMed[i] + 2.0*(dens1_real*cosphi - dens1_imag*sinphi);  
      
      vrad[l] = vr0(Rmed[i]) + 2.0*(vr1_real*cosphi - vr1_imag*sinphi); 
      
      vtheta[l] = sqrt(pow(vphi0(Rmed[i]),2.0) - \
                          Rmed[i]*GLOBAL_AxiSGAccr0[i+IMIN])  + 2.0*(vt1_real*cosphi - vt1_imag*sinphi); 
      
      if (Adiabatic) energ[l]  =  EnergyMed[i] + 2.0*(energ1_real*cosphi - energ1_imag*sinphi);
      
    }
  }
}


void DensityPerturbation(Vtheta, Vrad, Rho)//introduce m=1 pert without modifying total ang mom 
        PolarGrid *Vtheta, *Vrad, *Rho;
{
  int i, j, l, nr, ns, seed;
        real m;
        real rc, omega_p, nu, Q, Qsq, kappa, Omega, k_T, phi, k, test;
	real shat, vrhat, vthat; 
        real *vtheta, *vrad, *dens;
        nr = Vtheta -> Nrad;
        ns = Vtheta -> Nsec;
        vtheta = Vtheta -> Field;
        vrad   = Vrad -> Field;
        dens = Rho  -> Field; 

        srand((unsigned) (CPU_Rank+1)*time(0));
	
	m=1;
	rc = RMAX;
	omega_p = pow(rc, -1.5);

  for (i = 0; i < nr; i++) {
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      phi = 2.0*PI*(double)j/(double)ns;

      Omega = vtheta[l]/Rmed[i];
      kappa = Omega; 
 
      Q = sqrt(csq(Rmed[i]))*kappa/(PI*dens[l]);
      k_T = kappa*kappa/(2.0*PI*dens[l]);

      nu = m*(omega_p - Omega)/kappa; 

      Qsq = Q*Q; 

      test = Qsq*(1.0 - nu*nu);
      
      if(test < 1.0){
	k = (2.0/Qsq)*(1.0 - sqrt(1.0 - test)); //magnitude of dimensionless wavenum (long wave)
	k*=-k_T; //attach units and choose sign for trailing
  
	shat  = PERTMAX*dens[l]; 
	vrhat = -m*(omega_p - Omega)*shat/(k*dens[l]);
	vthat = (kappa*kappa)*shat/(2.0*Omega*dens[l]*k); 

	dens[l]  += shat*cos(-m*phi +k*Rmed[i]);
	vrad[l]  += vrhat*cos(-m*phi+k*Rmed[i]);
	vtheta[l]+= vthat*sin(-m*phi+k*Rmed[i]);
      }

    }
  }
  
}

void DensityScale(Rho, dt)
        PolarGrid *Rho;
	real      dt;
{
        int i, j, l, nr, ns;
        real *rho;
        nr = Rho -> Nrad;
        ns = Rho -> Nsec;
        rho = Rho -> Field;
#pragma omp parallel for private(j,l) 
  for (i = 0; i < nr; i++) {
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
        rho[l] += (SIGSCALE - 1.0)*SigmaMed[i]*dt/(NTOT*DT - TSCALE*2.0*M_PI);
    }
  }
}

void FillPolar1DArrays ()
{
  FILE *input, *output;
  int i,ii;
  real drrsep;
  float temporary;
  char InputName[256], OutputName[256];
  drrsep = (RMAX-RMIN)/(real)GLOBALNRAD;
  sprintf (InputName, "%s%s", OUTPUTDIR, "radii.dat");
  sprintf (OutputName, "%s%s", OUTPUTDIR, "used_rad.dat");
  input = fopen (InputName, "r");
  if (input == NULL) {
    mastererr ("Warning : no `radii.dat' file found. Using default.\n");
    if (LogGrid == YES) {
      for (i = 0; i <= GLOBALNRAD; i++) {
	Radii[i] = RMIN*exp((real)i/(real)GLOBALNRAD*log(RMAX/RMIN));
      }
    } else {
      for (i = 0; i <= GLOBALNRAD; i++) {
	Radii[i] = RMIN+drrsep*(real)(i);
      }
    }
  } else {
    mastererr ("Reading 'radii.dat' file.\n");
    for (i = 0; i <= GLOBALNRAD; i++) {
      fscanf (input, "%f", &temporary);
      Radii[i] = (real)temporary;
    }
  }
  for (i = 0; i < GLOBALNRAD; i++) {
    GlobalRmed[i] = 2.0/3.0*(Radii[i+1]*Radii[i+1]*Radii[i+1]-Radii[i]*Radii[i]*Radii[i]);
    GlobalRmed[i] = GlobalRmed[i] / (Radii[i+1]*Radii[i+1]-Radii[i]*Radii[i]);
  }
  for (i = 0; i < NRAD; i++) {
    ii = i+IMIN;
    Rinf[i] = Radii[ii];
    Rsup[i] = Radii[ii+1];
    Rmed[i] = 2.0/3.0*(Rsup[i]*Rsup[i]*Rsup[i]-Rinf[i]*Rinf[i]*Rinf[i]);
    Rmed[i] = Rmed[i] / (Rsup[i]*Rsup[i]-Rinf[i]*Rinf[i]);
    Surf[i] = PI*(Rsup[i]*Rsup[i]-Rinf[i]*Rinf[i])/(real)NSEC;
    InvRmed[i] = 1.0/Rmed[i];
    InvSurf[i] = 1.0/Surf[i];
    InvDiffRsup[i] = 1.0/(Rsup[i]-Rinf[i]);
    InvRinf[i] = 1.0/Rinf[i];
  }
  Rinf[NRAD]=Radii[NRAD+IMIN];
  for (i = 1; i < NRAD; i++) {
    InvDiffRmed[i] = 1.0/(Rmed[i]-Rmed[i-1]);
  }
  if (CPU_Master) {
    output = fopen (OutputName, "w");
    if (output == NULL) {
      mastererr ("Can't write %s.\nProgram stopped.\n", OutputName);
      prs_exit (1);
    }
    for (i = 0; i <= GLOBALNRAD; i++) {
      fprintf (output, "%.18g\n", Radii[i]);
    }
    fclose (output);
  }
  if (input != NULL) fclose (input);
}

void InitEuler (Vr, Vt, Rho, Energy)
     PolarGrid *Vr, *Vt, *Rho, *Energy;
{
  InitTransport ();
  InitViscosity ();
  RhoStar      = CreatePolarGrid(NRAD, NSEC, "RhoStar");
  RhoInt       = CreatePolarGrid(NRAD, NSEC, "RhoInt");
  VradNew      = CreatePolarGrid(NRAD, NSEC, "VradNew");
  VradInt      = CreatePolarGrid(NRAD, NSEC, "VradInt");
  VthetaNew    = CreatePolarGrid(NRAD, NSEC, "VthetaNew");
  VthetaInt    = CreatePolarGrid(NRAD, NSEC, "VthetaInt");
  EnergyNew    = CreatePolarGrid(NRAD, NSEC, "EnergyNew");
  EnergyInt    = CreatePolarGrid(NRAD, NSEC, "EnergyInt");
  TemperInt    = CreatePolarGrid(NRAD, NSEC, "TemperInt");
  Potential    = CreatePolarGrid(NRAD, NSEC, "Potential");
  Pressure     = CreatePolarGrid(NRAD, NSEC, "Pressure");
  SoundSpeed   = CreatePolarGrid(NRAD, NSEC, "SoundSpeed");
  Temperature  = CreatePolarGrid(NRAD, NSEC, "Temperature");
  Qplus        = CreatePolarGrid(NRAD, NSEC, "Qplus");
  InitComputeAccel ();
  /* Rho and Energy are already initialized: cf main.c */
  ComputeSoundSpeed (Rho, Energy);
  ComputePressureField0 (Rho, Energy);
  ComputeTemperatureField (Rho, Energy);
  InitGasVelocities (Vr, Vt);
}

real min2 (a,b)
     real a,b;
{
  if (b < a) return b;
  return a;
}

real max2 (a,b)
     real a,b;
{
  if (b > a) return b;
  return a;
}


void ActualiseGas (array, newarray)
     PolarGrid *array, *newarray;
{
  int i,j,l,ns,nr;
  real *old, *new;
  nr = array->Nrad;
  ns = array->Nsec;
  old= array->Field;
  new= newarray->Field;
#pragma omp parallel for private(j,l)
  for (i = 0; i <= nr; i++) {
    for (j = 0; j < ns; j++) {
      l = j+ns*i;
      old[l] = new[l];
    }
  }
}


void AlgoGas (force, Rho, Vrad, Vtheta, Energy, Label, sys)
     Force *force;
     PolarGrid *Rho, *Vrad, *Vtheta, *Energy, *Label;
     PlanetarySystem *sys;
{
  extern boolean ForcedCircular;
  int i, j, l, nr, ns;
  real *rho;
  real dt, dtemp=0.0;
  real OmegaNew, domega, mini;
  int gastimestepcfl;
  boolean Crashed=NO;
  FirstGasStepFLAG=1;
  gastimestepcfl = 1;
  nr = Vrad->Nrad;
  ns = Vrad->Nsec;
  rho = Rho->Field;
  /*Perform density perturbation*/
  if((PhysicalTime > PERTTIME*2.0*M_PI)){
        /*printf("here1 %d %d \n", CPU_Rank, PERTURBED);*/
        if(PERTURBED == NO){
        /*printf("here2 %d \n", CPU_Rank);*/
        AzimuthalAverage(Rho, Vrad, Vtheta, Energy);
	DensityPerturbation(Vtheta, Vrad, Rho);
        PERTURBED = YES;
        }
  }
  if (Adiabatic) {
    ComputeSoundSpeed (Rho, Energy);
    /* it is necesary to update computation of soundspeed if one uses
       alphaviscosity in FViscosity. It is not necesary in locally
       isothermal runs since cs is constant.  It is computed here for
       the needs of ConditionCFL. */
  }
  if (IsDisk == YES) {
    CommunicateBoundaries (Rho, Vrad, Vtheta, Energy, Label);
    if (SloppyCFL == YES)
      gastimestepcfl = ConditionCFL (Vrad, Vtheta, DT-dtemp, sys);
  }
  MPI_Allreduce (&gastimestepcfl, &GasTimeStepsCFL, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  dt = DT / (real)GasTimeStepsCFL;
  while (dtemp < 0.999999999*DT) {
    if(PhysicalTime < PLANET_ON*2.0*M_PI){
	MassTaper = 0.0;
    } else {
    MassTaper = (PhysicalTime-PLANET_ON*2.0*M_PI)/(MASSTAPER*2.0*M_PI);
    MassTaper = (MassTaper > 1.0 ? 1.0 : pow(sin(MassTaper*M_PI/2.0),2.0));
    }

 
    if (IsDisk == YES) {
      CommunicateBoundaries (Rho, Vrad, Vtheta, Energy, Label);
      if (SloppyCFL == NO) {
	gastimestepcfl = 1;
	gastimestepcfl = ConditionCFL (Vrad, Vtheta, DT-dtemp, sys);
	MPI_Allreduce (&gastimestepcfl, &GasTimeStepsCFL, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
	dt = (DT-dtemp)/(real)GasTimeStepsCFL;
      }
      AccreteOntoPlanets (Rho, Vrad, Vtheta, dt, sys);
    }
    dtemp += dt;
    DiskOnPrimaryAcceleration.x = 0.0;
    DiskOnPrimaryAcceleration.y = 0.0;
    if (Corotating == YES) GetPsysInfo (sys, MARK);
    if (IsDisk == YES) {
      /* Indirect term star's potential computed here */
      DiskOnPrimaryAcceleration   = ComputeAccel (force, Rho, 0.0, 0.0, 0.0, 0.0);
      /* Gravitational potential from star and planet(s) is computed
	 and stored here*/
      FillForcesArrays (sys, Rho, Energy);
      /* Planets' velocities are updated here from gravitationnal
	 interaction with disk */
      AdvanceSystemFromDisk (force, Rho, Energy, sys, dt);
    }
    /* Planets' positions and velocities are updated from
       gravitational interaction with star and other planets */
    AdvanceSystemRK5 (sys, dt);
    /* Below we correct vtheta, planet's position and velocities if we
       work in a frame non-centered on the star */
    if (Corotating == YES) {
      OmegaNew = GetPsysInfo(sys, GET) / dt;
      domega = OmegaNew-OmegaFrame;
      if (IsDisk == YES)
	CorrectVtheta (Vtheta, domega);
      OmegaFrame = OmegaNew;
    }
    RotatePsys (sys, OmegaFrame*dt);
    /* Now we update gas */
    if (IsDisk == YES) {
      ApplyBoundaryCondition (Vrad, Vtheta, Rho, Energy, dt);

//      DensityFloor(Rho);

      MPI_Barrier(MPI_COMM_WORLD);
      Crashed = DetectCrash (Rho);    /* test for negative density values */
      if(Crashed == YES) printf("density crash, post BC: %i \n", CPU_Rank);
      Crashed = DetectCrash (Energy);  /* test for negative energy values */
      if (Crashed == YES) {
	if (AlreadyCrashed == 0) {
	  timeCRASH=PhysicalTime;   /* if it appears to be the first crash */
	  fprintf (stdout,"\nCrash! at time %.12g\n", timeCRASH);
	  WriteDiskPolar (Rho, 999);    /* We write the HD arrays */
	  WriteDiskPolar (Vrad, 999);   /* in order to keep a track */
	  WriteDiskPolar (Vtheta, 999); /* of what happened */
	  WriteDiskPolar (Temperature, 999);
	}
	AlreadyCrashed++;
	masterprint ("c");
      } else {
	masterprint (">");
      }
      fflush (stdout);
      if (ZMPlus)
	compute_anisotropic_pressurecoeff (sys);
      /*    ADD SOURCE MASS */
      if(PhysicalTime > TSCALE*2.0*M_PI) DensityScale(Rho, dt);
     if(PhysicalTime < PLANET_ON*2.0*M_PI){
      ComputePressureField0 (Rho, Energy);
    	} else {
		ComputePressureField (Rho, Energy, sys);
    	}

     
//     Filter (Vrad, Vtheta, Rho, Energy, dt); 
     SubStep1 (Vrad, Vtheta, Rho, dt);

//     Filter (Vrad, Vtheta, Rho, Energy, dt);
     SubStep2 (Rho, Energy, dt);

     ActualiseGas (Vrad, VradNew);
     ActualiseGas (Vtheta, VthetaNew);
     ApplyBoundaryCondition (Vrad, Vtheta, Rho, Energy, dt);
     
      
      if (Adiabatic) {
	ComputeViscousTerms (Vrad, Vtheta, Rho);
	SubStep3 (Rho, dt);
	ActualiseGas (Energy, EnergyNew);
      }

//      Filter (Vrad, Vtheta, Rho, Energy, dt); 
      Transport (Rho, Vrad, Vtheta, Energy, Label, dt);

      ApplyBoundaryCondition (Vrad, Vtheta, Rho, Energy, dt);

      ComputeTemperatureField (Rho, Energy);
      mdcp = CircumPlanetaryMass (Rho, sys);
      exces_mdcp = mdcp - mdcp0;
     }
    PhysicalTime += dt;
  }
  masterprint ("\n");
}

void SubStep1 (Vrad, Vtheta, Rho, dt)
     PolarGrid *Vrad, *Vtheta, *Rho;
     real dt;
{
  int i, j, l, lim, ljm, ljp, nr, ns;
  boolean selfgravityupdate;
  extern boolean Evanescent;
  real *vrad, *vtheta, *rho;
  real *Pot, *Press;
  real *vradint, *vthetaint;
  real gradp, gradphi, vt2, dxtheta;
  real invdxtheta;
  real supp_torque=0.0;		/* for imposed disk drift */
  nr = Vrad->Nrad;
  ns = Vrad->Nsec;
  rho = Rho->Field;
  vrad = Vrad->Field;
  vtheta = Vtheta->Field;
  vradint = VradInt->Field;
  vthetaint = VthetaInt->Field;
  Pot = Potential->Field;
  Press = Pressure->Field;
  /* In this substep we take into account  */
  /* the source part of Euler equations  */
  /* We evolve velocities with pressure gradients */
  /* gravitational forces and curvature terms */
#pragma omp parallel private(j,l,lim,ljm,ljp,dxtheta,vt2,gradp,gradphi,invdxtheta,supp_torque)
  {
#pragma omp for nowait
    for (i = 1; i < nr; i++) {
      for (j = 0; j < ns; j++) {
	l = j+i*ns;
	lim = l-ns;
	ljp = l+1;
	if (j == ns-1) ljp = i*ns;
	gradp = (Press[l]-Press[lim])*2.0/(rho[l]+rho[lim])*InvDiffRmed[i];
	gradphi = (Pot[l]-Pot[lim])*InvDiffRmed[i];
	vt2 = vtheta[l]+vtheta[ljp]+vtheta[lim]+vtheta[ljp-ns];
	vt2 = vt2/4.0+Rinf[i]*OmegaFrame;
	vt2 = vt2*vt2;
	vradint[l] = vrad[l]+dt*(-gradp-gradphi+vt2*InvRinf[i]);
      }
    }
#pragma omp for
    for (i = 0; i < nr; i++) {
      supp_torque = IMPOSEDDISKDRIFT*0.5*pow(Rmed[i],-2.5+SIGMASLOPE);
      dxtheta = 2.0*PI/(real)ns*Rmed[i];
      invdxtheta = 1.0/dxtheta;
      for (j = 0; j < ns; j++) {
	l = j+i*ns;
	ljm = l-1;
	if (j == 0) ljm = i*ns+ns-1;
	gradp = (Press[l]-Press[ljm])*2.0/(rho[l]+rho[ljm])*invdxtheta;
	if ( ZMPlus )
	  gradp *= SG_aniso_coeff;
	gradphi = (Pot[l]-Pot[ljm])*invdxtheta;
	vthetaint[l] = vtheta[l]-dt*(gradp+gradphi);
	vthetaint[l] += dt*supp_torque;
      }
    }
  }
  if ( SelfGravity ) {
    selfgravityupdate = YES;
    compute_selfgravity(Rho, VradInt, VthetaInt, dt, selfgravityupdate);
  }
  ComputeViscousTerms (VradInt, VthetaInt, Rho);
  UpdateVelocitiesWithViscosity (VradInt, VthetaInt, Rho, dt);
   if ( !Evanescent ) 
    ApplySubKeplerianBoundary (VthetaInt);

}

void SubStep2 (Rho, Energy, dt)
     PolarGrid *Rho, *Energy;
     real dt;
{
  int i, j, l, lim, lip, ljm, ljp, nr, ns;
  real *vrad, *vtheta, *rho, *energy;
  real *vradnew, *vthetanew, *qt, *qr, *energyint;
  real dxtheta, invdxtheta;
  real dv;
  nr = Rho->Nrad;
  ns = Rho->Nsec;
  rho = Rho->Field;
  vrad = VradInt->Field;
  vtheta = VthetaInt->Field;
  qr = RhoInt->Field;
  vradnew = VradNew->Field;
  vthetanew = VthetaNew->Field;
  qt = TemperInt->Field;
  energy = Energy->Field;
  energyint = EnergyInt->Field;
#pragma omp parallel for private(j,dxtheta,l,lim,lip,ljm,ljp,dv)
  for (i = 0; i < nr; i++) {
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      lip = l+ns;
      ljp = l+1;
      if (j == ns-1) ljp = i*ns;
      dv = vrad[lip]-vrad[l];
      if (dv < 0.0)
        qr[l] = CVNR*CVNR*rho[l]*dv*dv;
      else 
        qr[l] = 0.0; 
      dv = vtheta[ljp]-vtheta[l];
      if (dv < 0.0)
        qt[l] = CVNR*CVNR*rho[l]*dv*dv;
      else
	qt[l] = 0.0;
    }
  }
#pragma omp parallel private(l,lim,lip,ljm,ljp,j,dxtheta,invdxtheta)
  {
#pragma omp for nowait
    for (i = 1; i < nr; i++) {
      for (j = 0; j < ns; j++) {
	l = j+i*ns;
	lim = l-ns;
	vradnew[l] = vrad[l]-dt*2.0/(rho[l]+rho[lim])*(qr[l]-qr[lim])*InvDiffRmed[i];
      }
    }
#pragma omp for
    for (i = 0; i < nr; i++) {
      dxtheta = 2.0*PI/(real)ns*Rmed[i];
      invdxtheta = 1.0/dxtheta;
      for (j = 0; j < ns; j++) {
	l = j+i*ns;
	ljm = l-1;
	if (j == 0) ljm = i*ns+ns-1;
	vthetanew[l] = vtheta[l]-dt*2.0/(rho[l]+rho[ljm])*(qt[l]-qt[ljm])*invdxtheta;
      }
    }
    /* If gas disk is adiabatic, we add artificial viscosity as a source */
    /* term for advection of thermal energy polargrid */
    if (Adiabatic) {
#pragma omp for nowait
      for (i = 0; i < nr; i++) {
	dxtheta = 2.0*PI/(real)ns*Rmed[i];
	invdxtheta = 1.0/dxtheta;
	for (j = 0; j < ns; j++) {
	  l = j+i*ns;
	  lip = l+ns;
	  ljp = l+1;
	  if (j == ns-1) ljp = i*ns;
	  energyint[l] = energy[l] -				\
	    dt*qr[l]*(vrad[lip]-vrad[l])*InvDiffRsup[i] -	\
	    dt*qt[l]*(vtheta[ljp]-vtheta[l])*invdxtheta;
	}
      }
    }
  }
}
	       
void SubStep3 (Rho, dt)
     PolarGrid *Rho;
     real dt;
{
  extern boolean Cooling;
  int i, j, l, nr, ns;
  int lip, li2p;
  real r, rip, ri2p, qpip, qpi2p;
  real *dens, *pres, *energy, *energynew;
  real *divergence, *Trr, *Trp, *Tpp, *qplus;
  real viscosity, energypred, num, den;
  nr = Rho->Nrad;
  ns = Rho->Nsec;
  dens = Rho->Field;
  pres = Pressure->Field;
  energy = EnergyInt->Field;
  energynew = EnergyNew->Field;
  divergence = DivergenceVelocity->Field;
  qplus = Qplus->Field;
  Trr = TAURR->Field;
  Trp = TAURP->Field;
  Tpp = TAUPP->Field;
  /* In this substep we take into account  */
  /* the source part of energy equation  */
  /* We evolve internal energy with  */
  /* compression/dilatation and heating terms */
#pragma omp parallel private(j,l,viscosity,energypred)
  {
#pragma omp for
    /* We calculate the heating source term Qplus from i=1 */
    for (i = 1; i < nr; i++) { /* Trp defined from i=1 */
      viscosity = FViscosity (Rmed[i]); /* Global_Bufarray holds cs */
      for (j = 0; j < ns; j++) {
	l = j+i*ns;
//        viscosity = FViscosity (Rmed[i]);
//        viscosity *= pow(ANU, 1.0 - dens[l]/Sigma(Rmed[i]));//visciosity reduction if density increases
	if (viscosity != 0.0) {
	  qplus[l] = 0.5/viscosity/dens[l]*(Trr[l]*Trr[l] +	\
					    Trp[l]*Trp[l] +	\
					    Tpp[l]*Tpp[l] );
	  qplus[l] += (2.0/9.0)*viscosity*dens[l]*divergence[l]*divergence[l];
	}
	else
	  qplus[l] = 0.0;
      }
    }
    /* We calculate the heating source term Qplus for i=0 */
    i = 0;
    r    = Rmed[i];
    rip  = Rmed[i+1];
    ri2p = Rmed[i+2];
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      lip = l+ns;
      li2p = lip+ns;
      qpip = qplus[lip];   /* qplus(i=1,j) */
      qpi2p = qplus[li2p]; /* qplus(i=2,j) */
      if (viscosity != 0.0) {
	/* power-law extrapolation */
	qplus[l] = qpip*exp( log(qpip/qpi2p) * log(r/rip) / log(rip/ri2p) );
      }
      else
	qplus[l] = 0.0;
    }
    /* Now we can update energy with source terms from i=0 */
    for (i = 0; i < nr; i++) {
      for (j = 0; j < ns; j++) {
	l = j+i*ns;
	if (!Cooling) {
	  num = dt*qplus[l] + energy[l];
	  den = 1.0+(ADIABATICINDEX-1.0)*dt*divergence[l];
	  energynew[l] = num/den;
	}
	if (Cooling) {
	  num = EnergyMed[i]*dt*dens[l]/SigmaMed[i] + CoolingTimeMed[i]*energy[l] + dt*CoolingTimeMed[i]*(qplus[l]-QplusMed[i]*dens[l]/SigmaMed[i]);
	  den = dt + CoolingTimeMed[i] + (ADIABATICINDEX-1.0)*dt*CoolingTimeMed[i]*divergence[l];
	  energynew[l] = num/den;
	}
      }
    }
  }
}
   		   
int ConditionCFL (Vrad, Vtheta, deltaT, sys)
     PolarGrid *Vrad, *Vtheta;
     real deltaT;
     PlanetarySystem *sys;
{
  static real Vresidual[MAX1D], Vmoy[MAX1D];
  int i, j, l, ns, nr, lip, ljp;
  real invdt1, invdt2, invdt3, invdt4, cs, newdt, dt, xp, yp, mp;
  int ideb, jdeb;
  real itdbg1, itdbg2, itdbg3, itdbg4, mdtdbg; /* debugging variables */
  real *vt, *vr, dxrad, dxtheta, dvr, dvt, viscr, visct;
  real *soundspeed;
  soundspeed = SoundSpeed->Field;
  ns = Vtheta->Nsec;
  nr = Vtheta->Nrad;
  vt = Vtheta->Field;
  vr = Vrad->Field;
  xp = sys->x[0];
  yp = sys->y[0];
  mp = sys->mass[0];

  newdt = 1e30;
  for (i = 0; i < nr; i++) {
    dxrad = Rsup[i]-Rinf[i];
    dxtheta=Rmed[i]*2.0*PI/(real)ns;
    Vmoy[i] = 0.0;
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      Vmoy[i] += vt[l];
    }
    Vmoy[i] /= (real)ns;
  }
  for (i = One_or_active; i < Max_or_active; i++) {
    dxrad = Rsup[i]-Rinf[i];
    dxtheta=Rmed[i]*2.0*PI/(real)ns;
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      if (FastTransport == YES)
	Vresidual[j] = vt[l]-Vmoy[i];  /* FARGO algorithm */
      else
	Vresidual[j] = vt[l];	       /* Standard algorithm */
    }
    Vresidual[ns]=Vresidual[0];
    for (j = 0; j < ns; j++) {
      l = j+i*ns;
      lip = l+ns;
      ljp = l+1;
      if (j == ns-1) ljp = i*ns;

      cs = soundspeed[l];
/*
      cs = pepsoundspeed(i,j,xp,yp,mp,ns);
*/

      invdt1 = cs/(min2(dxrad,dxtheta));
      invdt2 = fabs(vr[l])/dxrad;
      invdt3 = fabs(Vresidual[j])/dxtheta;
      dvr = vr[lip]-vr[l];
      dvt = vt[ljp]-vt[l];
      if (dvr >= 0.0) dvr = 1e-10;
      else dvr = -dvr;
      if (dvt >= 0.0) dvt = 1e-10;
      else dvt = -dvt;
      invdt4 = max2(dvr/dxrad,dvt/dxtheta);
      invdt4*= 4.0*CVNR*CVNR;
      dt = CFLSECURITY/sqrt(invdt1*invdt1+invdt2*invdt2+invdt3*invdt3+invdt4*invdt4);
      if (dt < newdt) {
	newdt = dt;
	if (debug == YES) {
	  ideb = i;
	  jdeb = j;
	  itdbg1 = 1.0/invdt1; itdbg2=1.0/invdt2; itdbg3=1.0/invdt3; itdbg4=1.0/invdt4;
	  mdtdbg = newdt;
	  viscr = dxrad/dvr/4.0/CVNR/CVNR;     
	  visct = dxtheta/dvt/4.0/CVNR/CVNR;
	}
      }  
    }
  }
  for (i = Zero_or_active; i < MaxMO_or_active; i++) {
    dt = 2.0*PI*CFLSECURITY/(real)NSEC/fabs(Vmoy[i]*InvRmed[i]-Vmoy[i+1]*InvRmed[i+1]);
    if (dt < newdt) newdt = dt;
  }
  if (deltaT < newdt) newdt = deltaT;
  if (debug == YES) {
    printf ("Timestep control information for CPU %d: \n", CPU_Rank);
    printf ("Most restrictive cell at i=%d and j=%d\n", ideb, jdeb);
    printf ("located at radius Rmed         : %g\n", Rmed[ideb]);
    printf ("Sound speed limit              : %g\n", itdbg1);
    printf ("Radial motion limit            : %g\n", itdbg2);
    printf ("Residual circular motion limit : %g\n", itdbg3);
    printf ("Viscosity limit                : %g\n", itdbg4);
    printf ("   Arise from r with limit     : %g\n", viscr);
    printf ("   and from theta with limit   : %g\n", visct);
    printf ("Limit time step for this cell  : %g\n", mdtdbg);
    printf ("Limit time step adopted        : %g\n", newdt);
    if (newdt < mdtdbg) {
      printf ("Discrepancy arise either from shear.\n");
      printf ("or from the imposed DT interval.\n");
    }
  }
  return (int)(ceil(deltaT/newdt));
}

void ComputeSoundSpeed (Rho, Energy)
     PolarGrid *Rho;
     PolarGrid *Energy;
{
  int i, j, l, nr, ns, ii;
  real *dens, *energ, *cs;
  nr = Rho->Nrad;
  ns = Rho->Nsec;
  dens = Rho->Field;
  energ = Energy->Field;
  cs = SoundSpeed->Field;

  for ( i = 0; i < nr; i++ ) {
    for ( j = 0; j < ns; j++ ) {
      l = i*ns + j;
      if (!Adiabatic){
	cs[l] = AspectRatio(Rmed[i])*sqrt(G*1.0/Rmed[i])*pow(Rmed[i], FLARINGINDEX);
      } else
	cs[l] = sqrt( ADIABATICINDEX*(ADIABATICINDEX-1.0)*energ[l]/dens[l] );
    }
  }
}

double pepsoundspeed( i ,j, xp,yp,mp, ns)
int i, j, ns;
double xp, yp, mp;
{
 double  x, y, angle, distanceS, distanceP, omegaS, omegaP, hrs, hrp, cs, plrad, rh, temp, dr;
        angle = (real)j/(real)ns*2.0*PI;
        x = Rmed[i] * cos(angle);
        y = Rmed[i] * sin(angle);
        distanceS = Rmed[i];
        hrs = ASPECTRATIO*distanceS;    
        omegaS=sqrt(1.0/pow(distanceS,3.0));  
 
/*      peplinski*/
       distanceP = sqrt( (x - xp)*(x - xp) + (y - yp)*(y - yp)
                          + pow(THICKNESSSMOOTHING*ASPECTRATIO*sqrt(xp*xp + yp*yp),2.0));

        omegaP=sqrt(mp/pow(distanceP,3.0));
        hrp = HPLANET*distanceP;
        cs = sqrt(omegaS*omegaS + omegaP*omegaP);
        cs *= hrs*hrp;
        cs /= pow((pow(hrs,EOSN) + pow(hrp,EOSN)),1.0/EOSN);

/*      GAUSSIAN ENVELOPE */
/*
        plrad = sqrt(xp*xp + yp*yp);
        rh = plrad*pow(mp/3.0, 1.0/3.0);
        dr = HEATWIDTH*rh;
  
        distanceP = sqrt( (x - xp)*(x - xp) + (y - yp)*(y - yp)); 
        cs = hrs*omegaS;
        temp = pow(distanceP/dr, 2.0);
        cs*= 1.0 + (HEATAMP - 1.0)*exp(-0.5*temp);
*/
        return cs;
}

void ComputePressureField (Rho, Energy, sys)
     PolarGrid *Rho;
     PolarGrid *Energy;
     PlanetarySystem *sys;
{
  int i, j, l, nr, ns, ii;
  real *dens, *pres, *energ, *cs;
  real peq, dp, cscustom, deltacs;
  real xp, yp, mp;
  nr = Rho->Nrad;
  ns = Rho->Nsec;
  dens = Rho->Field;
  pres = Pressure->Field;
  energ = Energy->Field;
  cs = SoundSpeed->Field;
  xp = sys->x[0];
  yp = sys->y[0];
  mp = sys->mass[0];

  for ( i = 0; i < nr; i++ ) {
    for ( j = 0; j < ns; j++ ) {
      l = i*ns + j;
      if (!Adiabatic) {

        cscustom = pepsoundspeed(i,j,xp,yp,mp,ns);
        deltacs = cscustom - cs[l];
        cscustom = cs[l] + deltacs*MassTaper;

	pres[l] = dens[l]*cscustom*cscustom; /* since SoundSpeed is not updated */
                                       /* from initialization, cs remains */ 
                                       /* axisymmetric */
      }
      else
	pres[l] = (ADIABATICINDEX-1.0)*energ[l];
    }
  }
}


void ComputePressureField0 (Rho, Energy)
     PolarGrid *Rho;
     PolarGrid *Energy;
{   
  int i, j, l, nr, ns, ii;
  real *dens, *pres, *energ, *cs;
  real peq, dp, cscustom;
  nr = Rho->Nrad;
  ns = Rho->Nsec;
  dens = Rho->Field;
  pres = Pressure->Field;
  energ = Energy->Field;
  cs = SoundSpeed->Field;

  for ( i = 0; i < nr; i++ ) {
    for ( j = 0; j < ns; j++ ) {
      l = i*ns + j;
      if (!Adiabatic) {
        pres[l] = dens[l]*cs[l]*cs[l]; /* since SoundSpeed is not updated */
                                       /* from initialization, cs remains */
                                       /* axisymmetric */
      }
      else
        pres[l] = (ADIABATICINDEX-1.0)*energ[l];
    }
  }
}




void ComputeTemperatureField (Rho, Energy)
     PolarGrid *Rho, *Energy;
{
  int i, j, l, nr, ns;
  real *dens, *pres, *energ, *temp;
  nr = Rho->Nrad;
  ns = Rho->Nsec;
  dens = Rho->Field;
  pres = Pressure->Field;
  energ = Energy->Field;
  temp = Temperature->Field;
  for ( i = 0; i < nr; i++ ) {
    for ( j = 0; j < ns; j++ ) {
      l = i*ns + j;
      if (!Adiabatic)
	temp[l] = MU/R* pres[l]/dens[l];
      else
	temp[l] = MU/R*(ADIABATICINDEX-1.0)*energ[l]/dens[l];
    }
  }
}

real CircumPlanetaryMass (Rho, sys)
     PolarGrid *Rho;
     PlanetarySystem *sys;
{
  int i, j, l, ns;
  real xpl, ypl;
  real dist, mdcplocal, mdcptotal;
  real *dens, *abs, *ord;
  ns = Rho->Nsec;
  dens = Rho->Field;
  abs = CellAbscissa->Field;
  ord = CellOrdinate->Field;
  xpl = sys->x[0];
  ypl = sys->y[0];
  mdcplocal = 0.0;
  mdcptotal = 0.0;
  if (FakeSequential && (CPU_Rank > 0)) 
    MPI_Recv (&mdcplocal, 1, MPI_DOUBLE, CPU_Rank-1, 0, MPI_COMM_WORLD, &stat);
  for ( i = Zero_or_active; i < Max_or_active; i++ ) {
    for ( j = 0; j < ns; j++ ) {
      l = i*ns + j;
      dist = sqrt ( (abs[l]-xpl)*(abs[l]-xpl) +		\
		    (ord[l]-ypl)*(ord[l]-ypl) );
      if ( dist < HillRadius ) {
	mdcplocal += Surf[i] * dens[l];
      /*if (isnan(mdcplocal))
         printf("%i %f %f %f \n", CPU_Rank, dist, Surf[i], dens[l])*/;
      }
    }
  }
  if (FakeSequential) {
     if (CPU_Rank < CPU_Number-1)
       MPI_Send (&mdcplocal, 1, MPI_DOUBLE, CPU_Rank+1, 0, MPI_COMM_WORLD);
  }
  else
    MPI_Allreduce (&mdcplocal, &mdcptotal, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  if (FakeSequential) {
    MPI_Bcast (&mdcplocal, 1, MPI_DOUBLE, CPU_Number-1, MPI_COMM_WORLD);
    mdcptotal = mdcplocal;
  }
  return mdcptotal;
}
