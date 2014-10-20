#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <cstdio>
#include <string>
using std::string;

#include "PatchPluto.H"
#include "LoHiSide.H"

// Flag everything as not defined or set
PatchPluto::PatchPluto()
{
  m_isDefined = false;
  m_isBoundarySet = false;
  m_isRiemannSet = false;
  m_isCurrentTimeSet = false;
  m_isCurrentBoxSet = false;
  m_isCurrentDlMinSet = false;
}

PatchPluto::~PatchPluto()
{
}

// Define this object and the boundary condition object
void PatchPluto::define(ProblemDomain& a_domain,
                          const Real&    a_dx,
                          const int&     a_level,
                          const int&     a_numGhost)
{

  // Store the domain and grid spacing
  m_domain = a_domain;
  m_dx = a_dx;
  m_level = a_level;
  m_numGhost = a_numGhost;

  m_isDefined = true;
}

// Factory method - this object is its own factory.  It returns a pointer
// to new PatchPluto object with its boundary condtions defined.
 PatchPluto* PatchPluto::new_patchPluto() const
{
   // Make the new object
   PatchPluto* retval =
     static_cast<PatchPluto*>(new PatchPluto());

   // Set the boundary and Riemann solver values
   retval->setBoundary(left_bound_side,
                       right_bound_side);
   retval->setRiemann(rsolver);

   // Return the new object
   return retval;
}

// Set the current physical time - used for time dependent boundary conditions
void PatchPluto::setCurrentTime(const Real& a_currentTime)
{
  m_currentTime = a_currentTime;
  m_isCurrentTimeSet = true;
}

// Set the current box for the updateSolution call
void PatchPluto::setCurrentBox(const Box& a_currentBox)
{
  m_currentBox = a_currentBox;
  m_isCurrentBoxSet = true;
}

// Set the current minimum cell size - used for GLM_Source
void PatchPluto::setCurrentDlMin(const Real& a_currentDlMin)
{
  m_currentDlMin = a_currentDlMin;
  m_isCurrentDlMinSet = true;
}

// Set the values corresponding to boundary conditions
void PatchPluto::setBoundary(const int *leftbound,
                             const int *rightbound)
{

 for (int i = 0; i < 3; i++)
    {
      left_bound_side[i] =  leftbound[i];
     right_bound_side[i] = rightbound[i];
    }

 m_isBoundarySet = true;

}

// Set the riemann solver to be used
void PatchPluto::setRiemann(Riemann_Solver *solver)
{

   rsolver = solver; 
   m_isRiemannSet = true;

}

// Set the grid[3] structure to be passed to updateState
// (or updateSolution and SWEEP in the next future)
void PatchPluto::setGrid(const Box&    a_box, struct GRID *grid)
{
  int    i, j, k, idim;
  int    np_int, np_tot, np_int_glob, np_tot_glob;
  double scrh, dxmin[3];
  struct GRID *G;

  for (int dir = 0; dir < SpaceDim; ++dir){  
    grid[dir].level       = m_level;
    grid[dir].np_tot      = a_box.size()[dir]+2*m_numGhost;
    grid[dir].np_int      = a_box.size()[dir];
    grid[dir].np_tot_glob = m_domain.domainBox().size()[dir]+2*m_numGhost;
    grid[dir].np_int_glob = m_domain.domainBox().size()[dir];
    grid[dir].nghost      = m_numGhost;
    grid[dir].beg         = a_box.loVect()[dir]+m_numGhost;
    grid[dir].end         = a_box.hiVect()[dir]+m_numGhost;
    grid[dir].gbeg        = m_domain.domainBox().loVect()[dir]+m_numGhost;
    grid[dir].gend        = m_domain.domainBox().hiVect()[dir]+m_numGhost;
    grid[dir].lbeg        = m_numGhost;
    grid[dir].lend        = grid[dir].np_int+m_numGhost-1;
   
    grid[dir].lbound = 0;
    grid[dir].rbound = 0;
   // Set flags to compute boundary conditions
   // is_gbeg (left boundary)
    Box tmp = a_box;
    tmp.shift(dir,-1);
    tmp &= m_domain;
    tmp.shift(dir,1);
    if (tmp != a_box) grid[dir].lbound = 1;
   // is_gend (right_boundary)
    tmp = a_box;
    tmp.shift(dir,1);
    tmp &= m_domain;
    tmp.shift(dir,-1);
    if (tmp != a_box) grid[dir].rbound = 1;

  }

  for (int dir = 0; dir < 3; ++dir){
    G = grid + dir;
    if (dir >= SpaceDim) {
        G->np_int = G->np_tot = G->np_int_glob = G->np_tot_glob = 1;
        G->nghost = G->beg = G->end = G->gbeg = G->gend =
        G->lbeg   = G->lend = G->lbound = G->rbound = 0;
    }
    np_tot_glob = G->np_tot_glob;
    np_int_glob = G->np_int_glob;
    np_tot = G->np_tot;
    np_int = G->np_int;

    G->x  = ARRAY_1D(np_tot, double);
    G->xr = ARRAY_1D(np_tot, double);
    G->xl = ARRAY_1D(np_tot, double);
    G->dx = ARRAY_1D(np_tot, double);

    for (i = 0; i < np_tot; i++){
    #if (GEOMETRY == SPHERICAL) && (LOGR == YES)
     if (dir == IDIR) {
      G->xl[i] = g_domBeg[IDIR]*exp(Real( i+grid[dir].beg-2*grid[dir].nghost )*m_dx);
      G->xr[i] = g_domBeg[IDIR]*exp(Real(i+1+grid[dir].beg-2*grid[dir].nghost)*m_dx);  
      G->x[i] = 0.5*(G->xr[i]+G->xl[i]);
      G->dx[i] = G->xr[i]-G->xl[i];
     } else {
    #endif
      G->xl[i]  = Real(i+grid[dir].beg-2*grid[dir].nghost)*m_dx;
      G->xl[i] += g_domBeg[dir];
      G->x[i]  = G->xl[i]+0.5*m_dx;
      G->xr[i] = G->xl[i]+m_dx;
      G->dx[i] = m_dx; 
    #if (LOGR == YES)
     }
    #endif
    }
  }

  CH_assert(m_isBoundarySet);
  grid[IDIR].lbound *=   left_bound_side[IDIR];
  grid[IDIR].rbound *=  right_bound_side[IDIR];
  grid[JDIR].lbound *=   left_bound_side[JDIR];
  grid[JDIR].rbound *=  right_bound_side[JDIR];
  grid[KDIR].lbound *=   left_bound_side[KDIR];
  grid[KDIR].rbound *=  right_bound_side[KDIR];
  MakeGeometry(grid);

// compute dl_min

 #if (GEOMETRY == CARTESIAN) || ((GEOMETRY == CYLINDRICAL) && (DIMENSIONS ==2))

  grid[IDIR].dl_min = m_dx;
  grid[JDIR].dl_min = m_dx;
  grid[KDIR].dl_min = m_dx;

 #else

  IBEG = grid[IDIR].lbeg; IEND = grid[IDIR].lend;
  JBEG = grid[JDIR].lbeg; JEND = grid[JDIR].lend;
  KBEG = grid[KDIR].lbeg; KEND = grid[KDIR].lend;

  for (idim = 0; idim < 3; idim++)  dxmin[idim] = 1.e30;

  for (i = IBEG; i <= IEND; i++) {
  for (j = JBEG; j <= JEND; j++) {
  for (k = KBEG; k <= KEND; k++) {

    scrh = Length_1(i, j, k, grid);
    dxmin[IDIR] = MIN (dxmin[IDIR], scrh);

    scrh = Length_2(i, j, k, grid);
    dxmin[JDIR] = MIN (dxmin[JDIR], scrh);

    scrh = Length_3(i, j, k, grid);
    dxmin[KDIR] = MIN (dxmin[KDIR], scrh);

  }}}

  grid[IDIR].dl_min = dxmin[IDIR];
  grid[JDIR].dl_min = dxmin[JDIR];
  grid[KDIR].dl_min = dxmin[KDIR];

 #endif

}

// Loop over the patches of a level to assign initial conditions

void PatchPluto::initialize(LevelData<FArrayBox>& a_U)
{

 CH_assert(m_isDefined);

 DataIterator dit = a_U.boxLayout().dataIterator();

 // Iterator of all grids in this level
 for (dit.begin(); dit.ok(); ++dit)
 {
   // Storage for current grid
   FArrayBox& U = a_U[dit()];
   // Set up initial condition in this patch
   startup(U);

 }
   
}

// Number of conserved variables
int PatchPluto::numConserved()
{
   return NVAR;
}

// Generate default names for the conserved variables, "variable#"
Vector<string> PatchPluto::ConsStateNames()
{
  Vector<string> retval;
  int nv;
  char varNameChar[80];

  for (nv = 0; nv < NVAR; nv++){
    if (nv == RHO) retval.push_back("Density");
    if (nv == MX1) retval.push_back("X-momentum");
    if (nv == MX2) retval.push_back("Y-momentum");
    if (nv == MX3) retval.push_back("Z-momentum");
    #if PHYSICS == MHD || PHYSICS == RMHD
     if (nv == BX1) retval.push_back("X-magnfield");
     if (nv == BX2) retval.push_back("Y-magnfield");
     if (nv == BX3) retval.push_back("Z-magnfield");
    #endif
    #if EOS != ISOTHERMAL
     if (nv == ENG) retval.push_back("energy-density");
    #endif

    #if NTRACER > 0   
     if (nv >= TRC && nv < (TRC + NTRACER)){
       sprintf(varNameChar,"rho*tracer%d",nv+1-TRC);
       retval.push_back(string(varNameChar));
     }
    #endif
    #if ENTROPY_SWITCH == YES
     if (nv == ENTR) retval.push_back("rho*entropy");
    #endif

    #ifdef GLM_MHD 
     if (nv == PSI_GLM) retval.push_back("psi_glm");
    #endif

    #if COOLING == SNEq
     if (nv == FNEUT) retval.push_back("rho*fneut");
    #endif

    #if COOLING == MINEq

     if (nv >= HI && nv < NFLX+NIONS){
       sprintf(varNameChar,"rho*fX%d",nv-NFLX);
       retval.push_back(string(varNameChar));
     }
/*
     if (nv == HI)    retval.push_back("rho*fHI");
     if (nv == HeI)   retval.push_back("rho*fHeI");
     if (nv == HeII)  retval.push_back("rho*fHeII");
     if (nv == CI)    retval.push_back("rho*fCI");
     if (nv == CII)   retval.push_back("rho*fCII");
     if (nv == CIII)  retval.push_back("rho*fCIII");
     if (nv == CIV)   retval.push_back("rho*fCIV");
     if (nv == CV)    retval.push_back("rho*fCV");
     if (nv == NI)    retval.push_back("rho*fNI");
     if (nv == NII)   retval.push_back("rho*fNII");
     if (nv == NIII)  retval.push_back("rho*fNIII");
     if (nv == NIV)   retval.push_back("rho*fNIV");
     if (nv == NV)    retval.push_back("rho*fNV");
     if (nv == OI)    retval.push_back("rho*fOI");
     if (nv == OII)   retval.push_back("rho*fOII");
     if (nv == OIII)  retval.push_back("rho*fOIII");
     if (nv == OIV)   retval.push_back("rho*fOIV");
     if (nv == OV)    retval.push_back("rho*fOV");
     if (nv == NeI)   retval.push_back("rho*fNeI");
     if (nv == NeII)  retval.push_back("rho*fNeII");
     if (nv == NeIII) retval.push_back("rho*fNeIII");
     if (nv == NeIV)  retval.push_back("rho*fNeIV");
     if (nv == NeV)   retval.push_back("rho*fNeV");
     if (nv == SI)    retval.push_back("rho*fSI");
     if (nv == SII)   retval.push_back("rho*fSII");
     if (nv == SIII)  retval.push_back("rho*fSIII");
     if (nv == SIV)   retval.push_back("rho*fSIV");
     if (nv == SV)    retval.push_back("rho*fSV");
*/
    #endif
  }
  return retval; 
}

// Generate default names for the primitive variables, "variable#"
Vector<string> PatchPluto::PrimStateNames()
{
  Vector<string> retval;
  int nv;
  static Output output;
  
  if (output.var_name == NULL){
    for (nv = 0; nv < NVAR; nv++){
      output.var_name = ARRAY_2D(64,128,char);
    }
    SetDefaultVarNames(&output);
  }

  for (nv = 0; nv < NVAR; nv++){
    retval.push_back(output.var_name[nv]);
  }
  return retval; 
}

// Number of flux variables
int PatchPluto::numFluxes()
{
  // In some computations there may be more fluxes than conserved variables
  return NVAR;
}

// Number of primitive variables
int PatchPluto::numPrimitives()
{
  // This doesn't equal the number of conserved variables because
  // auxiliary/redundant variable may be computed and stored
  return NVAR;
}

// Component index within the primitive variables of the pressure
int PatchPluto::pressureIndex()
{
  #if EOS != ISOTHERMAL
   return PR;
  #else 
   return -1;
  #endif

}

// Return true if everything is defined and setup
bool PatchPluto::isDefined() const
{
  return m_isDefined        &&
         m_isBoundarySet    &&
         m_isRiemannSet     &&
         m_isCurrentTimeSet &&
         m_isCurrentBoxSet  &&
         m_isCurrentDlMinSet;
}
