#ifndef PATCH_GRID_H
#define PATCH_GRID_H

extern "C" {
#include "pluto.h"
}

#include "MayDay.H"
#include "PiecewiseLinearFillPatch.H"
#include "FineInterp.H"

class PatchGrid
{
public:

 PatchGrid(){;}

  ~PatchGrid(){

    if (PGrid[0].x != NULL) { 
      FreeGrid(PGrid);  
    } 
  }

 struct GRID* getGrid()
 {
   return PGrid;
 }

protected:

 struct GRID PGrid[3];

private:


};

class PiecewiseLinearFillPluto : public PiecewiseLinearFillPatch
{
public:

  PiecewiseLinearFillPluto(){m_is_defined = false;}

 ~PiecewiseLinearFillPluto(){;}

  int get_geometry ()
  {
   #if GEOMETRY == CARTESIAN
     return 1;
   #endif

   #if GEOMETRY == CYLINDRICAL
    #if DIMENSIONS == 2
     return 2;
    #else
     MayDay::Error("Cylindrical geometry supported in Pluto/Chombo only in 2 dimensions");
    #endif
   #endif

   #if GEOMETRY == SPHERICAL
//    MayDay::Error("Spherical geometry not supported in Pluto/Chombo");
    #if DIMENSIONS >= 2
     #if LOGR == NO
      return 3;
     #else
      return 4;
     #endif
    #else
     MayDay::Error("Spherical geometry supported in Pluto/Chombo only in 2 dimensions");
    #endif
   #endif

   #if GEOMETRY == POLAR
     MayDay::Error("Only cartesian and 2D cylindrical/spherical geometries are supported in Pluto/Chombo");
   #endif
  }

};

class FineInterpPluto : public FineInterp
{
public:

  FineInterpPluto(){is_defined = false;}

 ~FineInterpPluto(){;}

  int get_geometry ()
  {
   #if GEOMETRY == CARTESIAN
     return 1;
   #endif

   #if GEOMETRY == CYLINDRICAL
    #if DIMENSIONS == 2
     return 2;
    #else
     MayDay::Error("Cylindrical geometry supported in Pluto/Chombo only in 2 dimensions");
    #endif
   #endif

   #if GEOMETRY == SPHERICAL
    #if DIMENSIONS >= 2
     #if LOGR == NO
      return 3;
     #else
      return 4;
     #endif
    #else
     MayDay::Error("Spherical geometry supported in Pluto/Chombo only in 2 dimensions");
    #endif
   #endif

   #if GEOMETRY == POLAR
     MayDay::Error("Only cartesian and 2D cylindrical/spherical geometries are supported in Pluto/Chombo");
   #endif
  }

};

#endif
