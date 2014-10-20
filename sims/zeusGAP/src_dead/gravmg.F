c=======================================================================
c
c                            Developed by
c                Laboratory for Computational Astrophysics
c                  University of California at San Diego
c
      subroutine grav3D_MG
c
c     3-D gravitational potential solver using James Bordner's
c     MGMPI package for Poisson's equation
c
      use real_prec
      use config
      use param
      use cons
      use domain
      use root
      use field
      use grid
      use bndry
      use mpiyes
      use mpipar
c
      implicit NONE
      return
      end
