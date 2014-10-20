c=======================================================================
c
c                            Developed by
c                Laboratory of Computational Astrophysics
c               University of Illinois at Urbana-Champaign
c
      subroutine msave(filename)
c
c  PURPOSE:  Writes [reads] all common block variables to the logical
c  unit 4.  Currently, the common blocks written [read] are:
c            common /gridvarr/ = real    grid     variables
c            common /gridvari/ = integer grid     variables
c            common /fieldr  / = real    field    variables
c            common /fieldi  / = integer field    variables
c            common /bndryr  / = real    boundary variables
c            common /bndryi  / = integer boundary variables
c            common /rootr   / = real    root     variables
c            common /rooti   / = integer root     variables
c            common /gravcomr/ = real    gravity  variables
c            common /gravcomi/ = integer gravity  variables
c
c  The following blocks are NOT written [read]:
c            common /frmlcomr/ = real    formal soln variables
c            common /frmlcomi/ = integer formal soln variables
c            common /mmntcomr/ = real    moment soln variables
c            common /mmntcomi/ = integer moment soln variables
c 
c  EXTERNALS: [none]
c
c  LOCALS:
c
c  MODIFIED: 26 Aug. 1996 by RAF for ZEUS-MP.
c  MODIFIED:  6 Jan. 1999 by JCH
c  MODIFIED:  7 Jan. 1999 by JCH (increased nrootr for root.h;
c                                 makes room for rad t-step var.)
c  MODIFIED: 18 May. 1999 by efh (increased nrootr due to tslice)
c  MODIFIED: December, 2005 by JCH; increased nfieldr for abun array
c-----------------------------------------------------------------------
      use real_prec
      use config
      use param
      use grid
      use mpiyes
      use mpipar
      use field
      use bndry
      use root
      use gravmod
      use restart_arrays
c
c
      implicit NONE
c
      character*15 :: filename
c
c=======================================================================
c
      open(unit=4,file=filename,status='unknown',form='unformatted')
c=======================================================================
c
      mgridr = 40*in + 36*jn + 24*kn
      ngridr = 27*in + 23*jn + 15*kn + 3
      ngridi = 13
      allocate(rlgrdvr(ngridr+mgridr))
      allocate(ntgrdvr(ngridi))
c
      nfieldr                  =           5*in*jn*kn
      if(xmhd        ) nfieldr = nfieldr + 3*in*jn*kn
      if(lrad  .ne. 0) nfieldr = nfieldr +   in*jn*kn
      if(xgrav       ) nfieldr = nfieldr +   in*jn*kn
      if(nspec .gt. 1) nfieldr = nfieldr +   in*jn*kn*nspec
      allocate(rlfldvr(nfieldr))
c
      nbdryr = 60*jn*kn + 60*in*jn + 60*in*kn + 6*nbvar
      nbdryi = 10*jn*kn + 10*in*jn + 10*in*kn + 6*nbvar
     .                  + 18
      allocate(rlbdryvr(nbdryr))
      allocate(ntbdryvr(nbdryi))
c
      call mapout
c
      write(4)  rl grd  vr , nt grd  vr
     &         ,rl fld  vr 
     &         ,rl bdry vr , nt bdry vr
     &         ,rl rt   vr , nt rt   vr , ch rt   vr
     &         ,rl grv  vr , nt grv  vr
      if (myid .eq. 0)
     &  write(2,"(/1x,'restart dump written at time=',1pe12.5,' cycle='
     & ,i6)") time,nhy
      close(unit=4)
c
      deallocate(rlgrdvr)
      deallocate(ntgrdvr)
      deallocate(rlfldvr)
      deallocate(rlbdryvr)
      deallocate(ntbdryvr)
c
      return
c
c-----------------------------  MGET  ----------------------------------
c
      entry mget(filename)
      open(unit=4,file=filename,status='old',form='unformatted')
c
      mgridr = 40*in + 36*jn + 24*kn
      ngridr = 27*in + 23*jn + 15*kn + 3
      ngridi = 13
      allocate(rlgrdvr(ngridr+mgridr))
      allocate(ntgrdvr(ngridi))
c
      nfieldr                  =           5*in*jn*kn
      if(xmhd        ) nfieldr = nfieldr + 3*in*jn*kn
      if(lrad  .ne. 0) nfieldr = nfieldr +   in*jn*kn
      if(xgrav       ) nfieldr = nfieldr +   in*jn*kn
      if(nspec .gt. 1) nfieldr = nfieldr +   in*jn*kn*nspec
      allocate(rlfldvr(nfieldr))
c
      nbdryr = 60*jn*kn + 60*in*jn + 60*in*kn + 6*nbvar
      nbdryi = 10*jn*kn + 10*in*jn + 10*in*kn + 6*nbvar
     .                  + 18
      allocate(rlbdryvr(nbdryr))
      allocate(ntbdryvr(nbdryi))
c
      read(4)   rl grd  vr , nt grd  vr
     &         ,rl fld  vr
     &         ,rl bdry vr , nt bdry vr
     &         ,rl rt   vr , nt rt   vr , ch rt   vr
     &         ,rl grv  vr , nt grv  vr
      close(unit=4)
c
      call mapin
c
      deallocate(rlgrdvr)
      deallocate(ntgrdvr)
      deallocate(rlfldvr)
      deallocate(rlbdryvr)
      deallocate(ntbdryvr)
c
      return
      end
