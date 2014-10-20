c=======================================================================
c
c    \\\\\\\\\\        B E G I N   P R O G R A M          //////////
c    //////////               Z E U S M P                 \\\\\\\\\\
c
c                            Developed by
c                Laboratory of Computational Astrophysics
c               University of Illinois at Urbana-Champaign
c
c=======================================================================
c
      program zeusmp
c
c PURPOSE
c   Main program for 3-D MPI version of ZEUS.
c
c AUTHOR
c   Robert A. Fiedler
c
c LAST MODIFIED by JCH, for F90
c   6/26/02.
c.......................................................................
c
c DECLARATIONS
c
      use real_prec
      use config
      use param
      use field
      use grid
      use root
      use scratch
#ifdef MPI_USED
      use mpiyes
#else
      use mpino
#endif
      use mpipar
      use gravmod
      use clockmod
      use impsoln
      use planet
      use gboundary
c
      implicit NONE
c

      real(rl4) :: cputime, wclock
c
      real(rl) :: zcs, etot, etot_glb
      real(rl) :: cpuall
c
      integer :: i , j , k, nwrite
      integer :: nx, ny, nz, snz, maxn
c

      ifsen(1) = 0
      ifsen(2) = 0
      ifsen(3) = 1
      ifsen(4) = 1
      ifsen(5) = 1
c
      myid_w    = 0
      myid      = 0
      nprocs_w  = 1
      nprocs    = 1
      coords(1) = 0
      coords(2) = 0
      coords(3) = 0
      reorder   = .false.
      totcgit = 0
      ncgcall = 0
c
c  call CONFIGURE
c
      call configure
c
c  Master writes greeting.
c
      if (myid_w .eq. 0) then
       call options
      endif
c
c Set up the problem: read input deck and possibly a restart file.
c
      call mstart
c
c Write out initial data dumps.
c
      call dataio( ifsen(2), ifsen(3), ifsen(4), ifsen(5), ifsen(6))
c
c Create best calculating plans for FFTW to perform forward
c and backward FFT.
c 
#ifdef FFT
      if(xgrvfft) then
        nx=ie-is+1
        ny=je-js+1
        nz=ke-ks+1
        call create_plan(ntiles(1)*nx,ntiles(2)*ny,ntiles(3)*nz)
      endif
#endif
c
c  Initialize cpu and wall clocks.  The quantities "cputime" and
c  "wclock" are the CPU and wall clock times (in seconds) in the main 
c  loop.
c
        wclock0 = 0.0
        cputime0 = 0.0
        call clocks (cputime, wclock)
        wclock0 = wclock
        cputime0 = cputime
c
      if (myid .eq. 0)
     &  write(6,"(/,' Set-up complete with ',i2,' warning(s):'
     &             ,' entering main loop...')") nwarn
c
      if(ldimen .eq. 1) nwrite = 1000
      if(ldimen .eq. 2) nwrite = 100
      if(ldimen .eq. 3) nwrite = 10
c
c--------------------------  start of main loop  -----------------------
c
c Execution ends when INTCHK returns a value of 1 for ifsen(1).
c
1000  continue

      nsub = 1

c
c Solve Poisson's equation for the gravitational potential.
c
      if(xgrav .or. xgrvfft) then
         if(time.ge.tsg) call gravity ! We can turn on SG at a later time if we wish. 
      endif                     ! xgrav
      
c     
c     Indirect potential. Get the integrals over the disk domain. 
c
      call indirect_disc_pot

c Get planet mass at time t.
      
      if(time .lt. planet_on*2d0*pi) then
         masstaper = 0d0
      else if((time .ge. planet_on*2d0*pi).and.
     &        (time .lt. (planet_on + switch_on)*2d0*pi)) then
         masstaper = (time - planet_on*2d0*pi)/(switch_on*2d0*pi)
         masstaper = (sin(masstaper*pi/2d0))**2
      else 
         masstaper = 1d0
      endif
      
      mplanet = planetmass*masstaper

c Evolve the planet motion from t to t+dt, where dt is one hydro 
c time step, using hydro info at t. Planet position at t is stored as
c xp, yp, zp for use in source step of hydro. 

      call update_planet_xyz

c
c     Preperation for Godon's non-reflecting BC. ***TO DELETE***
c

      if(godonbc.eq..true.) call godon_nrbc

      call MPI_BARRIER(comm3d, ierr )


c
c Evaluate all non-advective terms in the evolution equations.
c     
      if(lrad .eq. 0) then
         if(myid_w .eq. 0) then
            if(mod(nhy,nwrite) .eq. 0) then
               write(*,"('nhy, time, dt = ',i8,1p2d12.4)")
     .              nhy, time, dt
            endif               ! mod
         endif                  ! myid_w
      endif                     ! lrad
c     
      call srcstep
c
c Compute the advection of all field variables.
c
      if(ldimen .eq. 3) call transprt
      if(ldimen .eq. 2) call transprt_2D
      if(ldimen .eq. 1) call transprt_1D
c
c Update the step counter and current time value.
c
      nhy   = nhy   + 1
      time  = time  + dt

c
c Check the CPU time, number of steps, output times to determine if
c a stopping criterion has been met or output is desired.
c Also check for keyboard input, depending on the value of mbatch.
c
      call intchk( ifsen(2), ifsen(3), ifsen(4), ifsen(5), ifsen(6) )
c
c Compute new timestep
c
      call nudt
c
c Update the grid and related quantites.
c
      if(xvgrid) call newgrid
c
c Write out any desired output files now that everything has been
c updated.
c Skip dataio if the run is being terminated to avoid duplicate output.
c
      if (ifsen(1) .eq. 1) goto 2000
      if(myid.eq.0) call output_planetxyz
      call dataio( ifsen(2), ifsen(3), ifsen(4), ifsen(5), ifsen(6))
c
      goto 1000  !  Loop back to begin the next time step.
c
c--------------------------  end of main loop  -------------------------
c
c Terminate the run by making final dumps, write goodbyes
c
2000  continue
      call clocks (cputime, wclock)
#ifndef ARCH_CRAY
      tused = real(cputime)
#else
      tused = wclock
#endif
      ifsen(2) = 1
      ifsen(3) = 1
      ifsen(4) = 1
      ifsen(5) = 1
      ifsen(6) = 1
      if(myid.eq.0) call output_planetxyz  
      call dataio( ifsen(2), ifsen(3), ifsen(4), ifsen(5) , ifsen(6))
c
#ifdef MPI_USED
c
c Sum everyone's cputime (stored in tused) to get CPU time for all 
c processes.
c
      call MPI_REDUCE(tused, cpuall, 1, MPI_FLOAT,
     &                MPI_SUM, 0, comm3d, ierr )
#else /* MPI */
       cpuall = tused
#endif /* MPI */
      if (myid .eq. 0) then      
c
c Let's assume tused is user + system time on the master thread.
c One would also like to keep track of wall-clock time and the sum
c of CPU times used by each processor.
c
        zcs = real(nprocs_w*nhy*nx1z*nx2z*nx3z)/(tused+tiny)
        write(6,"(/' Execution terminated with ',i4,' warning(s)')") 
     &     nwarn
        write(6,"(/' Performance summary:')")
        write(6,"('  zone-cycles per cpu second =',1pe12.5)") zcs
        write(6,"('  Master CPU                 =',1pe12.5, ' sec')") 
     &     tused
        write(6,"('  Average CPU/node           =',1pe12.5, ' sec')") 
     &     cpuall/dble(nprocs_w)
        write(6,"('  Wall Clock                 =',1pe12.5, ' sec')") 
     &     wclock
        write(6,"()")
c
        if(xgrav .and. (.not. xsphgrv)) then
         if(lgeom .ne. 1) then
          write(6,"(/' GRAVITY SUMMARY:')")
          write(6,"('  Percentage of cycles with Phi updates: ',i3,
     .              '%')")int(100*float(ncgcall)/float(nhy))
          write(6,"('  Average number of iterations/update  : ',
     .              1pd12.4)")float(totcgit)/float(ncgcall)
         endif ! lgeom
        endif ! xgrav
        if(lrad .ne. 0) then
         write(6,"(/' RADIATION SUMMARY:')")
         write(6,"('  Average number of N-R iterations/cycle: ',
     .              1pd12.4)")totnrit/float(nhy)
         write(6,"('  Average number of CG  iterations/cycle: ',
     .              1pd12.4)")totlsit/float(nhy)
        endif ! lrad
c
        close(unit=2)
C        close(unit=3)
        close(unit=30)
       if(xtsl) then
        close(unit=31)
       endif ! xtsl
      endif
c
c Turn off MPI
c
#ifdef MPI_USED
      call MPI_FINALIZE ( ierr )
#endif
c
c=======================================================================
c
c    \\\\\\\\\\          E N D  P R O G R A M             //////////
c    //////////               Z E U S M P                 \\\\\\\\\\
c
c=======================================================================
c
      end



c------------------------------------------------------------------------c
c---------------------PLANET STUFF---------------------------------------c
c------------------------------------------------------------------------c



      subroutine update_planet_xyz
      use root
      use planet
      use mpiyes
      use mpipar
      
c
c     Compute the the planet position and velocities at t + dt, where t is the current time
c     and dt is the hydro time-step. 

      implicit none
      integer:: nok, nbad
      real*8 :: hillradius
      real*8 :: newtime
      external :: derivs, rkqs

      newtime = time + dt

      xp = planet_info(1) 
      yp = planet_info(2)
      zp = planet_info(3)

      vxp = planet_info(4)
      vyp = planet_info(5)
      vzp = planet_info(6)

      hillradius = sqrt(xp**2+yp**2+zp**2)*(planetmass/3d0)**(1d0/3d0)
      planetsoft = softeps*hillradius

c     
c     Is it time to output information by root CPU? 
c
c     dpforce(7:9) z-torque inside rp, outside rp and total,
c     dpforce(10:12) z-torque inside rp, outside rp and total, WITH TAPERING
c     dpforce(13) mass inside Hill sphere
c     
c     Output planet position and velocities to planetxy_dt.dat 
c
c     These are torques and forces at t (NOT t + dt)
     
      if ((time  .ge. (ttq+dttq)).or.(ifsen(1).eq.1)) then
         ttq = ttq + dttq
         call disc_planet_torques(xp, yp, zp)
         if(myid.eq.0) then 
            open(unit=926, file="planetxy_dt.dat", status="old"
     &           ,position="append")        
            open(unit=928, file="dptorque.dat", status="old"
     &           ,position="append")
            write(928,fmt="(1x,8(E22.15,1x))") time, dpforce(7),
     &           dpforce(8), dpforce(9), dpforce(10), dpforce(11), 
     &           dpforce(12), dpforce(13)
            write(926,fmt="(1x,7(E22.15,1x))") time, xp,
     &           yp, zp, vxp,
     &           vyp, vzp
            close(926)
            close(928)
         endif
      endif
      
      if(time.lt.release) then
c     
c     No migration yet so these are the analytic position and 
c     velocities at t + dt
         
         planet_info(1) =  planetrad*cos(planetomega*(newtime))
         planet_info(2) =  planetrad*sin(planetomega*(newtime))
         planet_info(3) =  0d0
         
         planet_info(4) = -planetrad*planetomega*sin(planetomega*
     &        (newtime))
         planet_info(5) =  planetrad*planetomega*cos(planetomega*
     &        (newtime))
         planet_info(6) =  0d0
      endif
      
      if(time.ge.release) then  
c
c     Planet migration allowed. Need to call Runge-Kutta integrator to update planet 
c     position and velocities. 
c
c     Should probably upgrade to RK78 integrator in the future. 

         call odeint(planet_info,6,time,newtime,1d-6,dt/5d0,0d0,nok,
     &        nbad, derivs, rkqs)
         
      endif
c
c     planet_info now updated with planet position and velocity at t+dt 
c
      return
      end
 


      subroutine forces_on_planet(plx, ply, plz)
      use planet
      implicit none
c
c     Obtaint the total forces on the planet.
c
      real*8 :: plx, ply, plz

c     Disk-on-planet forceds
      call disc_planet_torques(plx, ply, plz)

c     Other forces: star and indirect potential
      call other_forces_on_planet(plx, ply, plz) 

      totforce(1) = dpforce(1) + otherforces(1) ! Total force on planet, disc planet force WITHOUT TAPERING in Hill sphere
      totforce(2) = dpforce(2) + otherforces(2)
      totforce(3) = dpforce(3) + otherforces(3)

      totforce(4) = dpforce(4) + otherforces(1) ! Total force on planet, disc planet force WITH TAPERING in Hill sphere
      totforce(5) = dpforce(5) + otherforces(2)
      totforce(6) = dpforce(6) + otherforces(3)

      return
      end


      subroutine output_planetxyz 
      use root
      use planet
      
      implicit none
      
c     
c     Output planet information (only root CPU calls this routine).
c     The output format matches that of hdf and restart files provided by ZEUS.
c     These are to be used to restart a run (planetxy_restart.dat) and data analysis 

      if ((time  .ge. (thdf+dthdf)).or.(ifsen(1).eq.1)) then
         open(unit=924, file="planetxy_hdf.dat", status="old"
     &        ,position="append") 
         write(924,fmt="(1x,7(E22.15,1x))") time, planet_info(1),
     &        planet_info(2), planet_info(3), planet_info(4),
     &        planet_info(5), planet_info(6)
         close(924)
      endif
      
      if ((time .ge. (tdump+dtdump)).or.(ifsen(1).eq.1)) then
         open(unit=925, file="planetxy_restart.dat", status="old"
     &        ,position="append")
         write(925,fmt="(1x,7(E22.15,1x))") time, planet_info(1),
     &        planet_info(2), planet_info(3), planet_info(4),
     &        planet_info(5), planet_info(6)
         close(925)
      endif      
      return
      end 
      
