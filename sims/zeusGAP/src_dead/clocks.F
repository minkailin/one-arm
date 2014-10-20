c=======================================================================
c
       subroutine clocks (cputime, wclock)
c
c This routine obtains the CPU and wall-clock times in seconds for the 
c calling process since the times (cputime0, wclock0) at which the 
c clocks were initalized (passed through common block /clocks/). 
c
c With MPI, MPI_WALL time is used except on the EXEMPLAR.
c
c For systems without a CPU timer, cputime is set to zero.
c
c Written by RAF, last modified 3/25/96.
c......................................................................
c
      use real_prec
      use param
      use config
      use clockmod
c
      implicit NONE
c
      real(rl4) :: cputime, wclock
c
c
      real(rl) :: etime
c
c
      real(rl) :: MPI_WTIME, wall
c
c
c Wall clock time is easy to get with MPI:
c
      wall   = MPI_WTIME()
      wclock = real(wall) - wclock0
c
c Get the CPU time for this process/thread.
c
c For systems with the standard UNIX etime/itime routines.  Note that
c the output from itime is an integer array with (hours,minutes,seconds)
c and therefor is accurate only to the nearest second.  The output from
c function etime itself is the sum of user plus system times, which can
c be significantly longer than the user time alone (stored in 
c tarray(1)). 
c
c
      cputime = etime ( tarray ) - cputime0
c
      end
