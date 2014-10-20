program test
  implicit none
  include "mpif.h"
  call configure
end program test

subroutine configure
  implicit none
  include "mpif.h"
  integer ierr, myid_w, nprocs_w
  call MPI_INIT( ierr )
  call MPI_COMM_RANK( MPI_COMM_WORLD, myid_w  , ierr )
  call MPI_COMM_SIZE( MPI_COMM_WORLD, nprocs_w, ierr )
  print*, myid_w, nprocs_w
  call MPI_FINALIZE ( ierr )
  return
end subroutine configure

