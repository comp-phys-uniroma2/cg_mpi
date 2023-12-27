module mpi_globals
  use mpi 

  integer, public ::  my_comm, id, numprocs, mpierr, numchar
  logical, public ::  id0
  character(MPI_MAX_PROCESSOR_NAME), public :: hostname
  integer, public :: stats(MPI_STATUS_SIZE)

  public :: mpi_initialize

contains

 subroutine mpi_initialize
   my_comm = MPI_COMM_WORLD
   call MPI_INIT( mpierr )                                                 !#MPI#
   call MPI_COMM_SIZE( my_comm, numprocs, mpierr )                        !#MPI#
   call MPI_COMM_RANK( my_comm, id, mpierr )                              !#MPI#
   call MPI_GET_PROCESSOR_NAME(hostname,numchar, mpierr)
   if(id==0) id0=.true.
 end subroutine mpi_initialize

 subroutine mpi_close
   call MPI_FINALIZE(mpierr)
 end subroutine mpi_close

end module mpi_globals
