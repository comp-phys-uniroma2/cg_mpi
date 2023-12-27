program poisson
  use precision
  use sparsealg
  use matdef
  use adj_map
  use cg
  use mpi_globals
  use clock
  implicit none

  integer :: N 
  type(adjlist), allocatable :: graph(:)
  real(dp), allocatable :: phi(:), rhs(:)
  real(dp) :: dx, tol, L
  real(dp), parameter :: Pi = 4.0_dp*atan(1.0_dp)
  integer :: node, nnodes, nnz, i, j, k, nn, fd
  integer :: nrow_local, nrow_resto, rstart, rend
  type(rCSR) :: A_csr
  integer :: omp_get_max_threads

  call mpi_initialize()


  call mpi_reduce(nrow_local, N, 1, mpi_int, mpi_sum, 0, my_comm, mpierr)
  
  call mpi_bcast(N, 1, mpi_int, 0, my_comm, mpierr)

  if (id0) then
    print*,'num procs:',numprocs
    print*,'num threads:',omp_get_max_threads()
  end if

  !|---+---+---+---|
  !1   2   3   4   5
  
  L = 1.0_dp
  N = 40
  tol = 0.001
  dx = L/real(N-1,dp)
  
  nnodes = N*N*N
  
  call create_3D(graph, N, N, N)

  nrow_local = nnodes/numprocs
  nrow_resto = 0
  if (id == numprocs-1) then
    nrow_resto = mod(nnodes,numprocs)
  end if
  rstart = id * nrow_local + 1
  rend = (id+1)*nrow_local + nrow_resto
  
  write(*,*) 'node ',id,' rows:',rstart,rend

  nnz = count_nnz(graph,rstart,rend) 
  
  call create(A_csr, nnodes, nnz)

  k = 0
  A_csr%rowpnt(1) = 1
  do i = 1, rstart-1
     A_csr%rowpnt(i+1)=1 
  end do
  do i = rstart, rend
    k = k + 1
    A_csr%nzval(k) = -6.0_dp
    A_csr%colind(k) = i
    nn = size(graph(i)%neignode)
    do j = 1, nn
       k = k + 1
       A_csr%nzval(k) = 1.0_dp
       A_csr%colind(k) = graph(i)%neignode(j)
    end do
    A_csr%rowpnt(i+1) = k + 1
  end do
  do i = rend+1, nnodes
     A_csr%rowpnt(i+1)=k + 1 
  end do

 allocate(rhs(nnodes))
 allocate(phi(nnodes))

 rhs = 0.0_dp
 phi = 0.0_dp

 rhs(coo2node(N/2, N/2, N/2+10)) = 4.0_dp * Pi * 14.4 / dx
 rhs(coo2node(N/2, N/2, N/2-10)) = -4.0_dp * Pi * 14.4 / dx 

 if (id0) call set_clock()
 
 call conjgrads(A_csr, rhs, phi, phi, tol) 

 if (id0) call write_clock()
 
 call mpi_close()
 ! open(100, file = 'sol.dat'...)

 if (id0) then
   print*,'output solution'    
   open(newunit=fd, file='sol.dat')
   do i = 1, N
      do j = 1, N
        write(fd,*) phi(coo2node(N/2,i,j)) 
      end do
      write(fd,*)
   end do
   close(fd)
   print*,'done'    
 end if

 
contains
  

  function coo2node(i,j,k) result(node)
    integer, intent(in) :: i,j,k
    integer :: node

   node = (k-1)*N*N + (j-1)*N + i

  end function coo2node
 
 
 
end program poisson
 
 
