module cg
  use precision
  use matdef
  use sparsealg
  use mpi_globals
  implicit none
  private

  public :: conjgrads

contains

  subroutine conjgrads(A,b,x0,x,tol)
    type(rCSR), intent(in) :: A
    real(dp), intent(in) :: b(:)
    real(dp), intent(in) :: x0(:)
    real(dp), intent(inout) :: x(:)
    real(dp), intent(inout) :: tol

    real(dp), allocatable :: r(:), y(:), d(:)
    integer :: k, i, maxiter, nrow, rstart, rend, rovr, itmp
    integer :: nrow_local, nrow_resto
    real(dp) :: alfa, beta, norm1, norm2, sum_local, sum_t
    
    maxiter = size(b)
    nrow = size(b)
    if (id0)  print*,'maxiter=',maxiter
    allocate(r(size(b)))
    allocate(y(size(b)))
    allocate(d(size(b)))

    nrow_local = nrow/numprocs
    nrow_resto = 0
    if (id == numprocs-1) then
       nrow_resto = mod(nrow,numprocs)
    end if
    rstart = id * nrow_local + 1
    rend = (id+1)*nrow_local + nrow_resto

    if (id==0) then
      rovr = 0
      do k = rstart, rend
        itmp = maxval(A%colind(A%rowpnt(k):A%rowpnt(k+1)-1))
        if (itmp > rovr) then
          rovr = itmp
        end if
      end do
      rovr = rovr - nrow_local  
      print*,'overlap: ',rovr, exp(2.0_dp/3.0_dp*log(real(nrow,dp))) 
    end if
    
    call mpi_bcast(rovr, 1, MPI_INT, 0, my_comm, mpierr)
    
    call matvec(A,x0,y,rstart,rend)
    
    r = b - y
    d = r
    x = x0
    sum_local = dot_product(r(rstart:rend),r(rstart:rend))
    call mpi_reduce(sum_local, norm1, 1, MPI_DOUBLE, MPI_SUM, 0, my_comm, mpierr)
    call mpi_bcast(norm1, 1, MPI_DOUBLE, 0, my_comm, mpierr)

    do k = 1, maxiter
       call matvec(A,d,y,rstart,rend)
       sum_local = dot_product(d(rstart:rend),y(rstart:rend))
       call mpi_reduce(sum_local, sum_t, 1, MPI_DOUBLE, MPI_SUM, 0, my_comm, mpierr)     
       call mpi_bcast(sum_t, 1, MPI_DOUBLE, 0, my_comm, mpierr)
       alfa = norm1/sum_t
       x = x + alfa*d
       r = r - alfa*y
       sum_local = dot_product(r(rstart:rend),r(rstart:rend))
       call mpi_reduce(sum_local, norm2, 1, MPI_DOUBLE, MPI_SUM, 0, my_comm, mpierr)
       call mpi_bcast(norm2, 1, MPI_DOUBLE, 0, my_comm, mpierr)
       if (id0) then  
         write(*,'(a,I5,a,E10.3)') 'CG iteration:',k,'  err:',sqrt(norm2)
       end if 
       if (sqrt(norm2) < tol) then
          exit
       end if
       beta = norm2/norm1
       d = r + beta*d
       norm1 = norm2
       call mpi_barrier(my_comm, mpierr)

       call send_vector2(d, rstart, rend, rovr)
    end do

    tol = sqrt(norm2)
    
  end subroutine conjgrads
  

  subroutine send_vector(d,rstart,rend,rovr)
    real(dp) :: d(:)
    integer, intent(in) :: rstart, rend, rovr

    integer :: i
    
    do i = 0, numprocs-2 
      if (id == i) then
        call mpi_send(d(rend-rovr+1), rovr, MPI_DOUBLE, i+1, 0, my_comm, mpierr)
      else if (id == i+1) then
        call mpi_recv(d(rstart-rovr), rovr, MPI_DOUBLE, i,0, my_comm, stats, mpierr)
      end if
   end do
   
   do i = numprocs-1, 1, -1
     if (id == i) then
       call mpi_send(d(rstart), rovr, MPI_DOUBLE, i-1, 1, my_comm, mpierr)
     else if (id == i-1) then
       call mpi_recv(d(rend+1), rovr, MPI_DOUBLE, i, 1, my_comm, stats, mpierr)
     end if
   end do
    
  end subroutine send_vector

  subroutine send_vector2(d,rstart,rend,rovr)
    real(dp) :: d(:)
    integer, intent(in) :: rstart, rend, rovr

    integer :: i
    integer :: req(4)
    integer :: statsarr(MPI_STATUS_SIZE,4)

    req = MPI_REQUEST_NULL

    if (id > 0) then
      call mpi_isend(d(rstart), rovr, MPI_DOUBLE, id-1, 0, my_comm, req(1), mpierr)
    end if
    if (id < numprocs-1) then
      call mpi_isend(d(rend-rovr+1), rovr, MPI_DOUBLE, id+1, 1, my_comm, req(2), mpierr)
    end if
   
    if (id < numprocs-1) then
      call mpi_irecv(d(rend+1), rovr, MPI_DOUBLE, id+1, 0, my_comm, req(3), mpierr)
    end if
    if (id > 0) then
      call mpi_irecv(d(rstart-rovr), rovr, MPI_DOUBLE, id-1, 1, my_comm, req(4), mpierr)
    end if
    
    call mpi_waitall(4, req, statsarr, mpierr)

  end subroutine send_vector2
  
end module cg
