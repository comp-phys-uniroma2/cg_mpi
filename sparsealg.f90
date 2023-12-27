module sparsealg
  use precision
  use matdef
  implicit none
  private

  public :: matvec

contains

  subroutine matvec(A,x,y,rstart,rend)
    type(rCSR), intent(in) :: A
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: y(:)
    integer, intent(in) :: rstart, rend

    integer :: k, i, is, ie
    real(dp) :: ss

    !$OMP PARALLEL PRIVATE(i,k,ss,is,ie)
    !$OMP DO SCHEDULE(static)
    do k = rstart, rend
       ss = 0.0_dp 
       is=A%rowpnt(k)
       ie=A%rowpnt(k+1)-1
       do i = is, ie
         ss = ss + A%nzval(i)*x(A%colind(i))
       end do
       y(k) = ss
    end do
    !$OMP END DO
    !$OMP END PARALLEL
    
  end subroutine matvec

end module sparsealg
       

    
    
    
    
  
