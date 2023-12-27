module matdef
  use precision
  implicit none
  private

  public :: rCSR, rCOO
  public :: create
  
  type rCSR
     integer :: nrow
     integer :: nnz    
     real(dp), allocatable :: nzval(:)
     integer, allocatable :: colind(:)
     integer, allocatable :: rowpnt(:)    
  end type rCSR
     
  type rCOO
     integer :: nrow
     integer :: nnz
     real(dp), allocatable :: nzval(:)
     integer, allocatable :: colind(:)
     integer, allocatable :: rowind(:)     
  end type rCOO

  interface create
     module procedure create_rCSR
     module procedure create_rCOO
  end interface create
  

contains

  subroutine create_rCSR(M, nrow, nnz)
    type(rCSR), intent(inout) :: M
    integer, intent(in) :: nrow
    integer, intent(in) :: nnz

    M%nrow = nrow
    M%nnz = nnz
    
    allocate(M%nzval(nnz))
    allocate(M%colind(nnz))
    allocate(M%rowpnt(nrow+1))
    
  end subroutine create_rCSR

  subroutine create_rCOO(M, nrow, nnz)
    type(rCOO), intent(inout) :: M
    integer, intent(in) :: nrow
    integer, intent(in) :: nnz

    M%nrow = nrow
    M%nnz = nnz
    
    allocate(M%nzval(nnz))
    allocate(M%colind(nnz))
    allocate(M%rowind(nnz))
    
  end subroutine create_rCOO
  
  
end module matdef
