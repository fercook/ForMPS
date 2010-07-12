module MatrixProductState_Class

  use ErrorHandling
  use Constants
  use Matrix_Helper
  use MPSTensor_Class

  implicit none

  integer,parameter :: MAX_LENGTH = 200

!###############################
!#####  The class main object
!###############################
  type MatrixProductState
     private
     integer :: Length
     type(MPSTensor),Allocatable :: Tensors(:)
     logical :: initialized=.false.
     logical :: OpenBoundaryConditions=.true.
     integer :: CanonizedAtSite=-2
     integer :: BondDimension
   contains
     procedure :: delete => delete_MatrixProductState
     procedure :: isCanonized => is_full_MPS_Canonized
     procedure :: CanonizeAtSite => Canonize_MPS_At_Site
     procedure :: LCanonize => Canonize_MPS_At_Left
     procedure :: RCanonize => Canonize_MPS_At_Right
!     procedure :: Canonize => Canonize_MPS_At_Right
     procedure :: L2Norm => Norm_of_MatrixProductState
!!$TODO
!!$    ExpandBondTo
!!$    SaveToFile
!!$    ReadFromFile
!!$    ReplaceTensorBy
!!$    GetTensor
  end type MatrixProductState

  interface new_MatrixProductState
     module procedure new_MatrixProductState_Random_OBC
  end interface

  interface assignment (=)
     module procedure  new_MatrixProductState_from_Assignment
  end interface

contains

!######################################    constructor
  function new_MatrixProductState_Random_OBC(Length,spin,bond) result(this)
    type(MatrixProductState) :: this
    integer :: spin,bond
    integer :: loopIndex
    integer :: Length

    if(Length.gt.MAX_LENGTH) then
       call ThrowException('new_MatrixProductState_Random_OBC','Length larger than maximum',NoErrorCode,CriticalError)
    endif

    this%Length=Length
    this%BondDimension=bond
    allocate (this%Tensors(Length))

    this%Tensors(1) = new_MPSTensor(spin,1,bond)
    this%Tensors(Length) = new_MPSTensor(spin,bond,1)

    do loopIndex=2,Length-1
       this%Tensors(loopIndex) = new_MPSTensor(spin,bond,bond)
    enddo
    this%CanonizedAtSite=-2
    this%initialized=.true.

  end function new_MatrixProductState_Random_OBC




!######################################    delete
  integer function delete_MatrixProductState(this)
    class(MatrixProductState) :: this
    integer n, status

    do n=1,this%Length
       status=this%Tensors(n)%delete()
       if(WasThereError().or.status.ne.Normal) then
          call ThrowException('delete_MatrixProductState','Error in memory trying to delete tensors',NoErrorCode,status)
          delete_MatrixProductState=CriticalError
          return
       endif
    enddo

    deallocate(this%Tensors)
    this%Length=0
    this%initialized=.false.

    delete_MatrixProductState=Normal

  end function delete_MatrixProductState
    

  subroutine new_MatrixProductState_from_Assignment(lhs,rhs) 
    class(MatrixProductState),intent(out) :: lhs
    type(MatrixProductState),intent(in) :: rhs
    integer :: n

     if(.not.rhs%initialized) then
        call ThrowException('new_MatrixProductState_from_Assignment','Original MPS not initialized',NoErrorCode,CriticalError)
        return
     endif
    
    lhs%Length=rhs%Length
    lhs%OpenBoundaryConditions=rhs%OpenBoundaryConditions

!    The following check should be
!          if(lhs%initialized_) deallocate(lhs%Tensors)
!    But instead I check for allocated because of a bug in GFortran,
!    http://gcc.gnu.org/bugzilla/show_bug.cgi?id=43969
    if(allocated(lhs%Tensors)) deallocate(lhs%Tensors)
    
    allocate (lhs%Tensors(rhs%Length))

    do n=1,rhs%Length
       lhs%Tensors(n) = rhs%Tensors(n)
    enddo

    lhs%initialized=.true.

  end subroutine new_MatrixProductState_from_Assignment


  logical function is_full_MPS_Canonized(this)
    class(MatrixProductState),intent(in) :: this

     if(.not.this%initialized) then
        call ThrowException('is_full_MPS_Canonized','MPS not initialized',NoErrorCode,CriticalError)
        return
     endif
    if(this%CanonizedAtSite.le.-1) then 
       is_full_MPS_Canonized=.false.
    else
       is_full_MPS_Canonized=.true.
    endif
    return
  end function is_full_MPS_Canonized


  real(8) function Canonize_MPS_At_Site(this,RequestedSite)
    class(MatrixProductState),intent(inout) :: this
    type(MPSTensor) :: InnerMatrix
    integer,optional :: RequestedSite
    integer n,CanonizationSite

     if(.not.this%initialized) then
        call ThrowException('Canonize_MPS_At_Site','MPS not initialized',NoErrorCode,CriticalError)
        return
     endif
     if(RequestedSite.gt.(this%Length)+1.or.(RequestedSite.lt.0)) then
        call ThrowException('Canonize_MPS_At_Site','Canonization Site not correct',RequestedSite,MinorError)
        return
     endif

    if(present(RequestedSite)) then
       CanonizationSite=RequestedSite
    else
       CanonizationSite=0
    endif

    do n=this%Length,max(CanonizationSite+1,2),-1
       InnerMatrix=this%Tensors(n)%RCanonize()
       this%Tensors(n-1)=this%Tensors(n-1)*InnerMatrix
    enddo

    do n=1,min(CanonizationSite-1,this%Length-1),1
       InnerMatrix=this%Tensors(n)%LCanonize()
       this%Tensors(n+1)=InnerMatrix*this%Tensors(n+1)
    enddo

    if (.not.WasThereError()) then
       if(CanonizationSite.eq.0) then
          InnerMatrix=this%Tensors(1)%RCanonize()
          Canonize_MPS_At_Site=InnerMatrix%Norm()
       else if (CanonizationSite.eq.this%Length+1) then
          InnerMatrix=this%Tensors(this%Length)%LCanonize()
          Canonize_MPS_At_Site=InnerMatrix%Norm()
       else
          Canonize_MPS_At_Site=CanonizationSite
       endif
    else
       call ThrowException('new_MatrixProductState_from_Assignment','Original MPS not initialized',NoErrorCode,CriticalError)
    endif

    this%CanonizedAtSite=CanonizationSite
    return

  end function Canonize_MPS_At_Site

  real(8) function Canonize_MPS_At_Right(this)
    class(MatrixProductState),intent(inout) :: this
    
    Canonize_MPS_At_Right=Canonize_MPS_At_Site(this,this%Length+1)
    return
    
  end function Canonize_MPS_At_Right
  
  real(8) function Canonize_MPS_At_Left(this)
    class(MatrixProductState),intent(inout) :: this
    
    Canonize_MPS_At_Left=Canonize_MPS_At_Site(this,0)
    return
    
  end function Canonize_MPS_At_Left

  real(8) function  Norm_of_MatrixProductState(this)
    class(MatrixProductState),intent(inout) :: this
    integer :: n

    Norm_of_MatrixProductState=zero
    do n=1,this%Length
       Norm_of_MatrixProductState=Norm_of_MatrixProductState+this%Tensors(n)%Norm()
    enddo
    
  end function Norm_of_MatrixProductState

end module MatrixProductState_Class
