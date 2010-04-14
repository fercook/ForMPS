!!  This module describes an object with three legs as used by MPS algorithms.
!!  One leg is "special", it is the physical leg of the tensor.

module Tensor

  use ErrorHandling

  implicit none

  !Set maximum bond dimensions and spin
  integer,parameter :: MAX_spin = 2, MAX_D = 100

  type MPSTensor
     private
     integer spin,DLeft,DRight !D is bond dimension
     logical :: initialized=.false.
     complex(8) data(MAX_D,MAX_D,MAX_spin) !Now uses fixed max dimensions and internal variables to adjust
  end type MPSTensor

  interface operator (.print.)
     module procedure PrintMPSTensor
  end interface

 contains

   
!###########  new
   function new_MPSTensor (s,DL,DR) result (this)
     integer,intent(in) :: s,DL,DR
     type(MPSTensor) this

     if(s.gt.MAX_spin.or.DL.gt.MAX_D.or.DR.gt.MAX_D) then
        call ThrowException('new_MPSTensor','spin or bond dimension larger than maximum',NoErrorCode,CriticalError)
        return
     endif
     if(s.lt.1.or.DL.lt.1.or.DR.lt.1) then
        call ThrowException('new_MPSTensor','spin or bond dimension smaller than 1',NoErrorCode,CriticalError)
        return
     endif

     !initialize internal variables
     this%spin=s
     this%DLeft=DL
     this%DRight=DR
     !initialize data
     this%data=0.0d0
     !flip flag
     this%initialized=.true.

   end function new_MPSTensor
!###########


!########### delete
   integer function delete_MPSTensor (this) result(error)
     type(MPSTensor),intent(INOUT) :: this

     error=Warning

     if(.not.this%initialized) then
        call ThrowException('delete_MPSTensor','Trying to delete an uninitialized tensor',NoErrorCode,Warning)
        error = Warning
        return
     endif
     
     !Erase info
     this%spin=0
     this%DLeft=0
     this%DRight=0
     !Erase data
     this%data=0.0d0
     !Flip flag
     this%initialized=.false.     

     error=Normal

   end function delete_MPSTensor
!###########

!###########  print
   integer function PrintMPSTensor(this) result(error)
     type(MPSTensor),intent(IN) :: this
     integer i,j,k

     error = Warning
     
     if(.not.(this%initialized)) then
        call ThrowException('PrintMPSTensor','Tensor not initialized',NoErrorCode,Warning)
        return
     endif

     do i=1,this%spin
        print *,'State :',i
        do j=1,this%DLeft
           print *,(this%data(j,k,i),k=1,this%DRight)
        enddo
     enddo
     error=Normal
   end function PrintMPSTensor
!###########       

!###########       Accessor methods


end module Tensor
