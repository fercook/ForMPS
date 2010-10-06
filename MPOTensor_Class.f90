!!   Copyright 2010 Fernando M. Cucchietti
!
!    This file is part of ForMPS
!
!    ForMPS is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    ForMPS is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with ForMPS.  If not, see <http://www.gnu.org/licenses/>.

module MPOTensor_Class

    use ErrorHandling
    use Constants
    use Tensor_Class
    use Operator_Class
    use MPSTensor_Class

    implicit none
!    private

  type,public,extends(Tensor4) :: MPOTensor
     private
     integer :: spin,DLeft,DRight
   contains
     procedure,public :: getDRight => DRight_MPOTensor
     procedure,public :: getDLeft => DLeft_MPOTensor
     procedure,public :: getMaxBondDimension => Get_MaxBondDimensionMPOTensor
     procedure,public :: getSpin => Spin_MPOTensor
     procedure,public :: PrintDimensions => Print_MPOTensor_Dimensions
  end type MPOTensor

  interface new_MPOTensor
     module procedure new_MPOTensor_Random,new_MPOTensor_fromMPOTensor, &
          & new_MPOTensor_fromSplitData
  end interface

  interface assignment (=)
     module procedure new_MPOTensor_fromAssignment
  end interface

  interface operator (.applyTo.)
     module procedure Apply_MPO_to_MPS_Tensor
  end interface

contains

   function new_MPOTensor_Random (spin,DLeft,DRight) result (this)
     integer,intent(in) :: spin,DLeft,DRight
     type(MPOTensor) :: this

     this=new_Tensor(DLeft,DRight,spin,spin)
     !initialize internal variables
     this%spin=spin
     this%DLeft=DLeft
     this%DRight=DRight

   end function new_MPOTensor_Random

!##################################################################
   function new_MPOTensor_fromSplitData (originalData) result (this)
     complex(8),intent(in) :: originalData(:,:,:,:)
     complex(8),allocatable :: JoinedData(:,:,:,:)
     type(MPOTensor) this
     integer :: dims(4),spin,DLeft,Dright,n,m

     Dims=shape(originalData)

     if (dims(1).eq.dims(2)) then

        allocate (JoinedData(Dims(3),Dims(4),Dims(1),Dims(2)))
        do n=1,dims(1)
            do m=1,dims(2)
                JoinedData(:,:,n,m)=originalData(n,m,:,:)
            enddo
        enddo

        this=new_Tensor(JoinedData)
        !initialize internal variables
        this%spin=dims(1)
        this%DLeft=dims(3)
        this%DRight=dims(4)

      else
        call ThrowException('new_MPOTensor_with_SplitData','spin dimensions have different size',Dims(1)-Dims(2),CriticalError)
      endif

   end function new_MPOTensor_fromSplitData

!##################################################################
   function new_MPOTensor_fromMPOTensor (tensor) result (this)
     class(MPOTensor),intent(in) :: tensor
     type(MPOTensor) this

     this=new_Tensor(tensor)
     !initialize internal variables
     this%spin=tensor%spin
     this%DLeft=tensor%DLeft
     this%DRight=tensor%DRight

   end function new_MPOTensor_fromMPOTensor

!##################################################################
   subroutine new_MPOTensor_fromAssignment(lhs,rhs)
     class(MPOTensor),intent(out) :: lhs
     type(MPOTensor),intent(in) :: rhs

     lhs=new_Tensor(rhs)

     lhs%spin=rhs%spin
     lhs%DLeft=rhs%DLeft
     lhs%DRight=rhs%DRight

   end subroutine new_MPOTensor_fromAssignment

!##################################################################
!###########       Accessor methods
!##################################################################
   integer function spin_MPOTensor(this) result(s)
     class(MPOTensor),intent(IN) :: this  !!<<TYPE>>!!

    if(.not.(this%IsInitialized())) then
        call ThrowException('Spin','Tensor not initialized',NoErrorCode,Warning)
        return
     else
        s=this%spin
     endif

   end function spin_MPOTensor
!##################################################################

   integer function DLeft_MPOTensor(this) result(DL)
     class(MPOTensor),intent(IN) :: this   !!<<TYPE>>!!

     if(.not.(this%IsInitialized())) then
        call ThrowException('DLeft','Tensor not initialized',NoErrorCode,Warning)
        return
     else
        DL=this%DLeft
     endif

   end function DLeft_MPOTensor
!##################################################################
   integer function DRight_MPOTensor(this) result(DR)
     class(MPOTensor),intent(IN) :: this   !!<<TYPE>>!!

     if(.not.(this%IsInitialized())) then
        call ThrowException('DRight','Tensor not initialized',NoErrorCode,Warning)
        return
     else
        DR=this%DRight
     endif

   end function DRight_MPOTensor

   integer function Get_MaxBondDimensionMPOTensor(this) result(maxDimension)
     class(MPOTensor),intent(IN) :: this   !!<<TYPE>>!!

     if(.not.(this%IsInitialized())) then
        call ThrowException('Get_MaxBondDimensionMPOTensor','Tensor not initialized',NoErrorCode,Warning)
        return
     else
        maxDimension=max(this%DRight,this%DLeft)
     endif

   end function Get_MaxBondDimensionMPOTensor

    subroutine Print_MPOTensor_Dimensions(this,message)
        class(MPOTensor) :: this
        character*(*),optional :: message
        integer :: dims(4)

        if((this%IsInitialized())) then
            dims=this%GetDimensions()
            if (present(message) ) then
                write(*,'(A,", Spin: ",I3,", DLeft: ",I3,", DRight: ",I3)') ,message,dims(3),dims(1),dims(2)
            else
                write(*,'("Spin: ",I3,", DLeft: ",I3,", DRight: ",I3)') ,dims(3),dims(1),dims(2)
            endif
        else
            call ThrowException('Print_MPOTensor_Dimensions','Tensor not initialized',NoErrorCode,Warning)
            return
         endif

   end subroutine Print_MPOTensor_Dimensions

!##################################################################
!##################################################################
!##################################################################
!##################################################################

    function Apply_MPO_to_MPS_Tensor(anMPO,anMPS) result(this)
        class(MPOTensor),intent(IN) :: anMPO
        class(MPSTensor),intent(IN) :: anMPS
        type(MPSTensor) :: this

        this=new_MPSTensor(CompactBelow(anMPS,THIRD,anMPO,THIRD,FOURTH))

    end function Apply_MPO_to_MPS_Tensor

end module MPOTensor_Class
