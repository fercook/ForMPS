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

module PEPOTensor_Class

    use ErrorHandling
    use Constants
    use Tensor_Class
    use Operator_Class
    use PEPSTensor_Class

    implicit none
!    private

  type,public,extends(Tensor6) :: PEPOTensor
     private
     integer :: spin,DLeft,DRight,DUp,DDown
   contains
     procedure,public :: getDRight => DRight_PEPOTensor
     procedure,public :: getDLeft => DLeft_PEPOTensor
     procedure,public :: getDUp => DUp_PEPOTensor
     procedure,public :: getDDown => DDown_PEPOTensor
     procedure,public :: getMaxBondDimension => Get_MaxBondDimensionPEPOTensor
     procedure,public :: getSpin => Spin_PEPOTensor
     procedure,public :: PrintDimensions => Print_PEPOTensor_Dimensions
  end type PEPOTensor

  interface new_PEPOTensor
     module procedure new_PEPOTensor_Random,new_PEPOTensor_fromPEPOTensor, &
          & new_PEPOTensor_fromSplitData, new_PEPOTensor_fromTensor6_Transposed, &
          & new_PEPOTensor_withConstant
  end interface

  interface assignment (=)
     module procedure new_PEPOTensor_fromAssignment
  end interface

  interface operator (.applyTo.)
     module procedure Apply_PEPO_to_PEPS_Tensor
  end interface

contains

   function new_PEPOTensor_Random (spin,DLeft,DRight,DUp,DDown) result (this)
     integer,intent(in) :: spin,DLeft,DRight,DUp,DDown
     type(PEPOTensor) :: this

     this=new_Tensor(DLeft,DRight,DUp,DDown,spin,spin)
     !initialize internal variables
     this%spin=spin
     this%DLeft=DLeft
     this%DRight=DRight
     this%DUp=DUp
     this%DDown=DDown

   end function new_PEPOTensor_Random


   function new_PEPOTensor_withConstant(spin,DLeft,DRight,DUp,DDown,Const) result (this)
     integer,intent(in) :: spin,DLeft,DRight,DUp,DDown
     complex(8) :: const
     type(PEPOTensor) :: this

     this=new_Tensor(DLeft,DRight,DUp,DDown,spin,spin,Const)
     !initialize internal variables
     this%spin=spin
     this%DLeft=DLeft
     this%DRight=DRight
     this%DUp=DUp
     this%DDown=DDown

   end function new_PEPOTensor_withConstant


!##################################################################
   function new_PEPOTensor_fromSplitData (originalData) result (this)
     complex(8),intent(in) :: originalData(:,:,:,:,:,:)
     complex(8),allocatable :: JoinedData(:,:,:,:,:,:)
     type(PEPOTensor) this
     integer :: dims(6),spin,DLeft,Dright,DUp,DDown,n,m

     Dims=shape(originalData)

     if (dims(1).eq.dims(2)) then

        allocate (JoinedData(Dims(3),Dims(4),Dims(5),Dims(6),Dims(1),Dims(2)))
        do n=1,dims(1)
            do m=1,dims(2)
                JoinedData(:,:,:,:,n,m)=originalData(n,m,:,:,:,:)
            enddo
        enddo

        this=new_Tensor(JoinedData)
        !initialize internal variables
        this%spin=dims(1)
        this%DLeft=dims(3)
        this%DRight=dims(4)
        this%DUp=dims(5)
        this%DDown=dims(6)

      else
        call ThrowException('new_PEPOTensor_with_SplitData','spin dimensions have different size',Dims(1)-Dims(2),CriticalError)
      endif

   end function new_PEPOTensor_fromSplitData

!##################################################################
   function new_PEPOTensor_fromPEPOTensor (tensor) result (this)
     class(PEPOTensor),intent(in) :: tensor
     type(PEPOTensor) this

     this=new_Tensor(tensor)
     !initialize internal variables
     this%spin=tensor%spin
     this%DLeft=tensor%DLeft
     this%DRight=tensor%DRight
     this%DUp=tensor%DUp
     this%DDown=tensor%DDown

   end function new_PEPOTensor_fromPEPOTensor

!##################################################################
   subroutine new_PEPOTensor_fromAssignment(lhs,rhs)
     class(PEPOTensor),intent(out) :: lhs
     type(PEPOTensor),intent(in) :: rhs

     lhs=new_Tensor(rhs)

     lhs%spin=rhs%spin
     lhs%DLeft=rhs%DLeft
     lhs%DRight=rhs%DRight
     lhs%DUp=rhs%DUp
     lhs%DDown=rhs%DDown

   end subroutine new_PEPOTensor_fromAssignment

!##################################################################
   function new_PEPOTensor_fromTensor6_Transposed (tensor,whichDimIsSpinUp,whichDimIsSpinDown, &
          & whichDimIsLeft,whichDimIsRight,whichDimIsUp,whichDimIsDown) result (this)
     type(Tensor6),intent(in) :: tensor
     integer, intent(IN) :: whichDimIsSpinUp,whichDimIsSpinDown,whichDimIsLeft,whichDimIsRight,whichDimIsUp,whichDimIsDown
     type(PEPOTensor) this
     integer :: newDims(6)

     newDims=tensor%GetDimensions()
     if(newDims(whichDimIsSpinUp).eq.newDims(whichDimIsSpinDown)) then
        newDims(whichDimIsLeft)=1
        newDims(whichDimIsRight)=2
        newDims(whichDimIsUp)=3
        newDims(whichDimIsDown)=4
        newDims(whichDimIsSpinUp)=5
        newDims(whichDimIsSpinDown)=6
        this=TensorTranspose(tensor,newDims)
        !initialize internal variables
        newDims=this%GetDimensions()
        this%spin=newDims(5)
        this%DLeft=newDims(1)
        this%DRight=newDims(2)
        this%DUp=newDims(3)
        this%DDown=newDims(4)
     else
        call ThrowException('PEPOTensor from Tensor6','spin dimensions are not equal',NoErrorCode,Warning)
        return
     endif
   end function new_PEPOTensor_fromTensor6_Transposed

!##################################################################
!###########       Accessor methods
!##################################################################
   integer function spin_PEPOTensor(this) result(s)
     class(PEPOTensor),intent(IN) :: this  !!<<TYPE>>!!

    if(.not.(this%IsInitialized())) then
        call ThrowException('Spin','Tensor not initialized',NoErrorCode,Warning)
        return
     else
        s=this%spin
     endif

   end function spin_PEPOTensor
!##################################################################

   integer function DLeft_PEPOTensor(this) result(DL)
     class(PEPOTensor),intent(IN) :: this   !!<<TYPE>>!!

     if(.not.(this%IsInitialized())) then
        call ThrowException('DLeft_PEPO','Tensor not initialized',NoErrorCode,Warning)
        return
     else
        DL=this%DLeft
     endif

   end function DLeft_PEPOTensor
!##################################################################
   integer function DUp_PEPOTensor(this) result(DU)
     class(PEPOTensor),intent(IN) :: this   !!<<TYPE>>!!

     if(.not.(this%IsInitialized())) then
        call ThrowException('DRight','Tensor not initialized',NoErrorCode,Warning)
        return
     else
        DU=this%DUp
     endif

   end function DUp_PEPOTensor
!##################################################################
   integer function DDown_PEPOTensor(this) result(DD)
     class(PEPOTensor),intent(IN) :: this   !!<<TYPE>>!!

     if(.not.(this%IsInitialized())) then
        call ThrowException('DRight','Tensor not initialized',NoErrorCode,Warning)
        return
     else
        DD=this%DDown
     endif

   end function DDown_PEPOTensor
!##################################################################
   integer function DRight_PEPOTensor(this) result(DR)
     class(PEPOTensor),intent(IN) :: this   !!<<TYPE>>!!

     if(.not.(this%IsInitialized())) then
        call ThrowException('DRight','Tensor not initialized',NoErrorCode,Warning)
        return
     else
        DR=this%DRight
     endif

   end function DRight_PEPOTensor
!##################################################################

   integer function Get_MaxBondDimensionPEPOTensor(this) result(maxDimension)
     class(PEPOTensor),intent(IN) :: this   !!<<TYPE>>!!

     if(.not.(this%IsInitialized())) then
        call ThrowException('Get_MaxBondDimensionPEPOTensor','Tensor not initialized',NoErrorCode,Warning)
        return
     else
        maxDimension=max(this%DRight,this%DLeft,this%DUp,this%DDown)
     endif

   end function Get_MaxBondDimensionPEPOTensor

    subroutine Print_PEPOTensor_Dimensions(this,message)
        class(PEPOTensor) :: this
        character*(*),optional :: message
        integer :: dims(6)

        if((this%IsInitialized())) then
            if (present(message) ) write(*,'(A)') message
            dims=this%GetDimensions()
            write(*,'("Spin: ",I3)') dims(5)
            write(*,'("DLeft: ",I3,", DRight: ",I3)') dims(1),dims(2)
            write(*,'("DUp: ",I3,", DDown: ",I3)') dims(3),dims(4)
        else
            call ThrowException('Print_PEPOTensor_Dimensions','Tensor not initialized',NoErrorCode,Warning)
            return
         endif

   end subroutine Print_PEPOTensor_Dimensions

!##################################################################
!##################################################################
!##################################################################
!##################################################################

    function Apply_PEPO_to_PEPS_Tensor(aPEPO,aPEPS) result(this)
        class(PEPOTensor),intent(IN) :: aPEPO
        class(PEPSTensor),intent(IN) :: aPEPS
        type(PEPSTensor) :: this

        this=new_PEPSTensor(CompactBelow(aPEPS,FIFTH,aPEPO,FIFTH,SIXTH))

    end function Apply_PEPO_to_PEPS_Tensor

end module PEPOTensor_Class
