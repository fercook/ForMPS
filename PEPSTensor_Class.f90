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

module PEPSTensor_Class

  use ErrorHandling
  use Constants
  use Tensor_Class
  use Operator_Class

  implicit none

!  private

!###############################
!#####  The class main object
!###############################
  type,public,extends(Tensor5) :: PEPSTensor
     private
   contains
     procedure,public :: getDRight => DRight_PEPSTensor
     procedure,public :: getDLeft => DLeft_PEPSTensor
     procedure,public :: getDUp => DUp_PEPSTensor
     procedure,public :: getDDown => DDown_PEPSTensor
     procedure,public :: getMaxBondDimension => Get_MaxBondDimensionPEPSTensor
     procedure,public :: getSpin => Spin_PEPSTensor
     procedure,public :: PrintDimensions => Print_PEPSTensor_Dimensions
     procedure,public :: ApplyOperator => Apply_Operator_To_PEPS_Spin_Dimension
     procedure,public :: Collapse => Collapse_PEPS_Into_Tensor4
     procedure,public :: CompactBonds => CompactPEPSBondDimensions
  end type PEPSTensor

!###############################
!#####  Operators and methods
!###############################
  interface new_PEPSTensor
     module procedure new_PEPSTensor_Random,new_PEPSTensor_fromPEPSTensor, &
          & new_PEPSTensor_withConstant, new_PEPSTensor_with_SplitData, &
          & new_PEPSTensor_fromTensor5, new_PEPSTensor_fromTensor5_Transposed
  end interface

  interface assignment (=)
     module procedure new_PEPSTensor_fromAssignment
  end interface

!  interface operator (.times.)
!     module procedure PEPSTensor_times_matrix, matrix_times_PEPSTensor
!  end interface

  interface operator (.apply.)
     module procedure Apply_Operator_To_PEPS_Spin_Dimension
  end interface

  interface CollapsePEPS
     module procedure Collapse_PEPS_Into_Tensor4, Collapse_Two_PEPS_into_Tensor4
  end interface

!######################################################################################
!######################################################################################
!######################################################################################
!######################################################################################
!######################################################################################
!######################################################################################
!######################################################################################

 contains

!######################################################################################
!#####                           Creation operators
!######################################################################################
   function new_PEPSTensor_Random (spin,DLeft,DRight,Dup,DDown) result (this)
     integer,intent(in) :: spin,DLeft,DRight,Dup,DDown
     type(PEPSTensor) :: this

     this=new_Tensor(DLeft,DRight,Dup,DDown,spin)

   end function new_PEPSTensor_Random

!##################################################################
   function new_PEPSTensor_with_SplitData (originalDataUP,originalDataDOWN) result (this)
     complex(8),intent(in) :: originalDataUP(:,:,:,:),originalDataDOWN(:,:,:,:)
     complex(8),allocatable :: JoinedData(:,:,:,:,:)
     type(PEPSTensor) this
     integer :: upDims(4),downDims(4)

     upDims=shape(originalDataUP)
     downDims=shape(originalDataDOWN)

     if (upDims.equalvector.downDims) then
        allocate (JoinedData(upDims(1),upDims(2),upDims(3),upDims(4),2))
        JoinedData(:,:,:,:,1)=originalDataUP
        JoinedData(:,:,:,:,2)=originalDataDOWN

        this=new_Tensor(JoinedData)

      else
        call ThrowException('new_PEPSTensor_with_SplitData','Up and down data have different size',upDims(1)-downDims(2),CriticalError)
      endif
   end function new_PEPSTensor_with_SplitData

!##################################################################
   function new_PEPSTensor_withConstant (spin,DLeft,DRight,DUp,DDown,constant) result (this)
     integer,intent(in) :: spin,DLeft,DRight,DUp,DDown
     complex(8),intent(in) :: constant
     type(PEPSTensor) :: this

     this=new_Tensor(DLeft,DRight,DUp,DDown,spin,constant)

   end function new_PEPSTensor_withConstant

!##################################################################
   function new_PEPSTensor_fromPEPSTensor (tensor) result (this)
     class(PEPSTensor),intent(in) :: tensor
     type(PEPSTensor) this

     this=new_Tensor(tensor)

   end function new_PEPSTensor_fromPEPSTensor

!##################################################################
   subroutine new_PEPSTensor_fromAssignment(lhs,rhs)
     class(PEPSTensor),intent(out) :: lhs
     type(PEPSTensor),intent(in) :: rhs

     lhs=new_Tensor(rhs)

   end subroutine new_PEPSTensor_fromAssignment

!##################################################################
   function new_PEPSTensor_fromTensor5 (tensor) result (this)
     type(Tensor5),intent(in) :: tensor
     type(PEPSTensor) this
     integer :: dims(5)

     this=new_Tensor(tensor)
     dims=tensor%GetDimensions()

   end function new_PEPSTensor_fromTensor5

   function new_PEPSTensor_fromTensor5_Transposed(tensor, &
       &  whichDimIsSpin,whichDimIsLeft,whichDimIsRight,whichDimIsUp,whichDimIsDown ) result (this)
     type(Tensor5),intent(in) :: tensor
     integer, intent(IN) :: whichDimIsSpin,whichDimIsLeft,whichDimIsRight,whichDimIsUp,whichDimIsDown
     type(PEPSTensor) this
     integer :: newDims(5)

     newDims(whichDimIsLeft)=1
     newDims(whichDimIsRight)=2
     newDims(whichDimIsUp)=3
     newDims(whichDimIsDown)=4
     newDims(whichDimIsSpin)=5
     this=TensorTranspose(tensor,newDims)

   end function new_PEPSTensor_fromTensor5_Transposed



!##################################################################
!###########       Accessor methods
!##################################################################
   integer function spin_PEPSTensor(this) result(s)
     class(PEPSTensor),intent(IN) :: this
     integer :: dims(5)

    if(.not.(this%IsInitialized())) then
        call ThrowException('spin_PEPSTensor','Tensor not initialized',NoErrorCode,Warning)
        return
     else
        dims=this%GetDimensions()
        s=dims(5)
     endif

   end function spin_PEPSTensor
!##################################################################

   integer function DLeft_PEPSTensor(this) result(DL)
     class(PEPSTensor),intent(IN) :: this
     integer :: dims(5)

     if(.not.(this%IsInitialized())) then
        call ThrowException('DLeft_PEPS','Tensor not initialized',NoErrorCode,Warning)
        return
     else
        dims=this%GetDimensions()
        DL=dims(1)
     endif

   end function DLeft_PEPSTensor
!##################################################################
   integer function DRight_PEPSTensor(this) result(DR)
     class(PEPSTensor),intent(IN) :: this
     integer :: dims(5)

     if(.not.(this%IsInitialized())) then
        call ThrowException('DRight_PEPS','Tensor not initialized',NoErrorCode,Warning)
        return
     else
        dims=this%GetDimensions()
        DR=dims(2)
     endif

   end function DRight_PEPSTensor

!##################################################################

   integer function DUp_PEPSTensor(this) result(DU)
     class(PEPSTensor),intent(IN) :: this
     integer :: dims(5)

     if(.not.(this%IsInitialized())) then
        call ThrowException('DUp_PEPS','Tensor not initialized',NoErrorCode,Warning)
        return
     else
        dims=this%GetDimensions()
        DU=dims(3)
     endif

   end function DUp_PEPSTensor

!##################################################################
   integer function DDown_PEPSTensor(this) result(DD)
     class(PEPSTensor),intent(IN) :: this
     integer :: dims(5)

     if(.not.(this%IsInitialized())) then
        call ThrowException('DDown_PEPS','Tensor not initialized',NoErrorCode,Warning)
        return
     else
        dims=this%GetDimensions()
        DD=dims(4)
     endif

   end function DDown_PEPSTensor
!##################################################################

   integer function Get_MaxBondDimensionPEPSTensor(this) result(maxDimension)
     class(PEPSTensor),intent(IN) :: this
     integer :: dims(5)

     if(.not.(this%IsInitialized())) then
        call ThrowException('Get_MaxBondDimensionPEPSTensor','Tensor not initialized',NoErrorCode,Warning)
        return
     else
        dims=this%GetDimensions()
        maxDimension=maxval(dims)
     endif

   end function Get_MaxBondDimensionPEPSTensor

       subroutine Print_PEPSTensor_Dimensions(this,message)
        class(PEPSTensor) :: this
        character*(*),optional :: message
        integer :: dims(5)

        if((this%IsInitialized())) then
            if (present(message) ) write(*,'(A)') message
            dims=this%GetDimensions()
            write(*,'("Spin: ",I3)') ,dims(5)
            write(*,'("DLeft: ",I3,", DRight: ",I3)') dims(1),dims(2)
            write(*,'("DUp: ",I3,", DDown: ",I3)') dims(3),dims(4)
        else
            call ThrowException('Print_PEPSDims','Tensor not initialized',NoErrorCode,Warning)
            return
        endif

   end subroutine Print_PEPSTensor_Dimensions


!##################################################################
!#######################        Products by things
!##################################################################

   function Apply_Operator_To_PEPS_Spin_Dimension(this,anOperator) result(aTensor)
     class(PEPSTensor),intent(IN) :: this
     type(PEPSTensor) :: aTensor
     class(SpinOperator),intent(IN) :: anOperator
     integer :: opDims(2),tensorDims(5)

     opDims=anOperator%GetDimensions()
     if(this%IsInitialized()) then
        if(opDims(1).eq.this%getSpin().and.opDims(2).eq.this%getSpin()) then
           aTensor = new_Tensor(this*anOperator)
        else
           call ThrowException('Apply_Operator_To_Spin_Dimension','Operator is not of the right size',opDims(1)-opDims(2),CriticalError)
        endif
     else
        call ThrowException('Apply_Operator_To_Spin_Dimension','Tensor not initialized',NoErrorCode,CriticalError)
     endif

   end function Apply_Operator_To_PEPS_Spin_Dimension


!##################################################################
!##################################################################
!##################################################################

    function Collapse_PEPS_Into_Tensor4(this,anOperator) result (aTensor)
        class(PEPSTensor),intent(IN) :: this
        class(SpinOperator),intent(IN),optional :: anOperator
        type(Tensor4) :: aTensor

        if(this%IsInitialized()) then
            if(present(anOperator)) then
                aTensor= MirrorCompact(this, this.apply.anOperator, FIFTH)
            else
                aTensor= MirrorCompact(this, this, FIFTH)
            endif
        else
            call ThrowException('Collapse_PEPS_Into_Tensor4','Tensor not initialized',NoErrorCode,CriticalError)
        endif

    end function Collapse_PEPS_Into_Tensor4

!##################################################################

    function Collapse_Two_PEPS_into_Tensor4(upPEPS,downPEPS,anOperator) result (aTensor)
        class(PEPSTensor),intent(IN) :: upPEPS,downPEPS
        class(SpinOperator),intent(IN),optional :: anOperator
        type(Tensor4) :: aTensor

        if(upPEPS%IsInitialized().and.downPEPS%IsInitialized()) then
            if(present(anOperator)) then
                aTensor= MirrorCompact(upPEPS.apply.anOperator, downPEPS, FIFTH)
            else
                aTensor= MirrorCompact(upPEPS, downPEPS, FIFTH)
            endif
        else
            call ThrowException('Collapse_Two_PEPS_into_Tensor4','Tensor not initialized',NoErrorCode,CriticalError)
        endif

    end function Collapse_Two_PEPS_into_Tensor4

!##################################################################

    function CompactPEPSBondDimensions(aPEPS) result(aMatrix)
        class(PEPSTensor),intent(IN) :: aPEPS
        type(Tensor2) :: aMatrix

        if(aPEPS%Isinitialized()) then
            aMatrix=aPEPS%JoinIndices()
        else
            call ThrowException('CompactPEPSBondDimensions','Tensor not initialized',NoErrorCode,CriticalError)
        endif

    end function CompactPEPSBondDimensions

!##################################################################


end module PEPSTensor_Class
